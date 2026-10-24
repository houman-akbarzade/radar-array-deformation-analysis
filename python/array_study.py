"""Public C-band cylindrical-array deformation study (NumPy/SciPy port)."""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass, replace
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize_scalar


@dataclass(frozen=True)
class Config:
    c0: float = 299_792_458.0
    f_low: float = 5.2e9
    f_high: float = 5.8e9
    frequency: float = 5.8e9
    total_columns: int = 60
    active_columns: int = 20
    rows: int = 16
    seed: int = 3719
    trials: int = 500
    steer_az: float = 0.0
    analog_el: float = 0.0
    target_az: float = 1.0
    target_el: float = 4.0
    beam_el: float = 4.0
    calibration_el: float = 4.0
    rho_panel: float = .35
    rho_within: float = .995
    radial_ratio: float = 1e-4
    vertical_ratio: float = .02
    mean_tip_wavelengths: float = 0.0
    sigma_tip_wavelengths: float = .02
    phase_center_power: float = 1.0
    bending_shape: str = "uniform"
    mechanical_model: str = "base"
    wind_az: float = 0.0
    height_gradient: float = 0.0
    element_pattern: str = "isotropic"

    @property
    def reference_wavelength(self):
        return self.c0 / self.f_high

    @property
    def pitch(self):
        return 360.0 / self.total_columns

    @property
    def spacing(self):
        return .5 * self.reference_wavelength

    @property
    def radius(self):
        return self.spacing / (2 * np.sin(np.pi / self.total_columns))


def geometry(c: Config):
    if c.active_columns % 2 or c.active_columns > c.total_columns:
        raise ValueError("A centred, even active sector is required")
    az = (np.arange(c.active_columns) - (c.active_columns - 1) / 2) * c.pitch
    phi = np.tile(np.deg2rad(az), c.rows)
    z = np.repeat((np.arange(c.rows) - (c.rows - 1) / 2) * c.spacing, c.active_columns)
    panel = np.tile(np.arange(c.active_columns), c.rows)
    row = np.repeat(np.arange(c.rows), c.active_columns)
    xyz = np.column_stack((c.radius * np.cos(phi), c.radius * np.sin(phi), z))
    return dict(phi=phi, z=z, panel=panel, row=row, xyz=xyz, az=az)


def panel_correlation(c, g):
    # Covariance follows an open active arc, without first/last neighbour wrap.
    arc = c.radius * np.deg2rad(g["az"])
    d = np.abs(arc[:, None] - arc[None, :])
    if c.rho_panel == 0:
        return np.eye(c.active_columns)
    if c.rho_panel == 1:
        return np.ones((c.active_columns, c.active_columns))
    length = np.median(np.diff(arc)) / -np.log(c.rho_panel)
    return np.exp(-d / length)


def deformation(c, g, rng, panel_root=None, residual=None):
    if not (c.f_low <= c.frequency <= c.f_high):
        raise ValueError("Frequency outside the selected public band")
    xi = ((c.frequency - c.f_low) / (c.f_high - c.f_low)) ** c.phase_center_power
    if c.bending_shape == "uniform":
        gain = xi**2 * (6 - 4*xi + xi**2) / 3
    elif c.bending_shape == "tip_force":
        gain = xi**2 * (3 - xi) / 2
    else:
        raise ValueError("Unknown bending shape")
    k = panel_correlation(c, g)
    if panel_root is None:
        val, vec = np.linalg.eigh((k + k.T) / 2)
        if val.min() < -1e-10:
            raise ValueError("Invalid covariance")
        panel_root = vec @ np.diag(np.sqrt(np.maximum(val, 0)))
    p = panel_root @ rng.standard_normal(c.active_columns)
    if residual is None:
        residual = rng.standard_normal(g["phi"].size)
    q = np.sqrt(c.rho_within)*p[g["panel"]] + np.sqrt(1-c.rho_within)*residual
    scale = np.ones_like(q)
    if c.mechanical_model == "wind_parametric":
        scale = np.sin(np.deg2rad(c.wind_az)-g["phi"]) * (1+c.height_gradient*g["z"]/np.max(np.abs(g["z"])))
    elif c.mechanical_model != "base":
        raise ValueError("Unknown mechanical model")
    s = gain*scale*c.reference_wavelength*(c.mean_tip_wavelengths+c.sigma_tip_wavelengths*q)
    r = c.radial_ratio*s
    dz = c.vertical_ratio*s
    delta = np.column_stack((-np.sin(g["phi"])*s+np.cos(g["phi"])*r,
                             np.cos(g["phi"])*s+np.sin(g["phi"])*r, dz))
    actual = g["xyz"] + delta
    radius_change = np.hypot(actual[:, 0], actual[:, 1])-c.radius
    angle_change = np.angle(np.exp(1j*(np.arctan2(actual[:, 1],actual[:, 0])-g["phi"])))
    return dict(tangential=s, radial=r, vertical=dz, xyz=delta,
                radius_change=radius_change, angle_change=angle_change, gain=gain)


def taylor_window(n, nbar=5, sll=30):
    a = np.arccosh(10**(sll/20))/np.pi
    s2 = nbar**2/(a*a+(nbar-.5)**2)
    m = np.arange(1, nbar)
    factors = []
    for j in m:
        other = m[m != j]
        factors.append((-1)**(j+1)*np.prod(1-j*j/(s2*(a*a+(m-.5)**2))) /
                       (2*np.prod(1-j*j/other**2)))
    nvec = np.arange(n)
    v = 1+2*np.cos(2*np.pi*(nvec[:, None]-n/2+.5)*m/n) @ np.asarray(factors)
    return v/(1+2*np.sum(factors))


def weights(c, g, architecture):
    k = 2*np.pi*c.frequency/c.c0
    azwin = np.hamming(c.active_columns)
    elwin = taylor_window(c.rows)
    analog_el = c.beam_el if architecture == "A" else c.analog_el
    xyz = g["xyz"]
    horizontal = xyz[:, 0]*np.cos(np.deg2rad(c.steer_az))+xyz[:, 1]*np.sin(np.deg2rad(c.steer_az))
    combiner = azwin[g["panel"]]*np.exp(-1j*k*np.cos(np.deg2rad(analog_el))*horizontal)
    zrow = (np.arange(c.rows)-(c.rows-1)/2)*c.spacing
    digital = elwin*np.exp(-1j*k*zrow*np.sin(np.deg2rad(c.beam_el)))
    digital = np.repeat(digital, 2)/(np.sum(azwin)*np.sum(elwin))
    channel = 2*g["row"]+(g["panel"] >= c.active_columns//2)
    return combiner* digital[channel], combiner*digital[channel]*np.where(g["panel"] < c.active_columns//2,-1,1)


def responses(c, g, xyz, w, az, el):
    az = np.atleast_1d(np.asarray(az, dtype=float))
    el = np.broadcast_to(np.asarray(el,dtype=float),az.shape)
    ar, er = np.deg2rad(az),np.deg2rad(el)
    u = np.stack((np.cos(er)*np.cos(ar),np.cos(er)*np.sin(ar),np.sin(er)))
    a = np.exp(1j*(2*np.pi*c.frequency/c.c0)*(xyz @ u))
    if c.element_pattern == "cosine_front":
        a *= np.maximum(0,np.cos(g["phi"][:,None]-ar))*np.cos(er)
    elif c.element_pattern != "isotropic":
        raise ValueError("Unknown element pattern")
    return w[0] @ a, w[1] @ a


def calibration(c,g,w):
    az = c.steer_az+np.arange(-3,3.0001,.005)
    s,t = responses(c,g,g["xyz"],w,az,c.calibration_el)
    ratio = np.imag(t/s)
    valid = (np.abs(s) > .1*np.max(np.abs(s))) & np.isfinite(ratio)
    dv = np.diff(ratio)
    mid = int(np.argmin(np.abs(az-c.steer_az)))
    sign = np.sign(dv[min(mid,len(dv)-1)])
    if not valid[mid] or sign == 0:
        raise ValueError("Invalid calibration at boresight")
    lo=hi=mid
    while lo>0 and valid[lo-1] and sign*dv[lo-1]>1e-12: lo-=1
    while hi<len(az)-1 and valid[hi+1] and sign*dv[hi]>1e-12: hi+=1
    order=np.argsort(ratio[lo:hi+1]); x=ratio[lo:hi+1][order]; y=az[lo:hi+1][order]
    if hi-lo<10 or not (y.min() <= c.target_az <= y.max()):
        raise ValueError("Target outside monotone calibration interval")
    return dict(ratio=x,az=y,min_sum=.1*np.max(np.abs(s)))


def estimate(c,g,xyz,w,cal):
    s,t=responses(c,g,xyz,w,[c.target_az],c.target_el)
    v=np.imag(t[0]/s[0])
    if np.abs(s[0])<=cal["min_sum"] or not (cal["ratio"][0]<=v<=cal["ratio"][-1]):
        return np.nan
    return np.interp(v,cal["ratio"],cal["az"])


def peak(c,g,xyz,w):
    grid=c.steer_az+np.arange(-4,4.0001,.1)
    s,_=responses(c,g,xyz,w,grid,c.beam_el)
    i=int(np.argmax(np.abs(s)**2))
    if i in (0,len(grid)-1): return np.nan,np.nan,True
    def objective(a):
        x,_=responses(c,g,xyz,w,[a],c.beam_el)
        return -float(np.abs(x[0])**2)
    optimum=minimize_scalar(objective,bounds=(grid[i-1],grid[i+1]),method="bounded",
                            options={"xatol":1e-8})
    return optimum.x,-optimum.fun,False


def monte_carlo(c: Config):
    g=geometry(c); rng=np.random.default_rng(c.seed)
    wA=weights(c,g,"A"); wB=weights(c,g,"B")
    cal=calibration(c,g,wB)
    bB,pB,_=peak(c,g,g["xyz"],wB)
    val,vec=np.linalg.eigh(panel_correlation(c,g))
    root=vec@np.diag(np.sqrt(np.maximum(val,0)))
    beams=np.full((c.trials,4),np.nan); mono=beams.copy(); gains=beams.copy()
    edges=np.zeros_like(beams,dtype=bool)
    tangential=np.empty((c.trials,g["xyz"].shape[0]))
    az_change=np.empty_like(tangential)
    for j,w in enumerate((wA,wB)):
        b,p,edge=peak(c,g,g["xyz"],w)
        beams[:,j]=b-bB; mono[:,j]=estimate(c,g,g["xyz"],w,cal)-c.target_az
        gains[:,j]=10*np.log10(p/pB); edges[:,j]=edge
    example=None
    for i in range(c.trials):
        d=deformation(c,g,rng,root)
        if example is None: example=d
        tangential[i]=d["tangential"]
        az_change[i]=np.rad2deg(d["angle_change"])
        xyz=g["xyz"]+d["xyz"]
        for j,w in ((2,wA),(3,wB)):
            b,p,edge=peak(c,g,xyz,w)
            beams[i,j]=b-bB; mono[i,j]=estimate(c,g,xyz,w,cal)-c.target_az
            gains[i,j]=10*np.log10(p/pB); edges[i,j]=edge
    return dict(config=c,geometry=g,beam=beams,mono=mono,gain_db=gains,
                boundary=edges,example=example,cal=cal,reference_peak=bB,
                reference_power=pB,tangential=tangential,az_change=az_change)


def summary_rows(result):
    rows=[]
    for j,name in enumerate(("A_ideal_single","B_nominal_stack","C_deformed_single","D_deformed_stack")):
        for metric,key in (("beam_deg","beam"),("monopulse_deg","mono")):
            v=result[key][:,j]; x=v[np.isfinite(v)]
            rows.append(dict(scenario=name,metric=metric,bias=np.mean(x) if len(x) else np.nan,
                             std=np.std(x) if len(x) else np.nan,
                             rmse=np.sqrt(np.mean(x*x)) if len(x) else np.nan,
                             p025=np.percentile(x,2.5) if len(x) else np.nan,
                             p975=np.percentile(x,97.5) if len(x) else np.nan,
                             invalid_count=len(v)-len(x)))
    return rows


def write_csv(path,rows):
    with open(path,"w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def plot_baseline(o,out):
    c=o["config"];g=o["geometry"];d=o["example"]
    def save(name):
        plt.tight_layout();plt.savefig(out/(name+".png"),dpi=150);plt.close()
    fig=plt.figure(figsize=(7,5));ax=fig.add_subplot(projection="3d")
    for xyz,label in ((g["xyz"],"nominal"),(g["xyz"]+d["xyz"],"deformed")):
        ax.scatter(*xyz.T,s=4,label=label)
    ax.set(xlabel="X (m)",ylabel="Y (m)",zlabel="Z (m)");ax.legend();save("02_geometry")
    fig,axs=plt.subplots(2,3,figsize=(13,7));axs[0,0].plot(g["az"],d["tangential"].reshape(c.rows,c.active_columns).mean(axis=0)*1e3,"o-")
    axs[0,0].set(xlabel="Active panel azimuth (deg)",ylabel="Mean dS (mm)")
    for ax,(key,label,mult) in zip(axs.flat[1:],[("tangential","Tangential (mm)",1e3),("radial","Radial (mm)",1e3),("vertical","Vertical (mm)",1e3),("radius_change","Radius change (mm)",1e3),("angle_change","Element azimuth change (deg)",180/np.pi)]):
        ax.hist(d[key]*mult,bins=20);ax.set(xlabel=label,ylabel="Elements")
    save("03_displacement")
    az=c.steer_az+np.arange(-25,25.001,.05)
    fig,ax=plt.subplots(figsize=(8,5));fig2,ax2=plt.subplots(figsize=(8,5))
    for label,arch,xyz in (("A","A",g["xyz"]),("B","B",g["xyz"]),("C","A",g["xyz"]+d["xyz"]),("D","B",g["xyz"]+d["xyz"])):
        s,t=responses(c,g,xyz,weights(c,g,arch),az,c.target_el)
        ax.plot(az,20*np.log10(np.maximum(np.abs(s),1e-10)/np.sqrt(o["reference_power"])),label=label)
        valid=(np.abs(s)>o["cal"]["min_sum"]) & (az>=o["cal"]["az"].min()) & (az<=o["cal"]["az"].max())
        ax2.plot(az[valid],np.imag(t[valid]/s[valid]),label=label)
    ax.set(xlabel="Azimuth (deg)",ylabel="Response (dB relative to B peak)",ylim=(-60,5));ax.legend()
    fig.tight_layout();fig.savefig(out/"04_pattern.png",dpi=150);plt.close(fig)
    ax2.set(xlabel="Azimuth (deg)",ylabel="Imaginary difference / sum");ax2.legend()
    fig2.tight_layout();fig2.savefig(out/"05_monopulse_curve.png",dpi=150);plt.close(fig2)
    fig,axs=plt.subplots(1,2,figsize=(10,4))
    for ax,key,label in zip(axs,("beam","mono"),("Beam pointing error (deg)","Monopulse error (deg)")):
        v=o[key][:,2:];finite=v[np.isfinite(v)];edges=np.linspace(finite.min()-1e-9,finite.max()+1e-9,26)
        for j in range(2): ax.hist(v[:,j],bins=edges,histtype="step",label=("C","D")[j])
        ax.set(xlabel=label,ylabel="Count");ax.legend()
    save("06_errors")
    fig,ax=plt.subplots();ax.scatter(o["beam"][:,3],o["mono"][:,3],s=8)
    ax.set(xlabel="Beam error (deg)",ylabel="Monopulse error (deg)");save("07_error_relationship")
    fig,ax=plt.subplots()
    for j in (2,3):
        v=o["mono"][:,j];ok=np.isfinite(v)
        ax.plot(np.arange(1,len(v)+1),np.sqrt(np.cumsum(np.where(ok,v*v,0))/np.maximum(np.cumsum(ok),1)),label=("C","D")[j-2])
    ax.set(xlabel="Realizations",ylabel="Conditional monopulse RMSE (deg)");ax.legend();save("08_convergence")
    fig,ax=plt.subplots();x=np.arange(4);bars=np.array([[r["rmse"] for r in summary_rows(o)[i:i+2]] for i in (0,2,4,6)])
    ax.bar(x-.18,bars[:,0],.36,label="Beam");ax.bar(x+.18,bars[:,1],.36,label="Monopulse")
    ax.set(xticks=x,xticklabels=list("ABCD"),ylabel="RMSE (deg)");ax.legend();save("09_scenarios")


def frequency_study(c,out):
    rows=[]; displacement=[]
    for f in np.linspace(c.f_low,c.f_high,5):
        cc=replace(c,frequency=float(f));o=monte_carlo(cc)
        for row in summary_rows(o): rows.append(dict(frequency_hz=f,**row))
        xi=(f-c.f_low)/(c.f_high-c.f_low)
        gain=xi**2*(6-4*xi+xi**2)/3 if c.bending_shape=="uniform" else xi**2*(3-xi)/2
        displacement.append(dict(frequency_hz=f,phase_center_gain=gain,
            mean_tangential_m=float(o["tangential"].mean()),
            mean_element_sd_m=float(o["tangential"].std(axis=0).mean()),
            mean_azimuth_shift_deg=float(o["az_change"].mean()),
            mean_element_azimuth_sd_deg=float(o["az_change"].std(axis=0).mean())))
    write_csv(out/"frequency_error_summary.csv",rows)
    write_csv(out/"frequency_displacement_summary.csv",displacement)
    fig,axs=plt.subplots(1,2,figsize=(10,4))
    for metric,ax in (("beam_deg",axs[0]),("monopulse_deg",axs[1])):
        for scenario in sorted(set(r["scenario"] for r in rows)):
            subset=[r for r in rows if r["metric"]==metric and r["scenario"]==scenario]
            ax.plot([r["frequency_hz"]/1e9 for r in subset],[r["rmse"] for r in subset],"o-",label=scenario[0])
        ax.set(xlabel="Frequency (GHz)",ylabel=metric+" RMSE");ax.legend()
    fig.tight_layout();fig.savefig(out/"12_frequency_error_comparison.png",dpi=150);plt.close(fig)


def parameter_sweep(c,out,study=False):
    # Common seed at each point gives paired comparisons across the sweep.
    if study:
        grids={"sigma_tip_wavelengths":[0,.005,.01,.02,.04,.06,.08,.10],
               "mean_tip_wavelengths":[-.06,-.03,-.01,0,.01,.03,.06],
               "target_el":[0,2,4,6,8]}
        frequencies=(5.2e9,5.5e9,5.8e9)
    else:
        grids={"sigma_tip_wavelengths":[0,.005,.01,.02,.03,.05,.07,.10],
               "frequency":np.linspace(c.f_low,c.f_high,5),
               "rho_within":[0,.5,.9,.995,1],"rho_panel":[0,.35,.7,.95,1],
               "target_az":[-1.5,-1,0,.5,1,1.5],"target_el":[0,2,4,6,8]}
        frequencies=(c.frequency,)
    for parameter,values in grids.items():
        rows=[]
        for f in frequencies:
            for value in values:
                kw={parameter:float(value),"frequency":float(f)}
                if parameter=="target_el": kw.update(beam_el=float(value),calibration_el=float(value))
                v=replace(c,**kw)
                o=monte_carlo(v)
                for result in summary_rows(o):
                    rows.append(dict(parameter=parameter,value=value,frequency_hz=f,**result))
        write_csv(out/(parameter+"_summary.csv"),rows)
        fig,axs=plt.subplots(1,2,figsize=(10,4))
        for metric,ax in zip(("beam_deg","monopulse_deg"),axs):
            for f in frequencies:
                for scenario in ("B_nominal_stack","D_deformed_stack"):
                    series=[r for r in rows if r["frequency_hz"]==f and r["metric"]==metric and r["scenario"]==scenario]
                    ax.plot([r["value"] for r in series],[r["rmse"] for r in series],
                            "--" if scenario.startswith("B") else "-o",label=f"{scenario[0]} {f/1e9:.2f} GHz")
            ax.set(xlabel=parameter,ylabel=metric+" RMSE (deg)");ax.legend(fontsize=7)
        fig.tight_layout();fig.savefig(out/(parameter+"_rmse.png"),dpi=150);plt.close(fig)


def plot_receive_patterns(o,out):
    c=o["config"];g=o["geometry"];d=o["example"]
    az=np.arange(-12,12.001,.5);el=np.arange(-12,18.001,.5)
    azm,elm=np.meshgrid(az,el)
    fig=plt.figure(figsize=(12,9))
    for j,(label,arch,xyz) in enumerate((("A ideal single","A",g["xyz"]),
        ("B nominal stack","B",g["xyz"]),("C deformed single","A",g["xyz"]+d["xyz"]),
        ("D deformed stack","B",g["xyz"]+d["xyz"])),start=1):
        elevations=(c.beam_el,) if arch=="A" else (-4,0,4,8,12)
        values=[]
        for beam in elevations:
            cc=replace(c,beam_el=beam)
            s,_=responses(cc,g,xyz,weights(cc,g,arch),azm.ravel(),elm.ravel())
            values.append(20*np.log10(np.maximum(np.abs(s),1e-10)/np.sqrt(o["reference_power"])))
        db=np.maximum.reduce(values).reshape(azm.shape)
        ax=fig.add_subplot(2,2,j,projection="3d")
        ax.plot_surface(azm,elm,np.maximum(db,-35),cmap="viridis",linewidth=0)
        ax.set(xlabel="Azimuth (deg)",ylabel="Elevation (deg)",zlabel="Relative dB",title=label,zlim=(-35,3))
    fig.tight_layout();fig.savefig(out/"10_receive_patterns_3d.png",dpi=130);plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode",choices=("baseline","frequency","sensitivity","study","all"),default="baseline")
    parser.add_argument("--trials",type=int,default=500)
    parser.add_argument("--output",type=Path,default=Path("output_python"))
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    c=Config(trials=args.trials)
    if args.mode in ("baseline","all"):
        o=monte_carlo(c);write_csv(args.output/"scenario_summary.csv",summary_rows(o))
        plot_baseline(o,args.output)
        plot_receive_patterns(o,args.output)
        print("B reference peak (deg):",o["reference_peak"])
        print("D invalid estimates:",sum(~np.isfinite(o["mono"][:,3])))
    if args.mode in ("frequency","all"): frequency_study(c,args.output)
    if args.mode in ("sensitivity","all"): parameter_sweep(c,args.output)
    if args.mode in ("study","all"): parameter_sweep(c,args.output,study=True)
    if args.mode in ("study","all"):
        v=replace(c,frequency=c.f_high,target_el=6,beam_el=6,calibration_el=6,
                  sigma_tip_wavelengths=.08,mean_tip_wavelengths=.03)
        target=args.output/"strong_elevation_case";target.mkdir(exist_ok=True)
        o=monte_carlo(v);write_csv(target/"scenario_summary.csv",summary_rows(o))
        plot_baseline(o,target)
    print("Results:",args.output.resolve())


if __name__=="__main__": main()
