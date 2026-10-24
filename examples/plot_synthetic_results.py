"""Create independent Python figures from the synthetic model (no MATLAB required)."""
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src" / "python"))
from run_demo import Config, geometry, weights, response, displacement, run

def main():
    c=Config()
    g=geometry(c)
    rng=np.random.default_rng(c.seed)
    delta,_=displacement(c,g,rng)
    w=weights(c,g,"fixed_analog")
    az=np.linspace(-8,8,801)
    s0,d0=response(c,g["xyz"],w,az,c.beam_el_deg)
    s1,d1=response(c,g["xyz"]+delta,w,az,c.beam_el_deg)
    out=Path("results/figures"); out.mkdir(parents=True,exist_ok=True)
    fig,ax=plt.subplots(figsize=(8,4.5))
    for s,label in [(s0,"Nominal"),(s1,"Deformed (one synthetic draw)")]:
        ax.plot(az,20*np.log10(np.maximum(abs(s),1e-12)/max(abs(s0))),label=label)
    ax.set(xlabel="Azimuth (deg)",ylabel="Sum response (dB)",title="Synthetic array response")
    ax.grid(alpha=.25);ax.legend();fig.tight_layout()
    fig.savefig(out/"sum_pattern.png",dpi=200);plt.close(fig)
    fig,ax=plt.subplots(figsize=(8,4.5))
    ax.plot(az,np.imag(d0/s0),label="Nominal")
    ax.plot(az,np.imag(d1/s1),label="Deformed (one synthetic draw)")
    ax.axvline(c.target_az_deg,color="gray",ls="--",label="Target azimuth")
    ax.set(xlabel="Azimuth (deg)",ylabel="Im(Difference / Sum)",title="Synthetic monopulse ratio")
    ax.grid(alpha=.25);ax.legend();fig.tight_layout()
    fig.savefig(out/"monopulse_ratio.png",dpi=200);plt.close(fig)
    _,_,rows=run(c)
    names=("deformed_steered","deformed_fixed")
    fig,ax=plt.subplots(figsize=(8,4.5))
    for name in names:
        errors=[r["monopulse_error_deg"] for r in rows if r["scenario"]==name and np.isfinite(r["monopulse_error_deg"])]
        if errors: ax.hist(errors,bins=18,alpha=.55,label=name)
    ax.set(xlabel="Monopulse error (deg)",ylabel="Count",title="Synthetic Monte Carlo errors")
    ax.grid(alpha=.2);ax.legend();fig.tight_layout()
    fig.savefig(out/"monopulse_error_histogram.png",dpi=200);plt.close(fig)
    print("Saved three synthetic figures in",out)

if __name__=="__main__":main()
