"""Reproducible portfolio infographic. Quantitative curves use the public Python model.
Usage from repo root: python figures/make_linkedin_figure_v2.py --summary figures/scenario_summary.csv
The bent sheet is a model-based geometric schematic, not a structural solver output.
"""
import argparse, csv, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, Arc
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

ROOT=Path(__file__).resolve().parents[1]
MODEL_DIR=ROOT/'python' if (ROOT/'python/array_study.py').exists() else ROOT/'public_array_study/python'
sys.path.insert(0,str(MODEL_DIR))
from array_study import Config, geometry, deformation, weights, responses, calibration, estimate, peak
NAVY='#0b2745'; BLUE='#0875bd'; ORANGE='#e3882e'; BG='#f2f6fa'; INK='#15314b'; MUTED='#526b80'

def gain(t): return t*t*(6-4*t+t*t)/3

def card(fig,x,y,w,h,title):
    fig.patches.append(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.006,rounding_size=.013',
      facecolor='white',edgecolor='#d4e1e9',lw=1.1,transform=fig.transFigure,zorder=-10))
    fig.text(x+.016,y+h-.034,title,fontsize=15,fontweight='bold',color=INK,va='top')

def line_style(ax):
    ax.grid(alpha=.18,lw=.65);ax.spines[['top','right']].set_visible(False)
    ax.tick_params(labelsize=8,colors=MUTED)

def bent_panel(ax, c, a, height, width, s, scale, active, deformed):
    er=np.array([np.cos(a),np.sin(a)]);et=np.array([-np.sin(a),np.cos(a)])
    # Root fixed. Shape along radial panel follows the mechanical gain function.
    tt=np.linspace(0,1,6);zz=np.linspace(-height/2,height/2,4)
    polygons=[]
    for k in range(len(tt)-1):
        for m in range(len(zz)-1):
            vertices=[]
            for t,z in ((tt[k],zz[m]),(tt[k+1],zz[m]),(tt[k+1],zz[m+1]),(tt[k],zz[m+1])):
                p=(c.radius+width*t)*er+(scale*s*gain(t) if deformed else 0)*et
                vertices.append([p[0],p[1],z])
            polygons.append(vertices)
    color=ORANGE if active and deformed else BLUE if active else '#d8e1e9'
    edge='#aebdcc' if not active else '#194e76'
    ax.add_collection3d(Poly3DCollection(polygons,facecolors=color,edgecolors=edge,
       linewidths=.13 if active else .08,alpha=.82 if active else .14,zsort='average'))
    if active:
        # Illustrative row markers on the panel-mean tip locus; RF retains row residuals.
        zrows=(np.arange(c.rows)-(c.rows-1)/2)*c.spacing
        for z in zrows:
            p=(c.radius+width)*er+(scale*s if deformed else 0)*et
            ax.scatter(p[0],p[1],z,s=5,color='#fff4d8' if deformed else '#c4e8ff',depthshade=False)

def draw_array(ax,c,g,d,deformed,scale,width):
    height=(c.rows-1)*c.spacing
    az=2*np.pi*(np.arange(c.total_columns)+.5)/c.total_columns
    active_ids=set(range(50,60))|set(range(0,10)) # 20 contiguous positions centred on +x.
    for j,a in enumerate(az):
        if j in active_ids:
            p=(j-50) if j>=50 else j+10
            # Map physical sorted azimuth -57..+57 degrees to row-major active index.
            s=np.mean(d['tangential'].reshape(c.rows,c.active_columns)[:,p])
        else:s=0.0
        bent_panel(ax,c,a,height,width,s,scale,j in active_ids,deformed)
    ax.view_init(elev=23,azim=12)
    lim=c.radius+width+.035
    ax.set(xlim=(-lim,lim),ylim=(-lim,lim),zlim=(-height*.68,height*.68))
    ax.set_box_aspect((1,1,.95),zoom=1.04);ax.set_axis_off()

def zoom_panel(ax,c,d,scale,width):
    # Use the strongest panel in the selected realization so bending is legible.
    vals=d['tangential'].reshape(c.rows,c.active_columns).mean(axis=0)
    idx=int(np.argmax(np.abs(vals)));s=vals[idx]
    tt=np.linspace(0,1,120);y=scale*s*gain(tt)/width
    sign=1 if s>=0 else -1
    ax.fill_between(tt,0,y,color='#cae2f5',alpha=.55)
    ax.plot(tt,np.zeros_like(tt),color='#879bad',lw=2,ls='--',label='Nominal radial line')
    ax.plot(tt,y,color=BLUE,lw=3,label='Bent panel centreline')
    ax.scatter([0,1],[0,y[-1]],c=[NAVY,ORANGE],s=[45,48],zorder=5)
    slope=scale*s*(4/3)/width
    ax.plot([.77,1.13],[y[-1],y[-1]],color=MUTED,lw=1.0,ls=':')
    ax.plot([.78,1.1],[y[-1]+slope*(.78-1),y[-1]+slope*(1.1-1)],color=ORANGE,lw=1.4,ls=':')
    ang=np.degrees(np.arctan(slope))
    ax.add_patch(Arc((1,y[-1]),.36,.36,theta1=min(ang,0),theta2=max(ang,0),color=ORANGE,lw=1.5))
    ax.annotate('free edge',xy=(1,y[-1]),xytext=(.42,y[-1]+sign*.28),fontsize=8,
                arrowprops=dict(arrowstyle='->',color=ORANGE),color=INK)
    ax.annotate('fixed edge',xy=(0,0),xytext=(.04,-sign*.38),fontsize=8,
                arrowprops=dict(arrowstyle='->',color=NAVY),color=INK)
    ax.text(.63,y[-1]+sign*.18,r'$\alpha_{\mathrm{display}}$',fontsize=10,color=ORANGE)
    ax.set(xlim=(-.07,1.15),ylim=(-1.05,1.05),xlabel='Position from root to tip (normalized)',ylabel='Tangential offset / panel depth')
    ax.tick_params(labelsize=7);line_style(ax)
    return idx

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--summary',type=Path,required=True)
    ap.add_argument('--output',type=Path,default=Path('figures/array_deformation_linkedin.png'))
    args=ap.parse_args();rows=list(csv.DictReader(args.summary.open()))
    def metric(m):return next(r for r in rows if r['scenario']=='D_deformed_stack' and r['metric']==m)
    c=Config();g=geometry(c);w=weights(c,g,'B');cal=calibration(c,g,w)
    rng=np.random.default_rng(c.seed)
    for i in range(500):
        d=deformation(c,g,rng);xyz=g['xyz']+d['xyz'];est=estimate(c,g,xyz,w,cal)
        if np.isfinite(est) and .10<=abs(est-c.target_az)<=.14:break
    else:raise RuntimeError('No suitable realization')
    az=np.linspace(-5,5,1001);sn,dn=responses(c,g,g['xyz'],w,az,c.target_el)
    sd,dd=responses(c,g,xyz,w,az,c.target_el)
    ratio_nom=np.imag(dn/sn);ratio_def=np.imag(dd/sd)
    actual_ratio=np.interp(c.target_az,az,ratio_def)
    # The panel surface follows g(t); RF curves use the unchanged model coordinates.
    display_scale=30;panel_depth=2*c.reference_wavelength
    plt.rcParams.update({'font.family':'DejaVu Sans'})
    fig=plt.figure(figsize=(18,12),facecolor=BG)
    fig.patches.append(FancyBboxPatch((0,.906),1,.094,boxstyle='square,pad=0',
      facecolor=NAVY,edgecolor='none',transform=fig.transFigure,zorder=-10))
    fig.text(.035,.970,'When antenna geometry moves',color='white',fontsize=28,fontweight='bold',va='top')
    fig.text(.035,.928,'Correlated panel bending  •  sum/difference response  •  monopulse angle error',color='#c5d7e8',fontsize=12)
    card(fig,.026,.297,.459,.59,'Computed RF response: nominal B vs. deformed D')
    card(fig,.499,.509,.474,.378,'Cylindrical array: fixed roots, free tips')
    card(fig,.499,.297,.474,.198,'Monopulse mapping: 1.000° → 1.105°')
    card(fig,.026,.066,.947,.215,'Monte Carlo results  |  synthetic C-band example')
    a1=fig.add_axes([.086,.628,.366,.202]);a2=fig.add_axes([.086,.393,.366,.170]);a3=fig.add_axes([.555,.331,.370,.081])
    for ax,y0,y1,label,lim in ((a1,abs(sn)/max(abs(sn)),abs(sd)/max(abs(sn)),'|Σ| / nominal peak',(-5,5)),
      (a2,np.imag(dn)/max(abs(sn)),np.imag(dd)/max(abs(sn)),'Im(Δ) / nominal peak',(-5,5)),
      (a3,ratio_nom,ratio_def,'Im(Δ / Σ)',(-3,3))):
        ax.plot(az,y0,color=BLUE,lw=2,label='Nominal B')
        ax.plot(az,y1,color=ORANGE,lw=2,ls='--',label='Deformed D')
        ax.set(xlim=lim,ylabel=label);line_style(ax)
    a1.legend(fontsize=8,loc='upper right',ncol=2);a2.axhline(0,lw=.7,color=MUTED)
    a1.set_xlabel('Target azimuth (deg)',fontsize=9)
    a2.set_xlabel('Target azimuth (deg)',fontsize=9)
    a3.set_xlabel('Azimuth (deg)',fontsize=8)
    a3.scatter([c.target_az,est],[actual_ratio]*2,c=[ORANGE,BLUE],s=27,zorder=4)
    a3.plot([c.target_az,est],[actual_ratio]*2,lw=1,color=INK)
    a3.axvline(c.target_az,color=ORANGE,ls=':',lw=.8);a3.axvline(est,color=BLUE,ls=':',lw=.8)
    fig.text(.084,.344,'Blue: nominal B   •   Orange dashed: one deformed realization D',fontsize=8.8,color=INK)
    fig.text(.084,.319,'The sum beam changes little; the angle estimate can still shift.',fontsize=8.4,color=MUTED)
    # Two arrays occupy the main geometry card; the zoom is a separate inset.
    axn=fig.add_axes([.516,.595,.158,.258],projection='3d')
    axd=fig.add_axes([.675,.595,.158,.258],projection='3d')
    draw_array(axn,c,g,d,False,display_scale,panel_depth)
    draw_array(axd,c,g,d,True,display_scale,panel_depth)
    fig.text(.557,.578,'NOMINAL',fontweight='bold',fontsize=9,color=BLUE)
    fig.text(.715,.578,'DEFORMED',fontweight='bold',fontsize=9,color=ORANGE)
    z=fig.add_axes([.844,.615,.114,.188],facecolor='#f4f8fb')
    selected=zoom_panel(z,c,d,display_scale,panel_depth)
    z.set_title('Panel bending (top view)',fontsize=8.5,fontweight='bold',color=INK,pad=5)
    fig.text(.516,.543,'60 radial panels  •  20 active  •  16 elements/panel',fontsize=9,color=INK)
    fig.text(.516,.520,'Fixed root; the free edge moves as the panel bends (display ×30).',fontsize=8.3,color=MUTED)
    fig.text(.517,.425,f'One selected realization  •  angular error {est-c.target_az:+.3f}°',fontsize=8.7,color=MUTED)
    beam=metric('beam_deg');mono=metric('monopulse_deg')
    fig.text(.046,.213,'500 realizations  •  5.8 GHz  •  0.02 λref tip SD',fontsize=10,color=MUTED)
    fig.text(.046,.169,f'{float(beam["rmse"]):.3f}°',fontsize=30,fontweight='bold',color=BLUE)
    fig.text(.046,.131,'beam-peak shift RMSE',fontsize=9,color=INK)
    fig.text(.303,.169,f'{float(mono["rmse"]):.3f}°',fontsize=30,fontweight='bold',color=ORANGE)
    fig.text(.303,.131,'monopulse angle RMSE',fontsize=9,color=INK)
    fig.text(.555,.174,'Project assessment: acceptable modeled error.',fontsize=11,fontweight='bold',color=INK)
    fig.text(.555,.139,'Team-reported for the intended system; criteria are not public.',fontsize=9,color=MUTED)
    fig.lines.append(plt.Line2D([.046,.952],[.114,.114],transform=fig.transFigure,color='#dce5ed',lw=.8))
    fig.text(.046,.095,'MODEL',fontsize=8.2,fontweight='bold',color=BLUE)
    fig.text(.046,.078,'Correlated panel and element displacement',fontsize=8.1,color=INK)
    fig.text(.303,.095,'ESTIMATION',fontsize=8.2,fontweight='bold',color=ORANGE)
    fig.text(.303,.078,'Deformed ratio, nominal B calibration',fontsize=8.1,color=INK)
    fig.text(.555,.095,'CHECKS AND LIMIT',fontsize=8.2,fontweight='bold',color=INK)
    fig.text(.555,.078,'Generator statistics checked; physical validation pending',fontsize=8.1,color=INK)
    fig.text(.032,.033,'Bent sheets are schematic; RF model shifts phase centres without rotating element normals. No hardware/structural validation.',fontsize=9,color=MUTED)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(args.output,dpi=170,facecolor=BG);plt.close(fig)
    print(f'image={args.output} realization={i+1} estimate={est:.6f} zoom_panel={selected}')
if __name__=='__main__':main()
