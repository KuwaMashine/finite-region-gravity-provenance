#!/usr/bin/env python3
"""Build the normalized source-only H I radial profile from moment-0 pixels.

Input pixels are exported from the public source/moment-0 product, not from a
kinematic model.  Geometry uses the source-catalogue moment-0 ellipse
(ell_maj, ell_min, ell_pa).  No Vrot or fitted tilted-ring geometry is accepted.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd

REQ_PIX={'galaxy_id','dx_east_arcsec','dy_north_arcsec','moment0_weight'}
REQ_GEO={'galaxy_id','ell_maj_arcsec','ell_min_arcsec','ell_pa_deg','pixel_scale_arcsec'}
FORBIDDEN={'Vrot_model','vrot_kms','gobs','gbar','delta_ref','rar_residual','Vmax','Vflat','Inc_model','PA_model_g','SD_FO_model'}


def one_profile(pix:pd.DataFrame, geo:pd.Series):
    q=float(geo.ell_min_arcsec)/float(geo.ell_maj_arcsec)
    if not (np.isfinite(q) and 0<q<=1): raise ValueError('source ellipse must imply 0<ell_min/ell_maj<=1')
    ps=float(geo.pixel_scale_arcsec)
    if not (np.isfinite(ps) and ps>0): raise ValueError('pixel_scale_arcsec must be positive')
    x=pd.to_numeric(pix.dx_east_arcsec,errors='coerce').to_numpy(float)
    y=pd.to_numeric(pix.dy_north_arcsec,errors='coerce').to_numpy(float)
    w=pd.to_numeric(pix.moment0_weight,errors='coerce').to_numpy(float)
    ok=np.isfinite(x)&np.isfinite(y)&np.isfinite(w)
    x,y,w=x[ok],y[ok],w[ok]
    if len(w)<8: raise ValueError('fewer than 8 finite source-map pixels')
    # Position angle is taken in the standard sky convention: north through east.
    pa=np.deg2rad(float(geo.ell_pa_deg))
    major=x*np.sin(pa)+y*np.cos(pa)
    minor=x*np.cos(pa)-y*np.sin(pa)
    r=np.sqrt(major**2+(minor/q)**2)
    # The exported source pixel table must already be source/mask scoped. Negative
    # annular integrated weights are not clipped after the fact; they fail closed.
    rmax=float(np.max(r))
    nb=int(np.ceil(rmax/ps))
    if nb<8: raise ValueError('source radial support contains fewer than 8 pixel-scale annuli')
    edges=np.arange(nb+1,dtype=float)*ps
    idx=np.searchsorted(edges,r,side='right')-1
    idx=np.clip(idx,0,nb-1)
    sums=np.bincount(idx,weights=w,minlength=nb)
    counts=np.bincount(idx,minlength=nb)
    # Drop only trailing empty annuli. Interior annuli with no source-mask
    # pixels are retained with zero source weight; no interpolation is inserted.
    nz=np.flatnonzero(counts>0)
    if not len(nz): raise ValueError('no populated source annuli')
    last=int(nz[-1]); sums=sums[:last+1]; counts=counts[:last+1]; edges=edges[:last+2]
    if len(sums)<8: raise ValueError('fewer than 8 radial annuli across source support')
    if np.any(sums<0): raise ValueError('negative source annular moment-0 weight; no clipping is allowed')
    area=np.pi*(edges[1:]**2-edges[:-1]**2)
    dens=sums/area
    if not np.isfinite(dens).all() or dens.sum()<=0: raise ValueError('invalid source radial weight profile')
    centres=0.5*(edges[1:]+edges[:-1])
    return centres,dens,{'q_source':q,'pixel_scale_arcsec':ps,'n_annuli':len(centres),'rmax_arcsec':float(edges[-1])}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--pixels',type=Path,required=True)
    ap.add_argument('--geometry',type=Path,required=True)
    ap.add_argument('--out-profile',type=Path,required=True)
    ap.add_argument('--out-report',type=Path,required=True)
    a=ap.parse_args()
    px=pd.read_csv(a.pixels,dtype={'galaxy_id':str}); ge=pd.read_csv(a.geometry,dtype={'galaxy_id':str})
    for df,n,req in [(px,'pixels',REQ_PIX),(ge,'geometry',REQ_GEO)]:
        miss=req-set(df.columns)
        if miss: raise SystemExit(f'{n} missing columns {sorted(miss)}')
        leak=FORBIDDEN&set(df.columns)
        if leak: raise SystemExit(f'SOURCE LOCK: {n} contains kinematic/response columns {sorted(leak)}')
    if ge.galaxy_id.duplicated().any(): raise SystemExit('one geometry row per galaxy required')
    rows=[]; reps=[]; fails=[]
    for g in ge.sort_values('galaxy_id').itertuples(index=False):
        gid=str(g.galaxy_id); z=px[px.galaxy_id==gid]
        try:
            if z.empty: raise ValueError('no source-map pixels')
            r,d,meta=one_profile(z,pd.Series(g._asdict()))
            for rr,dd in zip(r,d): rows.append({'galaxy_id':gid,'rad_hi_source_arcsec':float(rr),'sigma_hi_source_weight':float(dd)})
            reps.append({'galaxy_id':gid,**meta})
        except Exception as e: fails.append({'galaxy_id':gid,'reason':str(e)})
    if not rows: raise SystemExit('no source profile passed: '+json.dumps(fails[:8],sort_keys=True))
    pd.DataFrame(rows).to_csv(a.out_profile,index=False,float_format='%.12g')
    report={'geometry':'source-catalogue moment-0 ellipse only; q=ell_min/ell_maj; PA north-through-east',
            'radial_bin_width':'one exported source-image pixel','response_products_used':False,
            'n_pass':len(reps),'n_fail':len(fails),'passes':reps,'failures':fails}
    a.out_report.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'n_pass':len(reps),'n_fail':len(fails)},sort_keys=True))
if __name__=='__main__': main()
