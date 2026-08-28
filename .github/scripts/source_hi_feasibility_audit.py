#!/usr/bin/env python3
"""Deterministic source-only H I feasibility audit for the frozen v4.99 profiler.

Selects 10 qflag=0 sources per final validation release by SHA256(name), downloads
only the resolved CADC source_data moment-0 and mask products, and applies the
exact byte-verified v4.99 `one_profile` function. The audit measures prospective
source-side attrition only; it does not query or read any kinematic/response data.
"""
from __future__ import annotations
import hashlib, importlib.util, json, math, urllib.request
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.wcs.utils import proj_plane_pixel_scales

ROOT=Path.cwd(); OUT=ROOT/'post-stage0/source-hi-feasibility-audit'; TMP=Path('/tmp/wallaby_source_hi_feasibility')
CAT=ROOT/'post-stage0/source-only-probe/source_catalogue_safe.csv'; RES=ROOT/'post-stage0/source-plane-resolution/SOURCE_PLANE_RESOLUTION.json'
ENG=ROOT/'post-stage0/frozen-v4.99-source-engine/source_moment0_profile.py'; BASE='https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/data/pub/WALLABY'
SAMPLE_PER_FIELD=10
RELEASES={'NGC4808':'NGC 4808 TR1','NGC5044':'NGC 5044 TR3','Vela':'Vela TR1'}

spec=importlib.util.spec_from_file_location('frozen_source_moment0_profile',ENG); frozen=importlib.util.module_from_spec(spec); spec.loader.exec_module(frozen)

def sha256_bytes(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def key(name): return hashlib.sha256(str(name).encode()).hexdigest()

def download(url,p):
    req=urllib.request.Request(url,headers={'User-Agent':'finite-region-gravity-hi-feasibility/1'})
    with urllib.request.urlopen(req,timeout=180) as r, Path(p).open('wb') as f:
        while True:
            b=r.read(1<<20)
            if not b: break
            f.write(b)
    if Path(p).stat().st_size==0: raise RuntimeError('zero-byte source product')

def filename_suffix(field):
    return {'NGC4808':'NGC_4808_TR1','NGC5044':'NGC_5044_TR3','Vela':'Vela_TR1'}[field]

def obs_from_name(name): return str(name).replace('WALLABY ','WALLABY_')

def pixel_inputs(mom,mask,header,row):
    m=np.asarray(mom).squeeze(); k=np.asarray(mask).squeeze()
    if m.ndim!=2 or k.ndim!=3 or tuple(k.shape[-2:])!=tuple(m.shape): raise ValueError(f'unexpected moment0/mask geometry {m.shape} {k.shape}')
    mask2=np.any(k>0,axis=0); yy,xx=np.nonzero(mask2 & np.isfinite(m))
    if len(xx)<8: raise ValueError('fewer than 8 finite source-map pixels')
    w=WCS(header).celestial; sc=w.pixel_to_world(xx.astype(float),yy.astype(float))
    center=SkyCoord(float(row.ra)*u.deg,float(row.dec)*u.deg,frame='icrs').transform_to(sc.frame); de,dn=center.spherical_offsets_to(sc)
    pix=pd.DataFrame({'galaxy_id':row['name'],'dx_east_arcsec':de.to_value(u.arcsec),'dy_north_arcsec':dn.to_value(u.arcsec),'moment0_weight':m[yy,xx].astype(float)})
    scales=np.abs(np.asarray(proj_plane_pixel_scales(w),float))*3600.; ps=float(np.mean(scales))
    geo=pd.Series({'galaxy_id':row['name'],'ell_maj_arcsec':6.*float(row.ell_maj),'ell_min_arcsec':6.*float(row.ell_min),'ell_pa_deg':float(row.ell_pa),'pixel_scale_arcsec':ps})
    return pix,geo,{'moment0_shape':list(m.shape),'mask_shape':list(k.shape),'n_mask_finite_pixels':int(len(xx)),
                    'moment0_inside_mask_sum':float(np.sum(m[yy,xx])),'moment0_inside_mask_negative_pixels':int(np.sum(m[yy,xx]<0)),
                    'pixel_scale_arcsec':ps,'wcs_frame':str(sc.frame.name)}

def signed_annulus_diag(pix,geo):
    # Diagnostic-only reproduction of the exact frozen binning, to quantify why a
    # negative-annulus rejection occurred. It never modifies or replaces the rule.
    q=float(geo.ell_min_arcsec)/float(geo.ell_maj_arcsec); ps=float(geo.pixel_scale_arcsec)
    x=pd.to_numeric(pix.dx_east_arcsec,errors='coerce').to_numpy(float); y=pd.to_numeric(pix.dy_north_arcsec,errors='coerce').to_numpy(float); w=pd.to_numeric(pix.moment0_weight,errors='coerce').to_numpy(float)
    ok=np.isfinite(x)&np.isfinite(y)&np.isfinite(w); x,y,w=x[ok],y[ok],w[ok]
    pa=np.deg2rad(float(geo.ell_pa_deg)); major=x*np.sin(pa)+y*np.cos(pa); minor=x*np.cos(pa)-y*np.sin(pa); r=np.sqrt(major**2+(minor/q)**2)
    nb=int(np.ceil(float(np.max(r))/ps)); edges=np.arange(nb+1,dtype=float)*ps; idx=np.clip(np.searchsorted(edges,r,side='right')-1,0,nb-1)
    sums=np.bincount(idx,weights=w,minlength=nb); counts=np.bincount(idx,minlength=nb); nz=np.flatnonzero(counts>0)
    if len(nz): sums=sums[:int(nz[-1])+1]
    pos=sums[sums>0]; scale=float(np.median(pos)) if len(pos) else None
    return {'n_annuli_to_last_populated':int(len(sums)),'n_negative_integrated_annuli':int(np.sum(sums<0)),'min_integrated_annulus_weight':float(np.min(sums)) if len(sums) else None,
            'median_positive_integrated_annulus_weight':scale,'min_over_median_positive':(float(np.min(sums))/scale if scale and scale>0 else None)}

def main():
    OUT.mkdir(parents=True,exist_ok=True); TMP.mkdir(parents=True,exist_ok=True)
    cat=pd.read_csv(CAT); res=json.loads(RES.read_text()); rows=[]; selections={}
    for field,release in RELEASES.items():
        z=cat[(cat.team_release==release)&(pd.to_numeric(cat.qflag,errors='coerce')==0)].copy()
        for c in ('ra','dec','dist_h','log_m_hi_corr','ell_maj','ell_min','ell_pa'): z=z[np.isfinite(pd.to_numeric(z[c],errors='coerce'))]
        z['selection_hash']=z['name'].map(key); z=z.sort_values(['selection_hash','name']).head(SAMPLE_PER_FIELD)
        selections[field]=z[['name','selection_hash']].to_dict('records'); suffix=filename_suffix(field); productID=res['fields'][field]['chosen']['productID']
        for rr in z.itertuples(index=False):
            name=str(rr.name); obs=obs_from_name(name); prefix=f'{obs}_{suffix}'; d=TMP/hashlib.sha256(name.encode()).hexdigest()[:12]; d.mkdir(exist_ok=True)
            rec={'field_id':field,'source_release':release,'galaxy_id':name,'selection_hash':key(name),'productID':productID,'response_products_used':False}
            try:
                mp=d/'mom0.fits'; kp=d/'mask.fits'; download(f'{BASE}/{prefix}_mom0.fits',mp); download(f'{BASE}/{prefix}_mask.fits',kp)
                with fits.open(mp,memmap=False) as h:
                    hh=h[0] if h[0].data is not None else h[1]; mom=np.asarray(hh.data).copy(); mh=hh.header.copy()
                with fits.open(kp,memmap=False) as h:
                    hh=h[0] if h[0].data is not None else h[1]; mask=np.asarray(hh.data).copy()
                pix,geo,struct=pixel_inputs(mom,mask,mh,pd.Series(rr._asdict())); rec.update(struct); rec['mom0_sha256']=sha256_bytes(mp); rec['mask_sha256']=sha256_bytes(kp)
                try:
                    centres,dens,meta=frozen.one_profile(pix,geo); rec.update({'status':'pass','n_profile_annuli':int(len(centres)),'min_profile_density':float(np.min(dens)),'profile_meta':meta})
                except Exception as e:
                    rec.update({'status':'profile_rejected','reason':str(e),'signed_annulus_diagnostic':signed_annulus_diag(pix,geo)})
            except Exception as e:
                rec.update({'status':'product_or_geometry_failure','reason':str(e)})
            rows.append(rec)
    df=pd.DataFrame([{k:(json.dumps(v,sort_keys=True) if isinstance(v,(dict,list)) else v) for k,v in r.items()} for r in rows]); df.to_csv(OUT/'source_hi_feasibility_rows.csv',index=False)
    summary={}
    for field in RELEASES:
        z=[r for r in rows if r['field_id']==field]; counts={s:sum(r['status']==s for r in z) for s in sorted(set(r['status'] for r in z))}
        summary[field]={'n_sampled':len(z),'status_counts':counts,'pass_fraction':sum(r['status']=='pass' for r in z)/len(z)}
    total_pass=sum(r['status']=='pass' for r in rows)
    report={'schema':'WALLABY_SOURCE_HI_FEASIBILITY_AUDIT_v1','scope':'deterministic source-only qflag=0 sample; exact frozen v4.99 one_profile; no response products queried',
            'selection_rule':f'within each final validation release, qflag=0 finite source rows sorted by SHA256(name), first {SAMPLE_PER_FIELD}',
            'sample_per_field':SAMPLE_PER_FIELD,'selections':selections,'field_summary':summary,'n_total':len(rows),'n_pass':total_pass,'overall_pass_fraction':total_pass/len(rows),
            'response_products_used':False,'frozen_profiler_sha256':sha256_bytes(ENG),
            'decision_note':'This audit measures whether the prospective negative-annulus/no-clipping rule is operationally viable. It does not authorize changing that rule.'}
    p=OUT/'SOURCE_HI_FEASIBILITY_AUDIT.json'; p.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    with (OUT/'SHA256SUMS.txt').open('w') as f:
        for q in [p,OUT/'source_hi_feasibility_rows.csv']: f.write(f'{sha256_bytes(q)}  {q.relative_to(ROOT)}\n')
    print(json.dumps(report,indent=2,sort_keys=True))
if __name__=='__main__': main()
