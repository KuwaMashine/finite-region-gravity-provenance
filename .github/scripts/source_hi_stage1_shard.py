#!/usr/bin/env python3
"""Build one deterministic shard of the WALLABY Stage-1 source-only H I export.

Only final Phase-2 source catalogue rows and public CADC source_data moment-0/mask
products are read. The H I radial profile is produced by the exact byte-verified
v4.99 one_profile implementation. No kinematic plane, rotation value or response
column is queried or accepted.
"""
from __future__ import annotations
import argparse, csv, hashlib, importlib.util, math, time, urllib.request
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.wcs.utils import proj_plane_pixel_scales

ROOT=Path.cwd()
CAT=ROOT/'post-stage0/source-only-probe/source_catalogue_safe.csv'
ENG=ROOT/'post-stage0/frozen-v4.99-source-engine/source_moment0_profile.py'
BASE='https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/data/pub/WALLABY'
RELEASES={
 'NGC 4808 TR1':('NGC4808','NGC_4808_TR1'),
 'NGC 5044 TR3':('NGC5044','NGC_5044_TR3'),
 'Vela TR1':('Vela','Vela_TR1'),
}
REQ_FINITE=('ra','dec','dist_h','log_m_hi_corr','ell_maj','ell_min','ell_pa')

spec=importlib.util.spec_from_file_location('frozen_source_moment0_profile',ENG)
frozen=importlib.util.module_from_spec(spec); spec.loader.exec_module(frozen)

def sha256(p:Path)->str:
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def shard_key(name:str,release:str,n:int)->int:
    d=hashlib.sha256((name+'\0'+release).encode()).digest()
    return int.from_bytes(d[:8],'big') % n

def download(url:str,p:Path):
    last=None
    for attempt in range(4):
        try:
            req=urllib.request.Request(url,headers={'User-Agent':'finite-region-gravity-stage1-hi/1'})
            with urllib.request.urlopen(req,timeout=240) as r, p.open('wb') as f:
                while True:
                    b=r.read(1<<20)
                    if not b: break
                    f.write(b)
            if p.stat().st_size<=0: raise RuntimeError('zero-byte product')
            return
        except Exception as e:
            last=e
            try: p.unlink()
            except FileNotFoundError: pass
            if attempt<3: time.sleep(2**attempt)
    raise last

def finite_row(r)->bool:
    if float(r.qflag)!=0.0: return False
    for c in REQ_FINITE:
        try:
            if not math.isfinite(float(getattr(r,c))): return False
        except Exception: return False
    return True

def project_pixels(mom,mask,hdr,r):
    m=np.asarray(mom).squeeze(); k=np.asarray(mask).squeeze()
    if m.ndim!=2 or k.ndim!=3 or tuple(k.shape[-2:])!=tuple(m.shape):
        raise ValueError(f'unexpected moment0/mask geometry {m.shape} vs {k.shape}')
    mask2=np.any(k>0,axis=0); yy,xx=np.nonzero(mask2 & np.isfinite(m))
    if len(xx)<8: raise ValueError('fewer than 8 finite source-map pixels')
    wcs=WCS(hdr).celestial
    sc=wcs.pixel_to_world(xx.astype(float),yy.astype(float))
    center=SkyCoord(float(r.ra)*u.deg,float(r.dec)*u.deg,frame='icrs').transform_to(sc.frame)
    de,dn=center.spherical_offsets_to(sc)
    pix=pd.DataFrame({'galaxy_id':str(r.name),'dx_east_arcsec':de.to_value(u.arcsec),
                      'dy_north_arcsec':dn.to_value(u.arcsec),'moment0_weight':m[yy,xx].astype(float)})
    scales=np.abs(np.asarray(proj_plane_pixel_scales(wcs),float))*3600.0
    ps=float(np.mean(scales))
    geo=pd.Series({'galaxy_id':str(r.name),'ell_maj_arcsec':6.0*float(r.ell_maj),
                   'ell_min_arcsec':6.0*float(r.ell_min),'ell_pa_deg':float(r.ell_pa),
                   'pixel_scale_arcsec':ps})
    struct={'moment0_shape':'x'.join(map(str,m.shape)),'mask_shape':'x'.join(map(str,k.shape)),
            'n_mask_finite_pixels':int(len(xx)),'pixel_scale_arcsec':ps,
            'moment0_inside_mask_sum':float(np.sum(m[yy,xx])),
            'moment0_inside_mask_negative_pixels':int(np.sum(m[yy,xx]<0))}
    return pix,geo,struct

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--shard-index',type=int,required=True); ap.add_argument('--shard-count',type=int,required=True); ap.add_argument('--out-dir',type=Path,required=True)
    a=ap.parse_args()
    if not (0<=a.shard_index<a.shard_count): raise SystemExit('invalid shard')
    a.out_dir.mkdir(parents=True,exist_ok=True); tmp=Path('/tmp')/f'wallaby_stage1_hi_{a.shard_index:03d}'; tmp.mkdir(exist_ok=True)
    df=pd.read_csv(CAT)
    df=df[df.team_release.isin(RELEASES)].copy()
    all_eligible=[]
    for r in df.itertuples(index=False):
        if finite_row(r) and shard_key(str(r.name),str(r.team_release),a.shard_count)==a.shard_index: all_eligible.append(r)
    statuses=[]; profiles=[]
    for ii,r in enumerate(all_eligible,1):
        name=str(r.name); release=str(r.team_release); field,suffix=RELEASES[release]; obs=name.replace('WALLABY ','WALLABY_'); prefix=f'{obs}_{suffix}'
        rec={'galaxy_id':name,'team_release':release,'field_id':field,'shard_index':a.shard_index,'status':'pending','response_products_used':False}
        d=tmp/hashlib.sha256((name+'\0'+release).encode()).hexdigest()[:16]; d.mkdir(exist_ok=True)
        try:
            mp=d/'mom0.fits'; kp=d/'mask.fits'
            mu=f'{BASE}/{prefix}_mom0.fits'; ku=f'{BASE}/{prefix}_mask.fits'
            download(mu,mp); download(ku,kp)
            with fits.open(mp,memmap=False) as h:
                hh=h[0] if h[0].data is not None else h[1]; mom=np.asarray(hh.data).copy(); hdr=hh.header.copy()
            with fits.open(kp,memmap=False) as h:
                hh=h[0] if h[0].data is not None else h[1]; mask=np.asarray(hh.data).copy()
            pix,geo,struct=project_pixels(mom,mask,hdr,r); rec.update(struct)
            rec.update({'mom0_url':mu,'mask_url':ku,'mom0_sha256':sha256(mp),'mask_sha256':sha256(kp),
                        'mom0_bytes':mp.stat().st_size,'mask_bytes':kp.stat().st_size})
            try:
                centres,dens,meta=frozen.one_profile(pix,geo)
                rec.update({'status':'pass','n_profile_annuli':int(len(centres)),'profile_q_source':float(meta['q_source']),
                            'profile_rmax_arcsec':float(meta['rmax_arcsec'])})
                for rr,dd in zip(centres,dens): profiles.append({'galaxy_id':name,'rad_hi_source_arcsec':float(rr),'sigma_hi_source_weight':float(dd)})
            except Exception as e:
                rec.update({'status':'profile_rejected','reason':str(e)})
        except Exception as e:
            rec.update({'status':'product_or_geometry_failure','reason':str(e)})
        statuses.append(rec)
        if ii%10==0: print(f'shard {a.shard_index}: {ii}/{len(all_eligible)}',flush=True)
    sf=['galaxy_id','team_release','field_id','shard_index','status','reason','response_products_used','n_profile_annuli','profile_q_source','profile_rmax_arcsec','moment0_shape','mask_shape','n_mask_finite_pixels','pixel_scale_arcsec','moment0_inside_mask_sum','moment0_inside_mask_negative_pixels','mom0_url','mask_url','mom0_sha256','mask_sha256','mom0_bytes','mask_bytes']
    with (a.out_dir/'status.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=sf,extrasaction='ignore'); w.writeheader();
        for r in statuses: w.writerow(r)
    with (a.out_dir/'profile.csv').open('w',newline='') as f:
        pf=['galaxy_id','rad_hi_source_arcsec','sigma_hi_source_weight']; w=csv.DictWriter(f,fieldnames=pf); w.writeheader();
        for r in profiles: w.writerow(r)
    summary={'shard_index':a.shard_index,'shard_count':a.shard_count,'n_eligible':len(all_eligible),'n_pass':sum(r['status']=='pass' for r in statuses),
             'n_profile_rejected':sum(r['status']=='profile_rejected' for r in statuses),'n_product_or_geometry_failure':sum(r['status']=='product_or_geometry_failure' for r in statuses),
             'frozen_profiler_sha256':'eb290a85d5d4600a03602d5ab3a4c7c8ae1a7e40d18c74cfdee70f2b10804306','response_products_used':False}
    (a.out_dir/'summary.json').write_text(__import__('json').dumps(summary,indent=2,sort_keys=True)+'\n')
    print(__import__('json').dumps(summary,sort_keys=True))
if __name__=='__main__': main()
