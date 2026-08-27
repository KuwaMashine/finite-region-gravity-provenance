#!/usr/bin/env python3
import csv, importlib.util, json, math, re, shutil, subprocess
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
import astropy.units as u

ROOT=Path.cwd(); TMP=Path('/tmp/wallaby_nersc_validation'); OUT=ROOT/'post-stage0/source-optical-nersc-validation'
BASE=ROOT/'.github/scripts/source_optical_crossfield_canary.py'
spec=importlib.util.spec_from_file_location('base_canary',BASE); base=importlib.util.module_from_spec(spec); spec.loader.exec_module(base)
NAME='WALLABY J130150+041953'; RELEASE='NGC 4808 TR1'; BRICK='1953p042'; SUB='195'
TMP.mkdir(parents=True,exist_ok=True); OUT.mkdir(parents=True,exist_ok=True)
src=[r for r in csv.DictReader((ROOT/'post-stage0/source-only-census/phase2_source_ids.csv').open()) if r['name']==NAME and r['team_release']==RELEASE]
if len(src)!=1: raise SystemExit('frozen source row not unique')
src=src[0]; ra=float(src['ra']); dec=float(src['dec']); coord=SkyCoord(ra*u.deg,dec*u.deg)

# Direct official DR10 coadds.
full={}
for b in 'grz':
    url=f'https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr10/south/coadd/{SUB}/{BRICK}/legacysurvey-{BRICK}-image-{b}.fits.fz'
    p=TMP/f'brick_{b}.fits.fz'; base.curl_retry(url,p)
    with fits.open(p,memmap=False) as h:
        hh=h[1] if len(h)>1 and h[1].data is not None else h[0]
        full[b]=(np.asarray(hh.data).copy(),hh.header.copy())

dust=f'https://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust?locstr={ra}%20{dec}'
base.curl_retry(dust,TMP/'dust.xml'); ebv=base.dust_ebv(TMP/'dust.xml')

def center_from_aux(p):
    m=re.search(r'center x:\s*([0-9.+-]+) pix, y:\s*([0-9.+-]+) pix',Path(p).read_text(errors='replace'))
    if not m: raise SystemExit('missing fitted center')
    return tuple(map(float,m.groups()))

def run_size(size):
    work=TMP/f'n{size}'; work.mkdir(exist_ok=True); sout=OUT/f'n{size}'; sout.mkdir(exist_ok=True)
    for b in 'grz':
        data,hdr=full[b]; w=WCS(hdr).celestial; x,y=w.world_to_pixel(coord)
        cut=Cutout2D(data,(x,y),(size,size),wcs=w,mode='strict')
        arr=np.asarray(cut.data,dtype='float32')
        if arr.shape!=(size,size) or not np.isfinite(arr).all(): raise SystemExit(f'{b} bad crop size={size}')
        fits.PrimaryHDU(arr,header=cut.wcs.to_header()).writeto(work/f'{b}.fits',overwrite=True)
    slug=base.slugify(NAME)
    cfg=work/'r_config.py'; cfg.write_text(
        "ap_process_mode='image'\n"+f"ap_image_file=r'{work/'r.fits'}'\n"+f"ap_name='{slug}_r_{size}'\n"+
        f"ap_pixscale={base.PIX_SCALE}\nap_zeropoint={base.ZEROPOINT}\nap_doplot=False\nap_isoclip=True\n"+
        f"ap_guess_center={{'x': {size/2:.1f}, 'y': {size/2:.1f}}}\n")
    base.run_autoprof(cfg,work)
    rp=work/f'{slug}_r_{size}.prof'; ra_=work/f'{slug}_r_{size}.aux'; shutil.copy2(rp,work/'r.prof'); shutil.copy2(ra_,work/'r.aux')
    rl=work/f'{slug}_r_{size}.log';
    if rl.exists(): shutil.copy2(rl,work/'r.log')
    for b in ('g','z'):
        cfg=work/f'{b}_config.py'; cfg.write_text(
            "ap_process_mode='forced image'\n"+f"ap_image_file=r'{work/f'{b}.fits'}'\n"+f"ap_name='{slug}_{b}_{size}'\n"+
            f"ap_pixscale={base.PIX_SCALE}\nap_zeropoint={base.ZEROPOINT}\nap_doplot=False\nap_isoclip=True\n"+f"ap_forcing_profile=r'{work/'r.prof'}'\n")
        base.run_autoprof(cfg,work)
        pp=work/f'{slug}_{b}_{size}.prof'; aa=work/f'{slug}_{b}_{size}.aux'; shutil.copy2(pp,work/f'{b}.prof'); shutil.copy2(aa,work/f'{b}.aux')
        ll=work/f'{slug}_{b}_{size}.log';
        if ll.exists(): shutil.copy2(ll,work/f'{b}.log')
    details=base.full_stellar_from_profiles(src,work,sout,ebv)
    cx,cy=center_from_aux(work/'r.aux'); safe=(min(cx,(size-1)-cx,cy,(size-1)-cy)-2.0)*base.PIX_SCALE
    details['fitted_center_pix']={'x':cx,'y':cy}; details['orientation_independent_safe_radius_arcsec']=safe
    details['common_aperture_fully_contained']=details['common_aperture_radius_arcsec']<=safe
    for fn in ('r.prof','r.aux','r.log','g.prof','g.aux','g.log','z.prof','z.aux','z.log'):
        p=work/fn
        if p.exists(): shutil.copy2(p,sout/fn)
    return details

r512=run_size(512); r1024=run_size(1024)
old=json.loads((ROOT/'post-stage0/source-optical-crossfield-canary/crossfield_canary_summary.json').read_text())
viewer=next(r for r in old['results'] if r['source_name']==NAME)
report={
  'scope':'source-only direct NERSC DR10 photometry validation; no validation kinematics queried',
  'source_name':NAME,'team_release':RELEASE,'brick':BRICK,'backend':'portal.nersc.gov DR10 south coadd image-g/r/z',
  'nersc_512':r512,'nersc_1024':r1024,
  'viewer_512_reference':{'R50_r_arcsec':viewer['R50_r_arcsec'],'log10_Mstar_adopted_median':viewer['log10_Mstar_adopted_median'],'common_aperture_radius_arcsec':viewer['common_aperture_radius_arcsec']},
  'delta_nersc512_minus_viewer512':{'R50_r_arcsec':r512['R50_r_arcsec']-viewer['R50_r_arcsec'],'log10_Mstar':r512['log10_Mstar_adopted_median']-viewer['log10_Mstar_adopted_median']},
  'delta_nersc1024_minus_nersc512':{'R50_r_arcsec':r1024['R50_r_arcsec']-r512['R50_r_arcsec'],'log10_Mstar':r1024['log10_Mstar_adopted_median']-r512['log10_Mstar_adopted_median']},
}
(OUT/'nersc_validation_summary.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
files=sorted(p for p in OUT.rglob('*') if p.is_file() and p.name!='SHA256SUMS.txt')
with (OUT/'SHA256SUMS.txt').open('w') as f:
    for p in files: f.write(f'{base.sha256(p)}  {p.relative_to(ROOT)}\n')
print(json.dumps(report,indent=2,sort_keys=True))
