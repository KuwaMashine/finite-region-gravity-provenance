#!/usr/bin/env python3
import csv
import importlib.util
import json
import math
import re
import shutil
from pathlib import Path

import numpy as np
from astropy.io import fits

ROOT=Path.cwd()
BASE=ROOT/'.github/scripts/source_optical_crossfield_canary.py'
spec=importlib.util.spec_from_file_location('canary_base',BASE)
base=importlib.util.module_from_spec(spec); spec.loader.exec_module(base)

NAME='WALLABY J133209-245132'
RELEASE='NGC 5044 TR3'
SIZE=1024
TMP=Path('/tmp/wallaby_frame_probe')
OUT=ROOT/'post-stage0/source-optical-frame-probe'
TMP.mkdir(parents=True,exist_ok=True); OUT.mkdir(parents=True,exist_ok=True)

rows=[r for r in csv.DictReader((ROOT/'post-stage0/source-only-census/phase2_source_ids.csv').open()) if r['name']==NAME and r['team_release']==RELEASE]
if len(rows)!=1: raise SystemExit(f'expected one frozen source row, got {len(rows)}')
src=rows[0]; ra=float(src['ra']); dec=float(src['dec'])

cut=f'https://www.legacysurvey.org/viewer/cutout.fits?ra={ra}&dec={dec}&layer=ls-dr10&pixscale={base.PIX_SCALE}&bands=grz&size={SIZE}'
dust=f'https://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust?locstr={ra}%20{dec}'
base.curl_retry(cut,TMP/'dr10_grz.fits'); base.curl_retry(dust,TMP/'dust.xml')
with fits.open(TMP/'dr10_grz.fits',memmap=False) as h:
    d=np.asarray(h[0].data); bands=str(h[0].header.get('BANDS',''))
if d.ndim!=3 or d.shape!=(3,SIZE,SIZE) or bands!='grz': raise SystemExit(f'unexpected cube {d.shape} bands={bands!r}')
if not np.isfinite(d).all(): raise SystemExit('nonfinite DR10 pixels')
for i,b in enumerate('grz'): fits.PrimaryHDU(d[i].astype('float32')).writeto(TMP/f'{b}.fits',overwrite=True)

slug=base.slugify(NAME)
r_cfg=TMP/'r_config.py'
r_cfg.write_text(
    "ap_process_mode='image'\n"
    f"ap_image_file=r'{TMP/'r.fits'}'\n"
    f"ap_name='{slug}_r_1024'\n"
    f"ap_pixscale={base.PIX_SCALE}\n"
    f"ap_zeropoint={base.ZEROPOINT}\n"
    "ap_doplot=False\n"
    "ap_isoclip=True\n"
    f"ap_guess_center={{'x': {SIZE/2:.1f}, 'y': {SIZE/2:.1f}}}\n"
)
base.run_autoprof(r_cfg,TMP)
rprof=TMP/f'{slug}_r_1024.prof'; raux=TMP/f'{slug}_r_1024.aux'
if not rprof.exists() or not raux.exists(): raise SystemExit('missing 1024 r outputs')
shutil.copy2(rprof,TMP/'r.prof'); shutil.copy2(raux,TMP/'r.aux')
rlog=TMP/f'{slug}_r_1024.log'
if rlog.exists(): shutil.copy2(rlog,TMP/'r.log')

for b in ('g','z'):
    cfg=TMP/f'{b}_config.py'
    cfg.write_text(
        "ap_process_mode='forced image'\n"
        f"ap_image_file=r'{TMP/f'{b}.fits'}'\n"
        f"ap_name='{slug}_{b}_1024'\n"
        f"ap_pixscale={base.PIX_SCALE}\n"
        f"ap_zeropoint={base.ZEROPOINT}\n"
        "ap_doplot=False\n"
        "ap_isoclip=True\n"
        f"ap_forcing_profile=r'{TMP/'r.prof'}'\n"
    )
    base.run_autoprof(cfg,TMP)
    prof=TMP/f'{slug}_{b}_1024.prof'; aux=TMP/f'{slug}_{b}_1024.aux'
    if not prof.exists() or not aux.exists(): raise SystemExit(f'missing forced {b} outputs')
    shutil.copy2(prof,TMP/f'{b}.prof'); shutil.copy2(aux,TMP/f'{b}.aux')
    log=TMP/f'{slug}_{b}_1024.log'
    if log.exists(): shutil.copy2(log,TMP/f'{b}.log')

ebv=base.dust_ebv(TMP/'dust.xml')
details=base.full_stellar_from_profiles(src,TMP,OUT,ebv)

m=re.search(r'center x:\s*([0-9.+-]+) pix, y:\s*([0-9.+-]+) pix',(TMP/'r.aux').read_text())
if not m: raise SystemExit('unable to read r-band fitted center')
cx,cy=map(float,m.groups())
safe_pix=min(cx,(SIZE-1)-cx,cy,(SIZE-1)-cy)
safe_arcsec=safe_pix*base.PIX_SCALE
contained=details['common_aperture_radius_arcsec'] <= safe_arcsec

old=json.loads((ROOT/'post-stage0/source-optical-crossfield-canary/crossfield_canary_summary.json').read_text())
oldrow=next(r for r in old['results'] if r['source_name']==NAME)
report={
    'scope':'source-only 1024-pixel frame-containment correction probe; no validation kinematics queried',
    'source_name':NAME,'team_release':RELEASE,'cutout_size_pix':SIZE,'pixscale_arcsec':base.PIX_SCALE,
    'cutout_side_arcsec':SIZE*base.PIX_SCALE,'fitted_center_pix':{'x':cx,'y':cy},
    'orientation_independent_safe_radius_arcsec':safe_arcsec,
    'common_aperture_fully_contained':contained,
    'full_stellar_1024':bool(contained and details['n_mlcr_estimates']==45 and details['n_retained']>=5),
    'details_1024':details,
    'comparison_512':{
        'common_aperture_radius_arcsec':oldrow['common_aperture_radius_arcsec'],
        'R50_r_arcsec':oldrow['R50_r_arcsec'],'Rd_star_kpc':oldrow['Rd_star_kpc'],
        'log10_Mstar_adopted_median':oldrow['log10_Mstar_adopted_median'],
        'orientation_independent_safe_radius_arcsec_from_aux':63.09222,
        'contained':False,
    },
    'delta_1024_minus_512':{
        'R50_r_arcsec':details['R50_r_arcsec']-oldrow['R50_r_arcsec'],
        'Rd_star_kpc':details['Rd_star_kpc']-oldrow['Rd_star_kpc'],
        'log10_Mstar_adopted_median':details['log10_Mstar_adopted_median']-oldrow['log10_Mstar_adopted_median'],
    },
    'decision':'512-pixel mass receipt is invalid for production because its adopted common aperture was not guaranteed contained; 1024 result is usable only if common_aperture_fully_contained=true.'
}
(OUT/'frame_probe_summary.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
for fn in ('r.prof','r.aux','r.log','g.prof','g.aux','g.log','z.prof','z.aux','z.log'):
    p=TMP/fn
    if p.exists(): shutil.copy2(p,OUT/fn)
files=sorted(p for p in OUT.iterdir() if p.is_file() and p.name!='SHA256SUMS.txt')
with (OUT/'SHA256SUMS.txt').open('w') as f:
    for p in files: f.write(f'{base.sha256(p)}  {p.relative_to(ROOT)}\n')
print(json.dumps(report,indent=2,sort_keys=True))
