#!/usr/bin/env python3
import csv, importlib.util, json, math, re, shutil
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
import astropy.units as u

ROOT=Path.cwd(); TMP=Path('/tmp/wallaby_nersc_frame_convergence'); OUT=ROOT/'post-stage0/source-optical-nersc-frame-convergence'
BASE=ROOT/'.github/scripts/source_optical_crossfield_canary.py'
spec=importlib.util.spec_from_file_location('base_canary',BASE); base=importlib.util.module_from_spec(spec); spec.loader.exec_module(base)
NAME='WALLABY J130150+041953'; RELEASE='NGC 4808 TR1'; BRICK='1953p042'; SUB='195'; SIZES=[512,768,1024,1280]
TMP.mkdir(parents=True,exist_ok=True); OUT.mkdir(parents=True,exist_ok=True)
src=[r for r in csv.DictReader((ROOT/'post-stage0/source-only-census/phase2_source_ids.csv').open()) if r['name']==NAME and r['team_release']==RELEASE]
if len(src)!=1: raise SystemExit('frozen source row not unique')
src=src[0]; ra=float(src['ra']); dec=float(src['dec']); coord=SkyCoord(ra*u.deg,dec*u.deg)

full={}
for b in 'grz':
    url=f'https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr10/south/coadd/{SUB}/{BRICK}/legacysurvey-{BRICK}-image-{b}.fits.fz'
    p=TMP/f'brick_{b}.fits.fz'; base.curl_retry(url,p)
    with fits.open(p,memmap=False) as h:
        hh=h[1] if len(h)>1 and h[1].data is not None else h[0]
        full[b]=(np.asarray(hh.data).copy(),hh.header.copy())

dust=f'https://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust?locstr={ra}%20{dec}'
base.curl_retry(dust,TMP/'dust.xml'); ebv=base.dust_ebv(TMP/'dust.xml')

def parse_aux(path):
    txt=Path(path).read_text(errors='replace')
    def one(pattern):
        m=re.search(pattern,txt); return None if not m else float(m.group(1))
    cm=re.search(r'center x:\s*([0-9.+-]+) pix, y:\s*([0-9.+-]+) pix',txt)
    return {
        'center_x_pix':None if not cm else float(cm.group(1)),
        'center_y_pix':None if not cm else float(cm.group(2)),
        'background_flux_per_pix':one(r'background:\s*([-+0-9.eE]+)'),
        'background_noise_flux_per_pix':one(r'noise:\s*([-+0-9.eE]+) flux/pix'),
        'fit_limit_semimajor_pix':one(r'fit limit semi-major axis:\s*([-+0-9.eE]+) pix'),
        'global_ellipticity':one(r'global ellipticity:\s*([-+0-9.eE]+)'),
        'size_pix':one(r'size:\s*([-+0-9.eE]+) pix'),
    }

def run_size(size):
    work=TMP/f'n{size}'; sout=OUT/f'n{size}'; work.mkdir(exist_ok=True); sout.mkdir(exist_ok=True)
    for b in 'grz':
        data,hdr=full[b]; w=WCS(hdr).celestial; x,y=w.world_to_pixel(coord)
        half=size/2.0
        margins=[x,(data.shape[1]-1)-x,y,(data.shape[0]-1)-y]
        if min(margins)<half: raise SystemExit(f'size {size} does not fit single brick; margins={margins}')
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
    aux=parse_aux(work/'r.aux'); cx=aux['center_x_pix']; cy=aux['center_y_pix']
    safe=(min(cx,(size-1)-cx,cy,(size-1)-cy)-2.0)*base.PIX_SCALE
    details.update({'cutout_size_pix':size,'cutout_side_arcsec':size*base.PIX_SCALE,'r_aux':aux,
                    'orientation_independent_safe_radius_arcsec':safe,
                    'common_aperture_fully_contained':details['common_aperture_radius_arcsec']<=safe})
    for fn in ('r.prof','r.aux','r.log','g.prof','g.aux','g.log','z.prof','z.aux','z.log'):
        p=work/fn
        if p.exists(): shutil.copy2(p,sout/fn)
    return details

results=[run_size(s) for s in SIZES]
steps=[]
for a,b in zip(results[:-1],results[1:]):
    steps.append({
      'from_size_pix':a['cutout_size_pix'],'to_size_pix':b['cutout_size_pix'],
      'delta_R50_r_arcsec':b['R50_r_arcsec']-a['R50_r_arcsec'],
      'delta_log10_Mstar':b['log10_Mstar_adopted_median']-a['log10_Mstar_adopted_median'],
      'delta_common_aperture_radius_arcsec':b['common_aperture_radius_arcsec']-a['common_aperture_radius_arcsec'],
      'delta_background_flux_per_pix':b['r_aux']['background_flux_per_pix']-a['r_aux']['background_flux_per_pix'],
      'delta_background_noise_flux_per_pix':b['r_aux']['background_noise_flux_per_pix']-a['r_aux']['background_noise_flux_per_pix'],
    })
report={
 'scope':'source-only direct NERSC DR10 finite-frame convergence ladder; no validation kinematics queried',
 'source_name':NAME,'team_release':RELEASE,'brick':BRICK,'sizes_pix':SIZES,
 'results':results,'successive_deltas':steps,
 'decision_rule':'choose the smallest tested fully-contained frame for which the next larger tested frame changes log10_Mstar by <=0.02 dex and R50 by <=0.5 arcsec; if no tested frame passes, do not freeze a production size from this ladder.'
}
# Evaluate rule prospectively from the ladder only.
chosen=None
for i in range(len(results)-1):
    r=results[i]; d=steps[i]
    if r['common_aperture_fully_contained'] and abs(d['delta_log10_Mstar'])<=0.02 and abs(d['delta_R50_r_arcsec'])<=0.5:
        chosen=r['cutout_size_pix']; break
report['smallest_converged_size_pix']=chosen
(OUT/'nersc_frame_convergence_summary.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
files=sorted(p for p in OUT.rglob('*') if p.is_file() and p.name!='SHA256SUMS.txt')
with (OUT/'SHA256SUMS.txt').open('w') as f:
    for p in files: f.write(f'{base.sha256(p)}  {p.relative_to(ROOT)}\n')
print(json.dumps(report,indent=2,sort_keys=True))
