#!/usr/bin/env python3
import csv, importlib.util, json, math, re, shutil
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from dl import authClient as ac
from dl import queryClient as qc
from reproject import reproject_interp

ROOT=Path.cwd(); TMP=Path('/tmp/wallaby_nersc_multibrick'); OUT=ROOT/'post-stage0/source-optical-nersc-multibrick-probe'
BASE=ROOT/'.github/scripts/source_optical_crossfield_canary.py'
spec=importlib.util.spec_from_file_location('base_canary',BASE); base=importlib.util.module_from_spec(spec); spec.loader.exec_module(base)
NAME='WALLABY J133209-245132'; RELEASE='NGC 5044 TR3'; SIZE=1024; PIX=base.PIX_SCALE
TMP.mkdir(parents=True,exist_ok=True); OUT.mkdir(parents=True,exist_ok=True)
src=[r for r in csv.DictReader((ROOT/'post-stage0/source-only-census/phase2_source_ids.csv').open()) if r['name']==NAME and r['team_release']==RELEASE]
if len(src)!=1: raise SystemExit('frozen source row not unique')
src=src[0]; ra0=float(src['ra']); dec0=float(src['dec'])

# Source-centred TAN target grid at the frozen 0.262 arcsec/pixel scale.
tw=WCS(naxis=2)
tw.wcs.crpix=[(SIZE+1)/2.0,(SIZE+1)/2.0]
tw.wcs.cdelt=np.array([-PIX/3600.0,PIX/3600.0])
tw.wcs.crval=[ra0,dec0]
tw.wcs.ctype=['RA---TAN','DEC--TAN']
tw.wcs.cunit=['deg','deg']
# Target pixel centres for unique geometrical brick ownership.
yy,xx=np.indices((SIZE,SIZE),dtype=float)
tra,tdec=tw.pixel_to_world_values(xx,yy)
# Query a slight bounding pad around all target pixel centres.
ramin,ramax=float(np.nanmin(tra))-1e-4,float(np.nanmax(tra))+1e-4
decmin,decmax=float(np.nanmin(tdec))-1e-4,float(np.nanmax(tdec))+1e-4
if ramax-ramin>180: raise SystemExit('RA wrap not implemented for this probe')

token=ac.login('anonymous')
sql=("SELECT brickname,ra,dec,ra1,ra2,dec1,dec2 FROM ls_dr10.bricks "
     f"WHERE ra2>{ramin} AND ra1<{ramax} AND dec2>{decmin} AND dec1<{decmax} ORDER BY brickname")
body=qc.query(token=token,sql=sql,fmt='csv')
bricks=list(csv.DictReader(body.splitlines()))
if not bricks: raise SystemExit(f'no bricks returned for target box: {body[:500]!r}')
for b in bricks:
    for k in ('ra','dec','ra1','ra2','dec1','dec2'): b[k]=float(b[k])

# Precompute unique owner mask from published brick cells.
owners=[]
owner_count=np.zeros((SIZE,SIZE),dtype=np.int16)
for b in bricks:
    own=(tra>=b['ra1'])&(tra<b['ra2'])&(tdec>=b['dec1'])&(tdec<b['dec2'])
    owners.append(own); owner_count+=own.astype(np.int16)
if not np.all(owner_count==1):
    vals,counts=np.unique(owner_count,return_counts=True)
    raise SystemExit(f'target grid does not have unique brick ownership: {dict(zip(vals.tolist(),counts.tolist()))}')

backend=[]; mosaics={}
for band in 'grz':
    outarr=np.full((SIZE,SIZE),np.nan,dtype=np.float32)
    brec=[]
    for b,own in zip(bricks,owners):
        if not np.any(own): continue
        brick=b['brickname']; sub=brick[:3]
        url=f'https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr10/south/coadd/{sub}/{brick}/legacysurvey-{brick}-image-{band}.fits.fz'
        p=TMP/f'{brick}_{band}.fits.fz'; base.curl_retry(url,p)
        with fits.open(p,memmap=False) as h:
            hh=h[1] if len(h)>1 and h[1].data is not None else h[0]
            data=np.asarray(hh.data); sw=WCS(hh.header).celestial
        repro,foot=reproject_interp((data,sw),tw,shape_out=(SIZE,SIZE),order='bilinear',return_footprint=True)
        bad=own & ((foot<=0)|(~np.isfinite(repro)))
        if np.any(bad): raise SystemExit(f'{brick} {band}: {int(bad.sum())} owned pixels lack finite reprojection')
        outarr[own]=repro[own].astype(np.float32)
        brec.append({'brickname':brick,'url':url,'owned_target_pixels':int(own.sum()),'source_shape':list(data.shape),'sha256':base.sha256(p)})
    if not np.isfinite(outarr).all(): raise SystemExit(f'{band}: final mosaic has nonfinite pixels')
    fits.PrimaryHDU(outarr,header=tw.to_header()).writeto(TMP/f'{band}.fits',overwrite=True)
    mosaics[band]={'min':float(np.min(outarr)),'max':float(np.max(outarr)),'finite_fraction':float(np.isfinite(outarr).mean()),'bricks':brec,'mosaic_sha256':base.sha256(TMP/f'{band}.fits')}

# Frozen AutoProf r geometry and forced g/z.
slug=base.slugify(NAME)
cfg=TMP/'r_config.py'; cfg.write_text(
    "ap_process_mode='image'\n"+f"ap_image_file=r'{TMP/'r.fits'}'\n"+f"ap_name='{slug}_r_multibrick'\n"+
    f"ap_pixscale={PIX}\nap_zeropoint={base.ZEROPOINT}\nap_doplot=False\nap_isoclip=True\n"+
    f"ap_guess_center={{'x': {SIZE/2:.1f}, 'y': {SIZE/2:.1f}}}\n")
base.run_autoprof(cfg,TMP)
rp=TMP/f'{slug}_r_multibrick.prof'; ra_=TMP/f'{slug}_r_multibrick.aux'; shutil.copy2(rp,TMP/'r.prof'); shutil.copy2(ra_,TMP/'r.aux')
rl=TMP/f'{slug}_r_multibrick.log';
if rl.exists(): shutil.copy2(rl,TMP/'r.log')
for band in ('g','z'):
    cfg=TMP/f'{band}_config.py'; cfg.write_text(
        "ap_process_mode='forced image'\n"+f"ap_image_file=r'{TMP/f'{band}.fits'}'\n"+f"ap_name='{slug}_{band}_multibrick'\n"+
        f"ap_pixscale={PIX}\nap_zeropoint={base.ZEROPOINT}\nap_doplot=False\nap_isoclip=True\n"+f"ap_forcing_profile=r'{TMP/'r.prof'}'\n")
    base.run_autoprof(cfg,TMP)
    pp=TMP/f'{slug}_{band}_multibrick.prof'; aa=TMP/f'{slug}_{band}_multibrick.aux'; shutil.copy2(pp,TMP/f'{band}.prof'); shutil.copy2(aa,TMP/f'{band}.aux')
    ll=TMP/f'{slug}_{band}_multibrick.log';
    if ll.exists(): shutil.copy2(ll,TMP/f'{band}.log')

dust=f'https://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust?locstr={ra0}%20{dec0}'
base.curl_retry(dust,TMP/'dust.xml'); ebv=base.dust_ebv(TMP/'dust.xml')
details=base.full_stellar_from_profiles(src,TMP,OUT,ebv)
# Explicit containment.
txt=(TMP/'r.aux').read_text(errors='replace'); m=re.search(r'center x:\s*([0-9.+-]+) pix, y:\s*([0-9.+-]+) pix',txt)
if not m: raise SystemExit('cannot parse fitted centre')
cx,cy=map(float,m.groups()); safe=(min(cx,(SIZE-1)-cx,cy,(SIZE-1)-cy)-2.0)*PIX
details['fitted_center_pix']={'x':cx,'y':cy}; details['orientation_independent_safe_radius_arcsec']=safe; details['common_aperture_fully_contained']=details['common_aperture_radius_arcsec']<=safe
if not details['common_aperture_fully_contained']: raise SystemExit('multibrick photometry remains frame-limited')
legacy=json.loads((ROOT/'post-stage0/source-optical-frame-probe/frame_probe_summary.json').read_text())['details_1024']
report={
 'scope':'source-only direct NERSC DR10 geometrically-owned multi-brick 1024 validation; no validation kinematics queried',
 'source_name':NAME,'team_release':RELEASE,'size_pix':SIZE,'pixscale_arcsec':PIX,
 'target_wcs_header':dict(tw.to_header()),'target_bounds_deg':{'ra_min':ramin,'ra_max':ramax,'dec_min':decmin,'dec_max':decmax},
 'n_intersecting_bricks':len(bricks),'bricks':bricks,'owner_count_unique':True,'mosaics':mosaics,'details_nersc_multibrick':details,
 'legacy_viewer_1024_reference':{'R50_r_arcsec':legacy['R50_r_arcsec'],'Rd_star_kpc':legacy['Rd_star_kpc'],'log10_Mstar_adopted_median':legacy['log10_Mstar_adopted_median'],'common_aperture_radius_arcsec':legacy['common_aperture_radius_arcsec']},
 'delta_nersc_minus_legacy':{'R50_r_arcsec':details['R50_r_arcsec']-legacy['R50_r_arcsec'],'Rd_star_kpc':details['Rd_star_kpc']-legacy['Rd_star_kpc'],'log10_Mstar':details['log10_Mstar_adopted_median']-legacy['log10_Mstar_adopted_median'],'common_aperture_radius_arcsec':details['common_aperture_radius_arcsec']-legacy['common_aperture_radius_arcsec']},
 'mosaic_rule':'Each target pixel is assigned only from the single published ls_dr10.bricks geometric cell containing that target-pixel sky coordinate; neighboring coadds are bilinearly reprojected to a source-centred TAN target grid before assignment; no overlap averaging.'
}
(OUT/'nersc_multibrick_probe_summary.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
for fn in ('r.prof','r.aux','r.log','g.prof','g.aux','g.log','z.prof','z.aux','z.log'):
    p=TMP/fn
    if p.exists(): shutil.copy2(p,OUT/fn)
files=sorted(p for p in OUT.rglob('*') if p.is_file() and p.name!='SHA256SUMS.txt')
with (OUT/'SHA256SUMS.txt').open('w') as f:
    for p in files: f.write(f'{base.sha256(p)}  {p.relative_to(ROOT)}\n')
print(json.dumps(report,indent=2,sort_keys=True))
