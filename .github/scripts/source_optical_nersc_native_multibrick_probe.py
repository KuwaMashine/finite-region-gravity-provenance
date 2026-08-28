#!/usr/bin/env python3
import csv, copy, importlib.util, json, math, re, shutil
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from dl import authClient as ac
from dl import queryClient as qc
from reproject import reproject_interp

ROOT=Path.cwd(); TMP=Path('/tmp/wallaby_nersc_native_multibrick'); OUT=ROOT/'post-stage0/source-optical-nersc-native-multibrick-probe'
BASE=ROOT/'.github/scripts/source_optical_crossfield_canary.py'
spec=importlib.util.spec_from_file_location('base_canary',BASE); base=importlib.util.module_from_spec(spec); spec.loader.exec_module(base)
NAME='WALLABY J133209-245132'; RELEASE='NGC 5044 TR3'; SIZE=1024; PIX=base.PIX_SCALE
TMP.mkdir(parents=True,exist_ok=True); OUT.mkdir(parents=True,exist_ok=True)
src=[r for r in csv.DictReader((ROOT/'post-stage0/source-only-census/phase2_source_ids.csv').open()) if r['name']==NAME and r['team_release']==RELEASE]
if len(src)!=1: raise SystemExit('frozen source row not unique')
src=src[0]; ra0=float(src['ra']); dec0=float(src['dec']); coord=SkyCoord(ra0*u.deg,dec0*u.deg)
token=ac.login('anonymous')

# Resolve the unique geometrical brick containing the source.
sql=f"SELECT brickname,ra,dec,ra1,ra2,dec1,dec2 FROM ls_dr10.bricks WHERE {ra0}>=ra1 AND {ra0}<ra2 AND {dec0}>=dec1 AND {dec0}<dec2"
rr=list(csv.DictReader(qc.query(token=token,sql=sql,fmt='csv').splitlines()))
if len(rr)!=1: raise SystemExit(f'expected one source brick, got {len(rr)}')
primary={k:(float(v) if k!='brickname' else v.strip()) for k,v in rr[0].items()}
pbrick=primary['brickname']; psub=pbrick[:3]

# Load primary r brick to define a target grid that is an exact integer-pixel extension of its native WCS.
purl=f'https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr10/south/coadd/{psub}/{pbrick}/legacysurvey-{pbrick}-image-r.fits.fz'
pr=TMP/f'{pbrick}_r.fits.fz'; base.curl_retry(purl,pr)
with fits.open(pr,memmap=False) as h:
    hh=h[1] if len(h)>1 and h[1].data is not None else h[0]
    pdata_r=np.asarray(hh.data).copy(); pw=WCS(hh.header).celestial
sx,sy=pw.world_to_pixel(coord)
sx=float(np.asarray(sx).reshape(-1)[0]); sy=float(np.asarray(sy).reshape(-1)[0])
x0=int(round(sx))-SIZE//2; y0=int(round(sy))-SIZE//2
tw=copy.deepcopy(pw); tw.wcs.crpix[0]-=x0; tw.wcs.crpix[1]-=y0
target_source_x,target_source_y=tw.world_to_pixel(coord)
target_source_x=float(np.asarray(target_source_x).reshape(-1)[0]); target_source_y=float(np.asarray(target_source_y).reshape(-1)[0])

# Determine all bricks intersecting the actual native-grid target footprint.
yy,xx=np.indices((SIZE,SIZE),dtype=float); tra,tdec=tw.pixel_to_world_values(xx,yy)
ramin,ramax=float(np.nanmin(tra))-1e-4,float(np.nanmax(tra))+1e-4
decmin,decmax=float(np.nanmin(tdec))-1e-4,float(np.nanmax(tdec))+1e-4
if ramax-ramin>180: raise SystemExit('RA wrap not implemented for this probe')
sql=("SELECT brickname,ra,dec,ra1,ra2,dec1,dec2 FROM ls_dr10.bricks "
     f"WHERE ra2>{ramin} AND ra1<{ramax} AND dec2>{decmin} AND dec1<{decmax} ORDER BY brickname")
bricks=list(csv.DictReader(qc.query(token=token,sql=sql,fmt='csv').splitlines()))
if not bricks: raise SystemExit('no intersecting bricks')
for b in bricks:
    for k in ('ra','dec','ra1','ra2','dec1','dec2'): b[k]=float(b[k])
owners=[]; owner_count=np.zeros((SIZE,SIZE),dtype=np.int16)
for b in bricks:
    own=(tra>=b['ra1'])&(tra<b['ra2'])&(tdec>=b['dec1'])&(tdec<b['dec2'])
    owners.append(own); owner_count+=own.astype(np.int16)
if not np.all(owner_count==1):
    vals,cnts=np.unique(owner_count,return_counts=True); raise SystemExit(f'nonunique owner map: {dict(zip(vals.tolist(),cnts.tolist()))}')

mosaics={}
for band in 'grz':
    outarr=np.full((SIZE,SIZE),np.nan,dtype=np.float32); brec=[]
    for b,own in zip(bricks,owners):
        if not np.any(own): continue
        brick=b['brickname']; sub=brick[:3]
        p=TMP/f'{brick}_{band}.fits.fz'
        url=f'https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr10/south/coadd/{sub}/{brick}/legacysurvey-{brick}-image-{band}.fits.fz'
        if not p.exists(): base.curl_retry(url,p)
        with fits.open(p,memmap=False) as h:
            hh=h[1] if len(h)>1 and h[1].data is not None else h[0]
            data=np.asarray(hh.data); sw=WCS(hh.header).celestial
        if brick==pbrick:
            # Exact native assignment: target pixel (x,y) is primary pixel (x+x0,y+y0).
            srcx=xx.astype(np.int64)+x0; srcy=yy.astype(np.int64)+y0
            valid=own & (srcx>=0)&(srcx<data.shape[1])&(srcy>=0)&(srcy<data.shape[0])
            if int(valid.sum())!=int(own.sum()): raise SystemExit(f'{band}: primary owner includes pixels outside primary array')
            outarr[own]=data[srcy[own],srcx[own]].astype(np.float32)
            mode='exact_native_integer_pixels'
        else:
            repro,foot=reproject_interp((data,sw),tw,shape_out=(SIZE,SIZE),order='bilinear',return_footprint=True)
            bad=own & ((foot<=0)|(~np.isfinite(repro)))
            if np.any(bad): raise SystemExit(f'{brick} {band}: {int(bad.sum())} owned pixels lack finite reprojection')
            outarr[own]=repro[own].astype(np.float32); mode='bilinear_neighbor_only'
        brec.append({'brickname':brick,'mode':mode,'owned_target_pixels':int(own.sum()),'url':url,'sha256':base.sha256(p)})
    if not np.isfinite(outarr).all(): raise SystemExit(f'{band}: final mosaic nonfinite')
    fits.PrimaryHDU(outarr,header=tw.to_header()).writeto(TMP/f'{band}.fits',overwrite=True)
    mosaics[band]={'finite_fraction':float(np.isfinite(outarr).mean()),'min':float(np.min(outarr)),'max':float(np.max(outarr)),'mosaic_sha256':base.sha256(TMP/f'{band}.fits'),'bricks':brec}

# Frozen AutoProf, with source coordinate as the deterministic initial centre on this native grid.
slug=base.slugify(NAME)
cfg=TMP/'r_config.py'; cfg.write_text(
    "ap_process_mode='image'\n"+f"ap_image_file=r'{TMP/'r.fits'}'\n"+f"ap_name='{slug}_r_native_multibrick'\n"+
    f"ap_pixscale={PIX}\nap_zeropoint={base.ZEROPOINT}\nap_doplot=False\nap_isoclip=True\n"+
    f"ap_guess_center={{'x': {target_source_x!r}, 'y': {target_source_y!r}}}\n")
base.run_autoprof(cfg,TMP)
rp=TMP/f'{slug}_r_native_multibrick.prof'; ra_=TMP/f'{slug}_r_native_multibrick.aux'; shutil.copy2(rp,TMP/'r.prof'); shutil.copy2(ra_,TMP/'r.aux')
rl=TMP/f'{slug}_r_native_multibrick.log';
if rl.exists(): shutil.copy2(rl,TMP/'r.log')
for band in ('g','z'):
    cfg=TMP/f'{band}_config.py'; cfg.write_text(
        "ap_process_mode='forced image'\n"+f"ap_image_file=r'{TMP/f'{band}.fits'}'\n"+f"ap_name='{slug}_{band}_native_multibrick'\n"+
        f"ap_pixscale={PIX}\nap_zeropoint={base.ZEROPOINT}\nap_doplot=False\nap_isoclip=True\n"+f"ap_forcing_profile=r'{TMP/'r.prof'}'\n")
    base.run_autoprof(cfg,TMP)
    pp=TMP/f'{slug}_{band}_native_multibrick.prof'; aa=TMP/f'{slug}_{band}_native_multibrick.aux'; shutil.copy2(pp,TMP/f'{band}.prof'); shutil.copy2(aa,TMP/f'{band}.aux')
    ll=TMP/f'{slug}_{band}_native_multibrick.log';
    if ll.exists(): shutil.copy2(ll,TMP/f'{band}.log')

dust=f'https://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust?locstr={ra0}%20{dec0}'
base.curl_retry(dust,TMP/'dust.xml'); ebv=base.dust_ebv(TMP/'dust.xml')
details=base.full_stellar_from_profiles(src,TMP,OUT,ebv)
txt=(TMP/'r.aux').read_text(errors='replace'); m=re.search(r'center x:\s*([0-9.+-]+) pix, y:\s*([0-9.+-]+) pix',txt)
if not m: raise SystemExit('cannot parse fitted centre')
cx,cy=map(float,m.groups()); safe=(min(cx,(SIZE-1)-cx,cy,(SIZE-1)-cy)-2.0)*PIX
details['fitted_center_pix']={'x':cx,'y':cy}; details['orientation_independent_safe_radius_arcsec']=safe; details['common_aperture_fully_contained']=details['common_aperture_radius_arcsec']<=safe
if not details['common_aperture_fully_contained']: raise SystemExit('native multibrick photometry remains frame-limited')
legacy=json.loads((ROOT/'post-stage0/source-optical-frame-probe/frame_probe_summary.json').read_text())['details_1024']
bilinear=json.loads((ROOT/'post-stage0/source-optical-nersc-multibrick-probe/nersc_multibrick_probe_summary.json').read_text())['details_nersc_multibrick']
report={
 'scope':'source-only direct NERSC DR10 native-primary-grid multi-brick 1024 validation; no validation kinematics queried',
 'source_name':NAME,'team_release':RELEASE,'size_pix':SIZE,'pixscale_arcsec':PIX,
 'primary_brick':primary,'primary_source_pixel_xy':[sx,sy],'target_integer_origin_in_primary_xy':[x0,y0],
 'source_pixel_in_target_xy':[target_source_x,target_source_y],'n_intersecting_bricks':len(bricks),'bricks':bricks,
 'target_wcs_header':dict(tw.to_header()),'owner_count_unique':True,'mosaics':mosaics,'details_native_multibrick':details,
 'legacy_viewer_1024_reference':{'R50_r_arcsec':legacy['R50_r_arcsec'],'Rd_star_kpc':legacy['Rd_star_kpc'],'log10_Mstar_adopted_median':legacy['log10_Mstar_adopted_median'],'common_aperture_radius_arcsec':legacy['common_aperture_radius_arcsec']},
 'source_centered_bilinear_nersc_reference':{'R50_r_arcsec':bilinear['R50_r_arcsec'],'Rd_star_kpc':bilinear['Rd_star_kpc'],'log10_Mstar_adopted_median':bilinear['log10_Mstar_adopted_median'],'common_aperture_radius_arcsec':bilinear['common_aperture_radius_arcsec']},
 'delta_native_minus_legacy':{'R50_r_arcsec':details['R50_r_arcsec']-legacy['R50_r_arcsec'],'Rd_star_kpc':details['Rd_star_kpc']-legacy['Rd_star_kpc'],'log10_Mstar':details['log10_Mstar_adopted_median']-legacy['log10_Mstar_adopted_median'],'common_aperture_radius_arcsec':details['common_aperture_radius_arcsec']-legacy['common_aperture_radius_arcsec']},
 'delta_native_minus_bilinear_nersc':{'R50_r_arcsec':details['R50_r_arcsec']-bilinear['R50_r_arcsec'],'Rd_star_kpc':details['Rd_star_kpc']-bilinear['Rd_star_kpc'],'log10_Mstar':details['log10_Mstar_adopted_median']-bilinear['log10_Mstar_adopted_median']},
 'mosaic_rule':'Target WCS is an exact integer-pixel extension of the containing DR10 brick WCS. Pixels geometrically owned by the containing brick are copied without interpolation; only neighboring-brick owned pixels are bilinearly reprojected. No overlap averaging.'
}
(OUT/'nersc_native_multibrick_probe_summary.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
for fn in ('r.prof','r.aux','r.log','g.prof','g.aux','g.log','z.prof','z.aux','z.log'):
    p=TMP/fn
    if p.exists(): shutil.copy2(p,OUT/fn)
files=sorted(p for p in OUT.rglob('*') if p.is_file() and p.name!='SHA256SUMS.txt')
with (OUT/'SHA256SUMS.txt').open('w') as f:
    for p in files: f.write(f'{base.sha256(p)}  {p.relative_to(ROOT)}\n')
print(json.dumps(report,indent=2,sort_keys=True))
