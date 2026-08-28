#!/usr/bin/env python3
"""One-object TARGET preflight through the exact frozen v4.99 Stage-1 source engine.

The preflight is intentionally source-only. It obtains only the CADC source_data
plane for one NGC 5044 TR3 source, constructs the frozen H I source profile,
normalizes an already validated native-grid DR10 optical profile, writes a
SOURCE_ONLY export manifest, then invokes the byte-verified v4.99 TARGET runner.
No kinematic plane/catalogue or response column is queried or read.
"""
from __future__ import annotations
import csv, hashlib, json, math, os, subprocess, sys, urllib.parse, urllib.request
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.wcs.utils import proj_plane_pixel_scales

ROOT=Path.cwd()
OUT=ROOT/'post-stage0/source-stage-target-preflight'
TMP=Path('/tmp/wallaby_source_stage_target_preflight')
ENG=ROOT/'post-stage0/frozen-v4.99-source-engine'
NAME='WALLABY J133209-245132'
OBS='WALLABY_J133209-245132'
TEAM='NGC 5044 TR3'
FIELD='NGC5044'
PLANE='source_data_NGC5044_TR3'
BASE='https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/data/pub/WALLABY'
TAP='https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/argus/sync'
OBJ='WALLABY_J133209-245132_NGC5044_TR3'
OPT=ROOT/'post-stage0/source-optical-nersc-native-multibrick-probe'
RECEIPT=ROOT/'post-stage0/STAGE0_EXTERNAL_TIMESTAMP_RECEIPT.json'
FORBIDDEN=('kinematic','kinmodel','rotcur','model_data','vrot','gobs','gbar','rar_residual','vmax','vflat')

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def download(url,p):
    req=urllib.request.Request(url,headers={'User-Agent':'finite-region-gravity-source-preflight/1'})
    with urllib.request.urlopen(req,timeout=180) as r, Path(p).open('wb') as f:
        while True:
            b=r.read(1<<20)
            if not b: break
            f.write(b)
    if Path(p).stat().st_size==0: raise RuntimeError(f'zero-byte download {url}')

def tap_source_artifacts():
    q=("SELECT o.observationID,p.productID,a.uri,a.productType,a.contentType,a.contentLength,a.contentChecksum "
       "FROM caom2.Observation AS o JOIN caom2.Plane AS p ON o.obsID=p.obsID "
       "JOIN caom2.Artifact AS a ON p.planeID=a.planeID "
       f"WHERE o.collection='WALLABY' AND o.observationID='{OBS}' AND p.productID='{PLANE}'")
    data=urllib.parse.urlencode({'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'csv','QUERY':q}).encode()
    req=urllib.request.Request(TAP,data=data,headers={'User-Agent':'finite-region-gravity-source-preflight/1'})
    with urllib.request.urlopen(req,timeout=180) as r: text=r.read().decode('utf-8-sig')
    rows=list(csv.DictReader(text.splitlines()))
    if not rows: raise RuntimeError('CADC source_data plane returned no artifacts')
    bad=[r for r in rows if any(x in (r.get('uri','')+' '+r.get('productID','')).lower() for x in FORBIDDEN)]
    if bad: raise RuntimeError('source-plane artifact firewall violation')
    return rows,text,q

def source_row():
    df=pd.read_csv(ROOT/'post-stage0/source-only-probe/source_catalogue_safe.csv')
    z=df[(df['name']==NAME)&(df['team_release']==TEAM)]
    if len(z)!=1: raise RuntimeError(f'expected one final source row, got {len(z)}')
    r=z.iloc[0]
    if float(r.qflag)!=0.0: raise RuntimeError('preflight source qflag is not zero')
    for c in ('ra','dec','dist_h','log_m_hi_corr','ell_maj','ell_min','ell_pa'):
        if not math.isfinite(float(r[c])): raise RuntimeError(f'nonfinite source field {c}')
    return r

def project_source_pixels(mom,mask,header,r):
    m=np.asarray(mom).squeeze(); k=np.asarray(mask).squeeze()
    if m.ndim!=2 or k.ndim!=3 or tuple(k.shape[-2:])!=tuple(m.shape):
        raise RuntimeError(f'unexpected moment0/mask shapes {m.shape} {k.shape}')
    mask2=np.any(k>0,axis=0); yy,xx=np.nonzero(mask2 & np.isfinite(m))
    if len(xx)<8: raise RuntimeError('fewer than eight finite moment0 pixels in source footprint')
    w=WCS(header).celestial
    sc=w.pixel_to_world(xx.astype(float),yy.astype(float))
    center=SkyCoord(float(r.ra)*u.deg,float(r.dec)*u.deg,frame='icrs')
    # spherical_offsets_to gives (d_lon, d_lat): east and north offsets.
    de,dn=center.spherical_offsets_to(sc)
    pix=pd.DataFrame({'galaxy_id':NAME,'dx_east_arcsec':de.to_value(u.arcsec),
                      'dy_north_arcsec':dn.to_value(u.arcsec),'moment0_weight':m[yy,xx].astype(float)})
    scales=np.abs(np.asarray(proj_plane_pixel_scales(w),float))*3600.0
    ps=float(np.mean(scales))
    geo=pd.DataFrame([{'galaxy_id':NAME,'ell_maj_arcsec':6.0*float(r.ell_maj),
                       'ell_min_arcsec':6.0*float(r.ell_min),'ell_pa_deg':float(r.ell_pa),
                       'pixel_scale_arcsec':ps}])
    return pix,geo,{'moment0_shape':list(m.shape),'mask_shape':list(k.shape),'projected_mask_finite_pixels':int(len(xx)),
                    'pixel_scale_arcsec_xy':[float(x) for x in scales],'pixel_scale_arcsec_mean':ps}

def read_prof(p):
    with Path(p).open() as f:
        f.readline(); return list(csv.DictReader(f))

def normalized_optical():
    P={b:read_prof(OPT/f'{b}.prof') for b in 'grz'}
    n=min(len(P[b]) for b in 'grz'); keep=n
    for i in range(n):
        rs=[float(P[b][i]['R']) for b in 'grz']
        if max(rs)-min(rs)>1e-8: raise RuntimeError(f'optical radius mismatch row {i}')
        es=[float(P[b][i]['SB_e']) for b in 'grz']
        if any((not math.isfinite(x)) or x>=0.22 for x in es): keep=i; break
    if keep<5: raise RuntimeError('fewer than five common optical rows')
    rows=[]
    for z in P['r'][:keep]:
        rows.append({'galaxy_id':NAME,'radius_arcsec':float(z['R']),'mu_r':float(z['SB']),
                     'muerr_r':float(z['SB_e']),'ellipticity':float(z['ellip'])})
    return pd.DataFrame(rows),{'n_common_retained':keep,'common_aperture_radius_arcsec':float(P['r'][keep-1]['R'])}

def main():
    OUT.mkdir(parents=True,exist_ok=True); TMP.mkdir(parents=True,exist_ok=True)
    r=source_row(); arts,art_csv,tap_query=tap_source_artifacts()
    (OUT/'cadc_source_artifacts.csv').write_text(art_csv)
    wanted={}
    for kind,suffix in [('mom0','_mom0.fits'),('mask','_mask.fits')]:
        matches=[a for a in arts if a.get('uri','').endswith(suffix)]
        if len(matches)!=1: raise RuntimeError(f'expected one {kind} artifact, got {len(matches)}')
        uri=matches[0]['uri']; fn=uri.split('/')[-1]
        if not fn.startswith(OBJ): raise RuntimeError(f'unexpected source artifact filename {fn}')
        p=TMP/f'{kind}.fits'; download(f'{BASE}/{fn}',p); wanted[kind]={'uri':uri,'filename':fn,'sha256':sha(p),'contentChecksum':matches[0].get('contentChecksum')}
    with fits.open(TMP/'mom0.fits',memmap=False) as h:
        hh=h[0] if h[0].data is not None else h[1]; mom=np.asarray(hh.data).copy(); mh=hh.header.copy()
    with fits.open(TMP/'mask.fits',memmap=False) as h:
        hh=h[0] if h[0].data is not None else h[1]; mask=np.asarray(hh.data).copy()
    pix,geo,histruct=project_source_pixels(mom,mask,mh,r)
    pix.to_csv(TMP/'pixels.csv',index=False,float_format='%.12g'); geo.to_csv(TMP/'geometry.csv',index=False,float_format='%.12g')
    subprocess.run([sys.executable,str(ENG/'source_moment0_profile.py'),'--pixels',str(TMP/'pixels.csv'),'--geometry',str(TMP/'geometry.csv'),
                    '--out-profile',str(OUT/'hi_profile.csv'),'--out-report',str(OUT/'hi_profile_report.json')],check=True)
    op,opmeta=normalized_optical(); op.to_csv(OUT/'optical_profile.csv',index=False,float_format='%.12g')
    optsum=json.loads((OPT/'nersc_native_multibrick_probe_summary.json').read_text())
    logm=float(optsum['details_native_multibrick']['log10_Mstar_adopted_median'])
    meta=pd.DataFrame([{'galaxy_id':NAME,'field_id':FIELD,'release_phase':2,'distance_mpc':float(r.dist_h),
                        'logMstar':logm,'logMHI':float(r.log_m_hi_corr),'ell_maj_beams':float(r.ell_maj)/5.0,
                        'robust_sample':True,'full_stellar':True}])
    meta.to_csv(OUT/'meta.csv',index=False,float_format='%.12g')
    normalized={'meta':OUT/'meta.csv','hi_profile':OUT/'hi_profile.csv','optical_profile':OUT/'optical_profile.csv'}
    manifest={
      'schema':'WALLABY_SOURCE_EXPORT_MANIFEST_v1','stage':'SOURCE_ONLY','join_key':'galaxy_id',
      'preflight_scope':'one-object TARGET interface validation only; not the Stage-1 population freeze',
      'validation_release_phase':2,'validation_fields':['NGC5044','NGC4808','Vela'],'response_products_loaded':False,
      'product_ids':{'source_catalogue':'post-stage0/source-only-probe/source_catalogue_safe.csv',
                     'cadc_plane':PLANE,'cadc_observation':OBS,'cadc_mom0':wanted['mom0']['uri'],'cadc_mask':wanted['mask']['uri'],
                     'optical_receipt':'post-stage0/source-optical-nersc-native-multibrick-probe'},
      'column_map':{'source_catalogue':{'name':'galaxy_id','dist_h':'distance_mpc','log_m_hi_corr':'logMHI','ell_maj/5':'ell_maj_beams'},
                    'moment0_mask':'WCS source pixels -> dx_east_arcsec,dy_north_arcsec,moment0_weight -> frozen source_moment0_profile',
                    'optical':'common grz SB_e<0.22; retained r R,SB,SB_e,ellip -> radius_arcsec,mu_r,muerr_r,ellipticity',
                    'stellar_mass':'native-grid DR10 45-MLCR median -> logMstar'},
      'normalized_files':{k:{'path':str(p.relative_to(ROOT)),'sha256':sha(p)} for k,p in normalized.items()},
      'tap_query':tap_query,'source_artifact_hashes':wanted,'hi_structure':histruct,'optical_normalization':opmeta,
    }
    (OUT/'SOURCE_EXPORT_MANIFEST.json').write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n')
    stage=OUT/'stage1'; stage.mkdir(exist_ok=True)
    subprocess.run([sys.executable,str(ENG/'run_source_stage.py'),'--mode','TARGET','--meta',str(OUT/'meta.csv'),
                    '--hi-profile',str(OUT/'hi_profile.csv'),'--optical-profile',str(OUT/'optical_profile.csv'),
                    '--initial-receipt',str(RECEIPT),'--export-manifest',str(OUT/'SOURCE_EXPORT_MANIFEST.json'),'--out-dir',str(stage)],check=True)
    rep={'schema':'WALLABY_STAGE1_TARGET_PREFLIGHT_v1','status':'pass','galaxy_id':NAME,'team_release':TEAM,
         'response_products_used':False,'source_artifacts':wanted,'hi_structure':histruct,'optical_normalization':opmeta,
         'exact_v499_engine_verification':'post-stage0/frozen-v4.99-source-engine/VERIFICATION.json',
         'stage0_receipt':'post-stage0/STAGE0_EXTERNAL_TIMESTAMP_RECEIPT.json',
         'stage1_payload_sha256':sha(stage/'STAGE1_TIMESTAMP_PAYLOAD.json'),
         'note':'This one-object payload is a TARGET interface preflight only and MUST NOT be used as the population Stage-1 timestamp receipt.'}
    (OUT/'PREFLIGHT_REPORT.json').write_text(json.dumps(rep,indent=2,sort_keys=True)+'\n')
    files=sorted(p for p in OUT.rglob('*') if p.is_file() and p.name!='SHA256SUMS.txt')
    with (OUT/'SHA256SUMS.txt').open('w') as f:
        for p in files: f.write(f'{sha(p)}  {p.relative_to(ROOT)}\n')
    print(json.dumps(rep,indent=2,sort_keys=True))

if __name__=='__main__': main()
