#!/usr/bin/env python3
"""Exact-v4.99 TARGET source-stage preflight on WALLABY J130150+041953.

This source was already used for the nuisance-aware 1024->1280 optical
convergence check. The script adds only source-plane CADC H I products, builds
the frozen signed moment-0 profile, normalizes the existing 1024 DR10 profile,
and invokes the byte-verified v4.99 TARGET source-stage engine. No response or
kinematic product is queried or read.
"""
from __future__ import annotations
import csv, hashlib, json, math, subprocess, sys, urllib.parse, urllib.request
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.wcs.utils import proj_plane_pixel_scales

ROOT=Path.cwd(); OUT=ROOT/'post-stage0/source-stage-target-preflight-ngc4808'; TMP=Path('/tmp/wallaby_target_preflight_ngc4808')
ENG=ROOT/'post-stage0/frozen-v4.99-source-engine'; RECEIPT=ROOT/'post-stage0/STAGE0_EXTERNAL_TIMESTAMP_RECEIPT.json'
NAME='WALLABY J130150+041953'; OBS='WALLABY_J130150+041953'; TEAM='NGC 4808 TR1'; FIELD='NGC4808'
OPT=ROOT/'post-stage0/source-optical-nersc-frame-convergence/n1024'
OPTSUM=ROOT/'post-stage0/source-optical-nersc-frame-convergence/nersc_frame_convergence_summary.json'
CAT=ROOT/'post-stage0/source-only-probe/source_catalogue_safe.csv'
BASE='https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/data/pub/WALLABY'; TAP='https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/argus/sync'
CANDIDATES=('source_data_NGC4808_TR1','source_data_NGC_4808_TR1','source_data_NGC4808_DR1','source_data_NGC_4808_DR1','source_data_N4808_TR1','source_data_N4808_DR1')
FORBIDDEN=('kinematic','kinmodel','rotcur','model_data','vrot','gobs','gbar','rar_residual','vmax','vflat')

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def download(url,p):
    req=urllib.request.Request(url,headers={'User-Agent':'finite-region-gravity-source-preflight-ngc4808/1'})
    with urllib.request.urlopen(req,timeout=180) as r, Path(p).open('wb') as f:
        while True:
            b=r.read(1<<20)
            if not b: break
            f.write(b)
    if Path(p).stat().st_size==0: raise RuntimeError(f'zero-byte download {url}')

def source_artifacts():
    attempts=[]
    for plane in CANDIDATES:
        q=("SELECT o.observationID,p.productID,a.uri,a.productType,a.contentType,a.contentLength,a.contentChecksum "
           "FROM caom2.Observation AS o JOIN caom2.Plane AS p ON o.obsID=p.obsID JOIN caom2.Artifact AS a ON p.planeID=a.planeID "
           f"WHERE o.collection='WALLABY' AND o.observationID='{OBS}' AND p.productID='{plane}'")
        data=urllib.parse.urlencode({'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'csv','QUERY':q}).encode()
        req=urllib.request.Request(TAP,data=data,headers={'User-Agent':'finite-region-gravity-source-preflight-ngc4808/1'})
        with urllib.request.urlopen(req,timeout=180) as r: text=r.read().decode('utf-8-sig')
        rows=list(csv.DictReader(text.splitlines())); attempts.append({'productID':plane,'n_artifacts':len(rows)})
        if not rows: continue
        if any(any(x in (z.get('uri','')+' '+z.get('productID','')).lower() for x in FORBIDDEN) for z in rows): raise RuntimeError('source-plane firewall violation')
        if any(z.get('productID')!=plane for z in rows): raise RuntimeError('CADC returned an unrequested plane')
        mom=[z for z in rows if z.get('uri','').endswith('_mom0.fits')]; mask=[z for z in rows if z.get('uri','').endswith('_mask.fits')]
        if len(mom)!=1 or len(mask)!=1: raise RuntimeError(f'{plane}: source plane lacks unique mom0/mask')
        prefix=mom[0]['uri'].split('/')[-1][:-len('_mom0.fits')]
        if mask[0]['uri'].split('/')[-1]!=prefix+'_mask.fits' or not prefix.startswith(OBS+'_'): raise RuntimeError('source artifact prefix mismatch')
        return plane,prefix,rows,text,q,attempts
    raise RuntimeError('no artifacts under closed source_data candidates: '+repr(attempts))

def source_row():
    df=pd.read_csv(CAT); z=df[(df.name==NAME)&(df.team_release==TEAM)]
    if len(z)!=1: raise RuntimeError(f'expected one source row, got {len(z)}')
    r=z.iloc[0]
    if float(r.qflag)!=0.0: raise RuntimeError('source qflag != 0')
    for c in ('ra','dec','dist_h','log_m_hi_corr','ell_maj','ell_min','ell_pa'):
        if not math.isfinite(float(r[c])): raise RuntimeError(f'nonfinite {c}')
    return r

def source_pixels(mom,mask,header,r):
    m=np.asarray(mom).squeeze(); k=np.asarray(mask).squeeze()
    if m.ndim!=2 or k.ndim!=3 or tuple(k.shape[-2:])!=tuple(m.shape): raise RuntimeError(f'unexpected shapes {m.shape} {k.shape}')
    mask2=np.any(k>0,axis=0); yy,xx=np.nonzero(mask2 & np.isfinite(m))
    if len(xx)<8: raise RuntimeError('fewer than eight finite source pixels')
    w=WCS(header).celestial; sc=w.pixel_to_world(xx.astype(float),yy.astype(float))
    center=SkyCoord(float(r.ra)*u.deg,float(r.dec)*u.deg,frame='icrs').transform_to(sc.frame)
    de,dn=center.spherical_offsets_to(sc)
    pix=pd.DataFrame({'galaxy_id':NAME,'dx_east_arcsec':de.to_value(u.arcsec),'dy_north_arcsec':dn.to_value(u.arcsec),'moment0_weight':m[yy,xx].astype(float)})
    scales=np.abs(np.asarray(proj_plane_pixel_scales(w),float))*3600.; ps=float(np.mean(scales))
    geo=pd.DataFrame([{'galaxy_id':NAME,'ell_maj_arcsec':6.*float(r.ell_maj),'ell_min_arcsec':6.*float(r.ell_min),'ell_pa_deg':float(r.ell_pa),'pixel_scale_arcsec':ps}])
    return pix,geo,{'moment0_shape':list(m.shape),'mask_shape':list(k.shape),'projected_mask_finite_pixels':int(len(xx)),
                    'pixel_scale_arcsec_xy':[float(x) for x in scales],'pixel_scale_arcsec_mean':ps,'moment0_wcs_frame':str(sc.frame.name),
                    'catalogue_center_input_frame':'icrs','catalogue_center_transformed_to_moment0_frame':True}

def read_prof(p):
    with Path(p).open() as f: f.readline(); return list(csv.DictReader(f))

def optical_profile():
    P={b:read_prof(OPT/f'{b}.prof') for b in 'grz'}; n=min(len(P[b]) for b in 'grz'); keep=n
    for i in range(n):
        rs=[float(P[b][i]['R']) for b in 'grz']; es=[float(P[b][i]['SB_e']) for b in 'grz']
        if max(rs)-min(rs)>1e-8: raise RuntimeError(f'optical radius mismatch row {i}')
        if any((not math.isfinite(x)) or x>=0.22 for x in es): keep=i; break
    if keep<5: raise RuntimeError('fewer than five common optical rows')
    out=pd.DataFrame([{'galaxy_id':NAME,'radius_arcsec':float(z['R']),'mu_r':float(z['SB']),'muerr_r':float(z['SB_e']),'ellipticity':float(z['ellip'])} for z in P['r'][:keep]])
    return out,{'source':'post-stage0/source-optical-nersc-frame-convergence/n1024','n_common_retained':keep,'common_aperture_radius_arcsec':float(P['r'][keep-1]['R'])}

def logmstar_1024():
    s=json.loads(OPTSUM.read_text()); rows=[z for z in s['results'] if int(z['cutout_size_pix'])==1024]
    if len(rows)!=1: raise RuntimeError('cannot resolve unique 1024 optical summary row')
    return float(rows[0]['log10_Mstar_adopted_median'])

def main():
    OUT.mkdir(parents=True,exist_ok=True); TMP.mkdir(parents=True,exist_ok=True); r=source_row()
    plane,prefix,arts,art_csv,tap_query,attempts=source_artifacts(); (OUT/'cadc_source_artifacts.csv').write_text(art_csv)
    wanted={}
    for kind,suffix in [('mom0','_mom0.fits'),('mask','_mask.fits')]:
        z=[a for a in arts if a.get('uri','').endswith(suffix)]; a=z[0]; fn=a['uri'].split('/')[-1]; p=TMP/f'{kind}.fits'; download(f'{BASE}/{fn}',p)
        wanted[kind]={'uri':a['uri'],'filename':fn,'sha256':sha(p),'contentChecksum':a.get('contentChecksum')}
    with fits.open(TMP/'mom0.fits',memmap=False) as h:
        hh=h[0] if h[0].data is not None else h[1]; mom=np.asarray(hh.data).copy(); mh=hh.header.copy()
    with fits.open(TMP/'mask.fits',memmap=False) as h:
        hh=h[0] if h[0].data is not None else h[1]; mask=np.asarray(hh.data).copy()
    pix,geo,histruct=source_pixels(mom,mask,mh,r); pix.to_csv(TMP/'pixels.csv',index=False,float_format='%.12g'); geo.to_csv(TMP/'geometry.csv',index=False,float_format='%.12g')
    subprocess.run([sys.executable,str(ENG/'source_moment0_profile.py'),'--pixels',str(TMP/'pixels.csv'),'--geometry',str(TMP/'geometry.csv'),
                    '--out-profile',str(OUT/'hi_profile.csv'),'--out-report',str(OUT/'hi_profile_report.json')],check=True)
    op,opmeta=optical_profile(); op.to_csv(OUT/'optical_profile.csv',index=False,float_format='%.12g')
    meta=pd.DataFrame([{'galaxy_id':NAME,'field_id':FIELD,'release_phase':2,'distance_mpc':float(r.dist_h),'logMstar':logmstar_1024(),
                        'logMHI':float(r.log_m_hi_corr),'ell_maj_beams':float(r.ell_maj)/5.,'robust_sample':True,'full_stellar':True}])
    meta.to_csv(OUT/'meta.csv',index=False,float_format='%.12g')
    norm={'meta':OUT/'meta.csv','hi_profile':OUT/'hi_profile.csv','optical_profile':OUT/'optical_profile.csv'}
    manifest={'schema':'WALLABY_SOURCE_EXPORT_MANIFEST_v1','stage':'SOURCE_ONLY','join_key':'galaxy_id','preflight_scope':'one-object TARGET interface validation only; not Stage-1 population freeze',
              'validation_release_phase':2,'validation_fields':['NGC5044','NGC4808','Vela'],'response_products_loaded':False,
              'product_ids':{'source_catalogue':'post-stage0/source-only-probe/source_catalogue_safe.csv','cadc_plane':plane,'cadc_observation':OBS,
                             'cadc_mom0':wanted['mom0']['uri'],'cadc_mask':wanted['mask']['uri'],'optical_receipt':'post-stage0/source-optical-nersc-frame-convergence/n1024'},
              'column_map':{'source_catalogue':{'name':'galaxy_id','dist_h':'distance_mpc','log_m_hi_corr':'logMHI','ell_maj/5':'ell_maj_beams'},
                            'moment0_mask':'WCS source pixels -> signed frozen source_moment0_profile','optical':'common grz SB_e<0.22 retained r R,SB,SB_e,ellip',
                            'stellar_mass':'1024 direct-NERSC 45-MLCR median -> logMstar'},
              'normalized_files':{k:{'path':str(p.relative_to(ROOT)),'sha256':sha(p)} for k,p in norm.items()},'tap_query':tap_query,
              'source_product_id_attempts':attempts,'source_artifact_hashes':wanted,'hi_structure':histruct,'optical_normalization':opmeta}
    (OUT/'SOURCE_EXPORT_MANIFEST.json').write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n')
    stage=OUT/'stage1'; stage.mkdir(exist_ok=True)
    subprocess.run([sys.executable,str(ENG/'run_source_stage.py'),'--mode','TARGET','--meta',str(OUT/'meta.csv'),'--hi-profile',str(OUT/'hi_profile.csv'),
                    '--optical-profile',str(OUT/'optical_profile.csv'),'--initial-receipt',str(RECEIPT),'--export-manifest',str(OUT/'SOURCE_EXPORT_MANIFEST.json'),
                    '--out-dir',str(stage)],check=True)
    pred=pd.read_csv(stage/'source_predictors.csv').iloc[0].to_dict(); prep=json.loads((stage/'source_prepare_report.json').read_text())
    rep={'schema':'WALLABY_STAGE1_TARGET_PREFLIGHT_v1','status':'pass','galaxy_id':NAME,'source_release':TEAM,'response_products_used':False,
         'cadc_source_plane':plane,'source_artifacts':wanted,'hi_structure':histruct,'optical_normalization':opmeta,
         'exact_v499_engine_verification':'post-stage0/frozen-v4.99-source-engine/VERIFICATION.json','stage0_receipt':'post-stage0/STAGE0_EXTERNAL_TIMESTAMP_RECEIPT.json',
         'predictor_row':pred,'source_prepare_report':prep,'stage1_payload_sha256':sha(stage/'STAGE1_TIMESTAMP_PAYLOAD.json'),
         'note':'One-object interface preflight only; its one-object tertile map is not a population Stage-1 receipt.'}
    (OUT/'PREFLIGHT_REPORT.json').write_text(json.dumps(rep,indent=2,sort_keys=True)+'\n')
    files=sorted(p for p in OUT.rglob('*') if p.is_file() and p.name!='SHA256SUMS.txt')
    with (OUT/'SHA256SUMS.txt').open('w') as f:
        for p in files: f.write(f'{sha(p)}  {p.relative_to(ROOT)}\n')
    print(json.dumps(rep,indent=2,sort_keys=True))

if __name__=='__main__': main()
