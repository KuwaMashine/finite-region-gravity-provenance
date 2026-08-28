#!/usr/bin/env python3
"""Second mechanical repair for the TARGET source-stage preflight.

Keeps the closed source_data productID probe from v2 and repairs only the WCS
frame mismatch: the catalogue centre is transformed into the moment-0 celestial
WCS frame before spherical east/north offsets are evaluated.
"""
from __future__ import annotations
import csv, importlib.util, urllib.parse, urllib.request
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.wcs.utils import proj_plane_pixel_scales

BASE_PATH=Path(__file__).with_name('source_stage_target_preflight.py')
spec=importlib.util.spec_from_file_location('source_stage_target_preflight_base',BASE_PATH)
base=importlib.util.module_from_spec(spec); spec.loader.exec_module(base)

CANDIDATES=(
    'source_data_NGC5044_TR3','source_data_NGC_5044_TR3',
    'source_data_NGC5044_DR3','source_data_NGC_5044_DR3',
    'source_data_N5044_TR3','source_data_N5044_DR3',
)

def source_artifacts_closed_probe():
    attempts=[]
    for plane in CANDIDATES:
        q=("SELECT o.observationID,p.productID,a.uri,a.productType,a.contentType,a.contentLength,a.contentChecksum "
           "FROM caom2.Observation AS o JOIN caom2.Plane AS p ON o.obsID=p.obsID "
           "JOIN caom2.Artifact AS a ON p.planeID=a.planeID "
           f"WHERE o.collection='WALLABY' AND o.observationID='{base.OBS}' AND p.productID='{plane}'")
        data=urllib.parse.urlencode({'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'csv','QUERY':q}).encode()
        req=urllib.request.Request(base.TAP,data=data,headers={'User-Agent':'finite-region-gravity-source-preflight/3'})
        with urllib.request.urlopen(req,timeout=180) as r: text=r.read().decode('utf-8-sig')
        rows=list(csv.DictReader(text.splitlines())); attempts.append({'productID':plane,'n_artifacts':len(rows)})
        if not rows: continue
        bad=[r for r in rows if any(x in (r.get('uri','')+' '+r.get('productID','')).lower() for x in base.FORBIDDEN)]
        if bad: raise RuntimeError('source-plane artifact firewall violation')
        if any(r.get('productID')!=plane for r in rows): raise RuntimeError('CADC returned an unrequested productID')
        mom=[r for r in rows if r.get('uri','').endswith('_mom0.fits')]
        mask=[r for r in rows if r.get('uri','').endswith('_mask.fits')]
        if len(mom)!=1 or len(mask)!=1: raise RuntimeError(f'{plane}: source plane lacks unique mom0/mask artifacts')
        prefix=mom[0]['uri'].split('/')[-1][:-len('_mom0.fits')]
        if mask[0]['uri'].split('/')[-1] != prefix+'_mask.fits': raise RuntimeError('mom0/mask source prefixes differ')
        if not prefix.startswith(base.OBS+'_'): raise RuntimeError(f'unexpected source artifact prefix {prefix}')
        base.PLANE=plane; base.OBJ=prefix
        return rows,text,q
    raise RuntimeError('no artifacts under closed source_data productID candidates: '+repr(attempts))

def project_source_pixels_frame_safe(mom,mask,header,r):
    m=np.asarray(mom).squeeze(); k=np.asarray(mask).squeeze()
    if m.ndim!=2 or k.ndim!=3 or tuple(k.shape[-2:])!=tuple(m.shape):
        raise RuntimeError(f'unexpected moment0/mask shapes {m.shape} {k.shape}')
    mask2=np.any(k>0,axis=0); yy,xx=np.nonzero(mask2 & np.isfinite(m))
    if len(xx)<8: raise RuntimeError('fewer than eight finite moment0 pixels in source footprint')
    w=WCS(header).celestial; sc=w.pixel_to_world(xx.astype(float),yy.astype(float))
    center_icrs=SkyCoord(float(r.ra)*u.deg,float(r.dec)*u.deg,frame='icrs')
    center=center_icrs.transform_to(sc.frame)
    de,dn=center.spherical_offsets_to(sc)
    pix=pd.DataFrame({'galaxy_id':base.NAME,'dx_east_arcsec':de.to_value(u.arcsec),
                      'dy_north_arcsec':dn.to_value(u.arcsec),'moment0_weight':m[yy,xx].astype(float)})
    scales=np.abs(np.asarray(proj_plane_pixel_scales(w),float))*3600.0; ps=float(np.mean(scales))
    geo=pd.DataFrame([{'galaxy_id':base.NAME,'ell_maj_arcsec':6.0*float(r.ell_maj),
                       'ell_min_arcsec':6.0*float(r.ell_min),'ell_pa_deg':float(r.ell_pa),
                       'pixel_scale_arcsec':ps}])
    return pix,geo,{'moment0_shape':list(m.shape),'mask_shape':list(k.shape),
                    'projected_mask_finite_pixels':int(len(xx)),
                    'pixel_scale_arcsec_xy':[float(x) for x in scales],
                    'pixel_scale_arcsec_mean':ps,
                    'moment0_wcs_frame':str(sc.frame.name),
                    'catalogue_center_input_frame':'icrs',
                    'catalogue_center_transformed_to_moment0_frame':True}

base.tap_source_artifacts=source_artifacts_closed_probe
base.project_source_pixels=project_source_pixels_frame_safe
base.main()
