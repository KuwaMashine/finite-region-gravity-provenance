#!/usr/bin/env python3
"""Freeze the exact normalized SOURCE_ONLY export manifest before TARGET source-stage execution.

This script consumes only already-committed source-side products. It does not read
or name any validation kinematic file. The manifest is the input-boundary object
required by the byte-verified v4.99 holdout engine.
"""
from __future__ import annotations
import hashlib, json
from pathlib import Path

ROOT=Path.cwd()
META=ROOT/'post-stage0/source-optical-stage1/meta.csv'
HI=ROOT/'post-stage0/source-hi-stage1/hi_profile.csv'
OPT=ROOT/'post-stage0/source-optical-stage1/optical_profile.csv'
HI_STATUS=ROOT/'post-stage0/source-hi-stage1/source_hi_status.csv'
OPT_AUDIT=ROOT/'post-stage0/source-optical-stage1/frame_audit.csv'
SOURCE_CAT=ROOT/'post-stage0/source-only-probe/source_catalogue_safe.csv'
FRAME_LOCK=ROOT/'post-stage0/SOURCE_OPTICAL_PRODUCTION_FRAME_LOCK.json'
OUT=ROOT/'post-stage0/SOURCE_EXPORT_MANIFEST.json'

def sha(p:Path)->str:
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def main():
    required=[META,HI,OPT,HI_STATUS,OPT_AUDIT,SOURCE_CAT,FRAME_LOCK]
    miss=[str(p) for p in required if not p.exists()]
    if miss: raise SystemExit(f'missing committed source-stage inputs: {miss}')
    obj={
      'schema':'WALLABY_SOURCE_EXPORT_MANIFEST_v1',
      'stage':'SOURCE_ONLY',
      'join_key':'galaxy_id',
      'validation_release_phase':2,
      'validation_fields':['NGC5044','NGC4808','Vela'],
      'response_products_loaded':False,
      'product_ids':{
        'wallaby_source_catalogue':{
          'service':'CANFAR/CIRADA TAP cirada.Wallaby_dr2_source_catalogue',
          'normalized_source_catalogue_receipt':str(SOURCE_CAT.relative_to(ROOT)),
          'normalized_source_catalogue_sha256':sha(SOURCE_CAT),
          'final_releases':['NGC 4808 TR1','NGC 5044 TR3','Vela TR1']
        },
        'wallaby_source_data_hi':{
          'service':'CADC public WALLABY source_data moment-0 and source-mask products',
          'per_object_url_and_sha256_receipt':str(HI_STATUS.relative_to(ROOT)),
          'per_object_url_and_sha256_receipt_sha256':sha(HI_STATUS)
        },
        'legacy_survey_dr10_optical':{
          'service':'NERSC DESI Legacy Imaging Surveys DR10 south coadd g/r/z images',
          'per_frame_brick_url_sha256_and_geometry_receipt':str(OPT_AUDIT.relative_to(ROOT)),
          'per_frame_brick_url_sha256_and_geometry_receipt_sha256':sha(OPT_AUDIT),
          'production_frame_lock':str(FRAME_LOCK.relative_to(ROOT)),
          'production_frame_lock_sha256':sha(FRAME_LOCK)
        },
        'galactic_extinction':{
          'service':'IRSA SFD dust query at frozen WALLABY source-catalogue sky centroid',
          'normalized_extinction_values_are_carried_in':'post-stage0/source-optical-stage1/optical_status.csv',
          'rule':'source coordinate only; Schlafly-Finkbeiner/Legacy g,r,z coefficients as frozen in Stage-0A optical implementation'
        }
      },
      'column_map':{
        'meta':{
          'galaxy_id':'source catalogue name',
          'field_id':'team_release -> frozen validation field id',
          'release_phase':'constant 2 for final validation releases',
          'distance_mpc':'dist_h',
          'logMstar':'median of frozen 45 MLCR log-mass estimates from DR10 common-aperture g/r/z photometry',
          'logMHI':'log_m_hi_corr',
          'ell_maj_beams':'ell_maj / 5.0 for 6 arcsec source-catalogue pixels and 30 arcsec beam',
          'robust_sample':'qflag==0 + finite source metadata + successful exact frozen source-moment0 profile construction',
          'full_stellar':'all-three-band DR10 + frozen AutoProf/common-aperture/MLCR + registered finite-frame convergence'
        },
        'hi_profile':{
          'galaxy_id':'source catalogue name',
          'rad_hi_source_arcsec':'exact frozen source_moment0_profile.py annulus centre',
          'sigma_hi_source_weight':'signed source moment-0 annular density; no pixel clipping'
        },
        'optical_profile':{
          'galaxy_id':'source catalogue name',
          'radius_arcsec':'AutoProf forced-profile R within registered common uncertainty aperture',
          'mu_r':'foreground-corrected r-band surface brightness profile input convention used by frozen Stage-1 estimator',
          'muerr_r':'r-band surface-brightness uncertainty',
          'ellipticity':'r-band AutoProf geometric ellipticity'
        }
      },
      'normalized_files':{
        'meta':{'path':str(META.relative_to(ROOT)),'sha256':sha(META)},
        'hi_profile':{'path':str(HI.relative_to(ROOT)),'sha256':sha(HI)},
        'optical_profile':{'path':str(OPT.relative_to(ROOT)),'sha256':sha(OPT)}
      },
      'source_side_evidence':{
        'hi_status_sha256':sha(HI_STATUS),
        'optical_frame_audit_sha256':sha(OPT_AUDIT),
        'frame_lock_sha256':sha(FRAME_LOCK)
      },
      'firewall':'No response/rotation product has been loaded. This manifest freezes exact normalized source inputs before TARGET source-stage execution.'
    }
    OUT.write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':'SOURCE_EXPORT_MANIFEST_READY','path':str(OUT),'sha256':sha(OUT),'response_products_loaded':False},sort_keys=True))
if __name__=='__main__': main()
