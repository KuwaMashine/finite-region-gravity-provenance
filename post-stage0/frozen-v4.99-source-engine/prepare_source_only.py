#!/usr/bin/env python3
"""Prepare gravity-blind WALLABY predictor products before opening response rows.

Inputs are normalized exports of public/source-side products only.  No Vrot,
gobs, RAR residual or flat-curve quantity is accepted by this script.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
from wallaby_source_geometry import optical_exponential_equivalent_scale, kappa_hi_from_source_profiles
from baryonic_gravity import normalize_hi_profile_to_mass
from holdout_integrity import verify_initial_receipt, verify_export_manifest

DEV_FIELDS={"Hydra","Norma","NGC4636"}
VAL_FIELDS={"NGC5044","NGC4808","Vela"}
REQ_META={"galaxy_id","field_id","release_phase","distance_mpc","logMstar","logMHI","ell_maj_beams","robust_sample","full_stellar"}
REQ_HI={"galaxy_id","rad_hi_source_arcsec","sigma_hi_source_weight"}
REQ_OPT={"galaxy_id","radius_arcsec","mu_r","muerr_r","ellipticity"}
FORBIDDEN={"Vrot_model","vrot_kms","gobs","gbar","delta_ref","rar_residual","Vmax","Vflat"}


def _bool(v):
    if isinstance(v,(bool,np.bool_)): return bool(v)
    if isinstance(v,(int,np.integer,float,np.floating)): return bool(v)
    return str(v).strip().lower() in {"1","true","t","yes","y"}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--meta',type=Path,required=True)
    ap.add_argument('--hi-profile',type=Path,required=True)
    ap.add_argument('--optical-profile',type=Path,required=True)
    ap.add_argument('--role',choices=['development','validation'],required=True)
    ap.add_argument('--mode',choices=['MOCK','TARGET'],default='TARGET')
    ap.add_argument('--initial-receipt',type=Path)
    ap.add_argument('--export-manifest',type=Path)
    ap.add_argument('--out-predictors',type=Path,required=True)
    ap.add_argument('--out-baryonic-hi',type=Path,required=True)
    ap.add_argument('--out-report',type=Path,required=True)
    a=ap.parse_args()
    if a.mode=='TARGET':
        if a.initial_receipt is None: raise SystemExit('TARGET LOCK: source-only load requires initial external timestamp receipt')
        if a.export_manifest is None: raise SystemExit('TARGET LOCK: source-only load requires source export manifest')
        verify_initial_receipt(a.initial_receipt, Path(__file__).resolve().parent/'PRETARGET_MANIFEST_SHA256.txt')
        verify_export_manifest(a.export_manifest, {'meta':a.meta,'hi_profile':a.hi_profile,'optical_profile':a.optical_profile}, 'SOURCE_ONLY')
    meta=pd.read_csv(a.meta,dtype={'galaxy_id':str}); hi=pd.read_csv(a.hi_profile,dtype={'galaxy_id':str}); op=pd.read_csv(a.optical_profile,dtype={'galaxy_id':str})
    for df,n,req in [(meta,'meta',REQ_META),(hi,'hi',REQ_HI),(op,'optical',REQ_OPT)]:
        miss=req-set(df.columns)
        if miss: raise SystemExit(f'{n} missing columns: {sorted(miss)}')
        leak=FORBIDDEN&set(df.columns)
        if leak: raise SystemExit(f'SOURCE LOCK: {n} contains forbidden response columns: {sorted(leak)}')
    allowed=DEV_FIELDS if a.role=='development' else VAL_FIELDS
    phase=1 if a.role=='development' else 2
    m=meta[(meta.release_phase.astype(int)==phase)&meta.field_id.astype(str).isin(allowed)].copy()
    m=m[m.robust_sample.map(_bool)&m.full_stellar.map(_bool)].copy()
    if m.galaxy_id.duplicated().any(): raise SystemExit('meta contains duplicate galaxy_id')
    preds=[]; rings=[]; failures=[]
    for row in m.sort_values('galaxy_id').itertuples(index=False):
        gid=str(row.galaxy_id)
        zh=hi[hi.galaxy_id==gid].sort_values('rad_hi_source_arcsec')
        zo=op[op.galaxy_id==gid].sort_values('radius_arcsec')
        try:
            if len(zh)<8: raise ValueError('fewer than 8 H I profile rings')
            if len(zo)<5: raise ValueError('fewer than 5 optical profile points')
            os=optical_exponential_equivalent_scale(zo.radius_arcsec,zo.mu_r,zo.muerr_r,zo.ellipticity)
            kh=kappa_hi_from_source_profiles(zh.rad_hi_source_arcsec,zh.sigma_hi_source_weight,gid,os.rd_star_arcsec)
            mhi=10.0**float(row.logMHI)
            bnorm,scale=normalize_hi_profile_to_mass(kh['baryonic_profile'],float(row.distance_mpc),mhi)
            cm=kh['crossfit_meta']
            preds.append({
                'galaxy_id':gid,'field_id':str(row.field_id),'release_phase':phase,
                'distance_mpc':float(row.distance_mpc),'logMstar':float(row.logMstar),'logMHI':float(row.logMHI),
                'ell_maj_beams':float(row.ell_maj_beams),'R_HI90':float(kh['R_HI90_arcsec']),'R_dstar':float(os.rd_star_arcsec),'R_HI90_arcsec':float(kh['R_HI90_arcsec']),
                'R50_r_arcsec':float(os.r50_arcsec),'R_dstar_arcsec':float(os.rd_star_arcsec),
                'kappa_HI':float(kh['kappa_HI']),'crossfit_phase_bit':int(cm['phase_bit']),
                'crossfit_dropped_source_index':cm['dropped_source_index'],
                'crossfit_n_source_used':int(cm['n_source_used']),'crossfit_n_paired_bins':int(cm['n_paired_bins']),
                'hi_mass_renorm_scale':float(scale),
            })
            for j,(rin,rout,sig) in enumerate(zip(bnorm.rin,bnorm.rout,bnorm.sigma)):
                rings.append({'galaxy_id':gid,'pair_index':j,'rin_arcsec':float(rin),'rout_arcsec':float(rout),'sigma_hi_msun_pc2':float(sig)})
        except Exception as e:
            failures.append({'galaxy_id':gid,'reason':str(e)})
    if not preds: raise SystemExit('no source-only predictor passed')
    pd.DataFrame(preds).to_csv(a.out_predictors,index=False,float_format='%.12g')
    pd.DataFrame(rings).to_csv(a.out_baryonic_hi,index=False,float_format='%.12g')
    rep={'role':a.role,'release_phase':phase,'allowed_fields':sorted(allowed),'n_meta_eligible':int(len(m)),'n_predictor_pass':len(preds),'n_fail':len(failures),'failures':failures,
         'firewall':'No response/rotation columns accepted; H I profile is explicitly source-derived; strata are frozen later from kappa_HI before kinematic response load.', 'mode':a.mode}
    a.out_report.write_text(json.dumps(rep,indent=2,sort_keys=True)+'\n')
    print(json.dumps(rep,sort_keys=True))
if __name__=='__main__': main()
