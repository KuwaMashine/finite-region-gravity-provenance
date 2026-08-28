#!/usr/bin/env python3
"""Source-only WALLABY beam/resolution design control.

This is not a target-outcome analysis and does not alter kappa_HI or its tertiles.
It carries forward the beam/inclination planning surrogate already frozen in the
12-August contract and makes the >=2-beam primary lane auditable before Vrot is
opened. Full cube-level injection/recovery remains a separate R2/Gate-C control.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
from holdout_integrity import sha256

# Frozen planning control for r_true=0.30 from the registered execution contract.
BEAMS=np.array([1.3,2.0,3.0,4.0,6.0],float)
R_REC=np.array([0.206,0.230,0.235,0.241,0.256],float)
R_TRUE=0.30
PRIMARY_BEAMS=2.0


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--predictors',type=Path,required=True)
    ap.add_argument('--out-json',type=Path,required=True)
    a=ap.parse_args()
    df=pd.read_csv(a.predictors,dtype={'galaxy_id':str})
    req={'galaxy_id','field_id','ell_maj_beams'}
    miss=req-set(df.columns)
    if miss: raise SystemExit(f'missing columns: {sorted(miss)}')
    b=pd.to_numeric(df.ell_maj_beams,errors='coerce').to_numpy(float)
    if (~np.isfinite(b)).any() or (b<=0).any(): raise SystemExit('ell_maj_beams must be finite positive')
    rec=np.interp(np.clip(b,BEAMS[0],BEAMS[-1]),BEAMS,R_REC)
    att=rec/R_TRUE
    primary=b>=PRIMARY_BEAMS
    report={
      'source_only':True,
      'predictor_input':str(a.predictors),
      'predictor_input_sha256':sha256(a.predictors),
      'n_galaxies':int(len(df)),
      'primary_beam_rule':'>=2.0 beams across major axis',
      'n_primary':int(primary.sum()),
      'n_below_primary':int((~primary).sum()),
      'beam_counts':{
        'min':float(np.min(b)),'median':float(np.median(b)),'max':float(np.max(b))
      },
      'registered_surrogate':{
        'r_true':R_TRUE,
        'beam_grid':BEAMS.tolist(),
        'median_recovered_r':R_REC.tolist(),
        'interpolation':'piecewise linear; clamp outside registered 1.3--6 beam grid',
        'median_attenuation_factor_in_source_sample':float(np.median(att)),
        'primary_sample_median_attenuation_factor':float(np.median(att[primary])) if np.any(primary) else None,
      },
      'decision_role':'design/robustness only; does not alter kappa_HI, tertile labels, Gate-A/B sign, or p-values',
      'full_cube_injection_status':'still required separately for R2/Gate-C; this source-only report is the frozen Gate-A/B resolution control',
      'response_products_used':False,
    }
    a.out_json.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(json.dumps(report,sort_keys=True))
if __name__=='__main__': main()
