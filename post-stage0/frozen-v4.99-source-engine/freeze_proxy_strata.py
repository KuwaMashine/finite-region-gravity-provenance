#!/usr/bin/env python3
"""Freeze source-only kappa_HI rank tertiles before any response product is opened."""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd

REQ = {"galaxy_id", "field_id", "R_HI90", "R_dstar"}

def sha256(p: Path) -> str:
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--predictors',type=Path,required=True)
    ap.add_argument('--out-csv',type=Path,required=True)
    ap.add_argument('--out-json',type=Path,required=True)
    a=ap.parse_args()
    df=pd.read_csv(a.predictors)
    miss=REQ-df.columns.to_series().index.to_list() if False else REQ-set(df.columns)
    if miss: raise SystemExit(f"missing columns: {sorted(miss)}")
    if df['galaxy_id'].duplicated().any(): raise SystemExit('one predictor row per galaxy is required')
    for c in ['R_HI90','R_dstar']:
        df[c]=pd.to_numeric(df[c],errors='coerce')
        if (~np.isfinite(df[c])).any() or (df[c]<=0).any(): raise SystemExit(f'{c} must be finite and positive')
    df['kappa_HI']=df['R_HI90']/df['R_dstar']
    df['galaxy_id']=df['galaxy_id'].astype(str)
    df['field_id']=df['field_id'].astype(str)
    df=df.sort_values(['kappa_HI','galaxy_id'],kind='mergesort').reset_index(drop=True)
    chunks=np.array_split(np.arange(len(df)),3)
    labels=np.empty(len(df),dtype=object)
    for inds,lab in zip(chunks,['low','middle','high']): labels[inds]=lab
    df['stratum']=labels
    out=df[['galaxy_id','field_id','R_HI90','R_dstar','kappa_HI','stratum']].copy()
    out.to_csv(a.out_csv,index=False,float_format='%.12g')
    # realized value ranges are audit metadata, not re-optimized cuts
    ranges={}
    for lab in ['low','middle','high']:
        z=out[out.stratum==lab]['kappa_HI']
        ranges[lab]={'n':int(len(z)),'min':float(z.min()) if len(z) else None,'max':float(z.max()) if len(z) else None}
    meta={
      'rule':'global rank tertiles; stable sort by (kappa_HI, galaxy_id); numpy.array_split sizes',
      'predictor_input':str(a.predictors),
      'predictor_input_sha256':sha256(a.predictors),
      'output_csv':str(a.out_csv),
      'output_csv_sha256':sha256(a.out_csv),
      'n_galaxies':int(len(out)),
      'ranges_audit_only':ranges,
      'response_products_used':False,
    }
    a.out_json.write_text(json.dumps(meta,indent=2,sort_keys=True)+'\n')
if __name__=='__main__': main()
