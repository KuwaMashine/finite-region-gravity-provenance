#!/usr/bin/env python3
"""Run the WALLABY validation source-only stage without opening kinematics.

TARGET mode requires the Stage-0 external timestamp receipt. The output payload
contains hashes to be externally timestamped at Stage 1 before Vrot is opened.
"""
from __future__ import annotations
import argparse, json, subprocess
from pathlib import Path
from holdout_integrity import sha256, verify_initial_receipt

HERE=Path(__file__).resolve().parent

def run(cmd):
    subprocess.run([str(x) for x in cmd],check=True)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--mode',choices=['MOCK','TARGET'],default='TARGET')
    ap.add_argument('--meta',type=Path,required=True)
    ap.add_argument('--hi-profile',type=Path,required=True)
    ap.add_argument('--optical-profile',type=Path,required=True)
    ap.add_argument('--initial-receipt',type=Path)
    ap.add_argument('--export-manifest',type=Path)
    ap.add_argument('--out-dir',type=Path,required=True)
    a=ap.parse_args(); a.out_dir.mkdir(parents=True,exist_ok=True)
    if a.mode=='TARGET':
        if a.initial_receipt is None: raise SystemExit('TARGET LOCK: Stage-0 receipt required before source-only load')
        if a.export_manifest is None: raise SystemExit('TARGET LOCK: source export manifest required before source-only load')
        verify_initial_receipt(a.initial_receipt,HERE/'PRETARGET_MANIFEST_SHA256.txt')
    pred=a.out_dir/'source_predictors.csv'; bhi=a.out_dir/'baryonic_hi_fold.csv'; prep=a.out_dir/'source_prepare_report.json'
    strata=a.out_dir/'proxy_strata.csv'; smeta=a.out_dir/'proxy_strata_metadata.json'; beam=a.out_dir/'source_beam_report.json'
    cmd=[HERE/'prepare_source_only.py','--meta',a.meta,'--hi-profile',a.hi_profile,'--optical-profile',a.optical_profile,
         '--role','validation','--mode',a.mode,'--out-predictors',pred,'--out-baryonic-hi',bhi,'--out-report',prep]
    if a.initial_receipt is not None: cmd += ['--initial-receipt',a.initial_receipt]
    if a.export_manifest is not None: cmd += ['--export-manifest',a.export_manifest]
    run(cmd)
    run([HERE/'freeze_proxy_strata.py','--predictors',pred,'--out-csv',strata,'--out-json',smeta])
    run([HERE/'source_beam_design_control.py','--predictors',pred,'--out-json',beam])
    manifest=HERE/'PRETARGET_MANIFEST_SHA256.txt'
    payload={
      'purpose':'External Stage-1 timestamp payload; create receipt before opening validation kinematics',
      'frozen_manifest_file_sha256':sha256(manifest),
      'proxy_strata_map_sha256':sha256(strata),
      'proxy_strata_metadata_sha256':sha256(smeta),
      'source_beam_report_sha256':sha256(beam),
      'source_predictors_sha256':sha256(pred),
      'baryonic_hi_fold_sha256':sha256(bhi),
      'source_prepare_report_sha256':sha256(prep),
      'response_products_used':False,
    }
    out=a.out_dir/'STAGE1_TIMESTAMP_PAYLOAD.json'; out.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':'SOURCE_STAGE_COMPLETE','stage1_payload':str(out),**payload},sort_keys=True))
if __name__=='__main__': main()
