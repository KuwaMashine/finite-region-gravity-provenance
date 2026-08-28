#!/usr/bin/env python3
"""Consolidate production optical shards and patch final Stage-1 source metadata."""
from __future__ import annotations
import argparse,csv,hashlib,json
from collections import Counter,defaultdict
from pathlib import Path

ROOT=Path.cwd(); PLAN=ROOT/'post-stage0/source-optical-stage1-preplan/optical_stage1_plan.csv'; META_HI=ROOT/'post-stage0/source-hi-stage1/source_meta_hi.csv'; OUT=ROOT/'post-stage0/source-optical-stage1'; OUT.mkdir(parents=True,exist_ok=True)
LOCK=ROOT/'post-stage0/SOURCE_OPTICAL_PRODUCTION_FRAME_LOCK.json'; VELA=ROOT/'post-stage0/source-vela-dr10-field-coverage/VELA_DR10_FIELD_COVERAGE.json'

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): h.update(b)
 return h.hexdigest()

def boolish(x): return str(x).strip().lower() in {'1','true','t','yes','y'}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--shards-dir',type=Path,required=True);a=ap.parse_args()
 plan=list(csv.DictReader(PLAN.open())); candidate={r['galaxy_id']:r for r in plan if r['dr10_optical_field_status']=='candidate'}; vela={r['galaxy_id']:r for r in plan if r['dr10_optical_field_status']=='no_dr10_field_coverage'}
 if len(candidate)+len(vela)!=len(plan): raise SystemExit('unexpected optical preplan status')
 status_files=sorted(a.shards_dir.rglob('status.csv')); profile_files=sorted(a.shards_dir.rglob('optical_profile.csv')); audit_files=sorted(a.shards_dir.rglob('frame_audit.csv'))
 if not status_files: raise SystemExit('no optical shard status files')
 statuses=[]
 for p in status_files: statuses+=list(csv.DictReader(p.open()))
 ids=[r['galaxy_id'] for r in statuses]; dup=[g for g,n in Counter(ids).items() if n!=1]
 if dup: raise SystemExit(f'duplicate optical status ids {dup[:10]}')
 got=set(ids); exp=set(candidate)
 if got!=exp: raise SystemExit(f'optical shard coverage mismatch missing={len(exp-got)} extra={len(got-exp)}')
 unresolved=[r for r in statuses if r['status']=='unresolved_network']; impl=[r for r in statuses if r['status']=='implementation_failure']
 if unresolved or impl: raise SystemExit(f'cannot freeze optical stage unresolved_network={len(unresolved)} implementation_failure={len(impl)}')
 profiles=[]
 for p in profile_files: profiles+=list(csv.DictReader(p.open()))
 pc=Counter(r['galaxy_id'] for r in profiles); passed={r['galaxy_id'] for r in statuses if r['status']=='pass' and boolish(r['full_stellar'])}
 bad=[g for g in passed if pc[g]<5]; extra=set(pc)-passed
 if bad or extra: raise SystemExit(f'optical profile consistency failure short={bad[:10]} extra={list(extra)[:10]}')
 # Explicit Vela false rows are source-side classifications fixed by the field coverage receipt.
 for gid,p in vela.items():
  statuses.append({'galaxy_id':gid,'field_id':'Vela','team_release':'Vela TR1','spatial_shard':'','status':'no_dr10_field_coverage','full_stellar':False,
                   'reason':'complete padded final Vela TR1 source rectangle has zero DR10 Tractor objects; no substitute survey permitted','logMstar':'','response_products_used':False})
 statuses.sort(key=lambda r:(r.get('field_id',''),r['galaxy_id']))
 profiles.sort(key=lambda r:(r['galaxy_id'],float(r['radius_arcsec'])))
 with (OUT/'optical_status.csv').open('w',newline='') as f:
  fields=sorted({k for r in statuses for k in r});w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(statuses)
 with (OUT/'optical_profile.csv').open('w',newline='') as f:
  fields=['galaxy_id','radius_arcsec','mu_r','muerr_r','ellipticity'];w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(profiles)
 audits=[]
 for p in audit_files: audits+=list(csv.DictReader(p.open()))
 if audits:
  fields=sorted({k for r in audits for k in r}); audits.sort(key=lambda r:(r.get('galaxy_id',''),float(r.get('size_pix') or 0)))
  with (OUT/'frame_audit.csv').open('w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(audits)
 # Patch the complete 1760-row H I metadata without changing H I robust_sample.
 meta=list(csv.DictReader(META_HI.open())); sm={r['galaxy_id']:r for r in statuses}
 if len(meta)!=1760: raise SystemExit(f'H I metadata row count {len(meta)} != 1760')
 for r in meta:
  gid=r['galaxy_id']; s=sm.get(gid)
  if s is None:
   r['full_stellar']=False; r['logMstar']=''
  else:
   yes=(s['status']=='pass' and boolish(s['full_stellar'])); r['full_stellar']=yes; r['logMstar']=s.get('logMstar','') if yes else ''
 meta.sort(key=lambda r:(r['field_id'],r['galaxy_id']))
 req=['galaxy_id','field_id','release_phase','distance_mpc','logMstar','logMHI','ell_maj_beams','robust_sample','full_stellar']
 with (OUT/'meta.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=req);w.writeheader();w.writerows([{k:r[k] for k in req} for r in meta])
 # Final source-only counts.
 byfield=defaultdict(Counter)
 for r in statuses: byfield[r['field_id']][r['status']]+=1
 robust={r['galaxy_id'] for r in meta if boolish(r['robust_sample'])}; final_full={r['galaxy_id'] for r in meta if boolish(r['robust_sample']) and boolish(r['full_stellar'])}
 if final_full!=passed: raise SystemExit(f'final full-stellar set mismatch metadata={len(final_full)} optical_pass={len(passed)}')
 report={'schema':'WALLABY_STAGE1_SOURCE_OPTICAL_FULL_v1','scope':'source-only optical processing of committed H I survivors; no response products',
         'n_hi_robust':len(robust),'n_optical_candidates_ngc_fields':len(candidate),'n_vela_hi_survivors_forced_optical_false':len(vela),'n_full_stellar':len(final_full),
         'n_optical_profile_rows':len(profiles),'status_by_field':{k:dict(v) for k,v in sorted(byfield.items())},
         'production_frame_lock_sha256':sha(LOCK),'vela_coverage_receipt_sha256':sha(VELA),'response_products_used':False,
         'next':'Create exact source export manifest for meta.csv, committed H I profile and optical_profile.csv; run byte-verified v4.99 source engine in TARGET mode; freeze predictor/HI-fold/strata/beam receipts before response opening.'}
 (OUT/'SOURCE_OPTICAL_STAGE1_REPORT.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
 files=[OUT/'meta.csv',OUT/'optical_status.csv',OUT/'optical_profile.csv',OUT/'SOURCE_OPTICAL_STAGE1_REPORT.json']
 if (OUT/'frame_audit.csv').exists(): files.append(OUT/'frame_audit.csv')
 with (OUT/'SHA256SUMS.txt').open('w') as f:
  for p in files:f.write(f'{sha(p)}  {p.relative_to(ROOT)}\n')
 print(json.dumps(report,indent=2,sort_keys=True))
if __name__=='__main__':main()
