#!/usr/bin/env python3
"""Consolidate all deterministic H I source shards into the Stage-1 source export."""
from __future__ import annotations
import argparse,csv,hashlib,json,math
from collections import Counter,defaultdict
from pathlib import Path
import pandas as pd

ROOT=Path.cwd(); CAT=ROOT/'post-stage0/source-only-probe/source_catalogue_safe.csv'
OUT=ROOT/'post-stage0/source-hi-stage1'
RELEASES={'NGC 4808 TR1':'NGC4808','NGC 5044 TR3':'NGC5044','Vela TR1':'Vela'}
REQ_FINITE=('ra','dec','dist_h','log_m_hi_corr','ell_maj','ell_min','ell_pa')

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): h.update(b)
 return h.hexdigest()

def eligible(r):
 try:
  if float(r['qflag'])!=0.0:return False
  return all(math.isfinite(float(r[c])) for c in REQ_FINITE)
 except:return False

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--shards-dir',type=Path,required=True);a=ap.parse_args();OUT.mkdir(parents=True,exist_ok=True)
 cat=pd.read_csv(CAT,dtype=str); final=cat[cat.team_release.isin(RELEASES)].copy()
 if len(final)!=1760: raise SystemExit(f'final-release count {len(final)} != 1760')
 expected={str(r['name']) for _,r in final.iterrows() if eligible(r)}
 status_files=sorted(a.shards_dir.rglob('status.csv')); profile_files=sorted(a.shards_dir.rglob('profile.csv'))
 if not status_files or not profile_files: raise SystemExit('missing shard artifacts')
 statuses=[]
 for p in status_files: statuses += list(csv.DictReader(p.open()))
 ids=[r['galaxy_id'] for r in statuses]
 dup=[k for k,v in Counter(ids).items() if v!=1]
 if dup: raise SystemExit(f'duplicate shard status ids: {dup[:8]}')
 got=set(ids)
 if got!=expected: raise SystemExit(f'shard coverage mismatch missing={len(expected-got)} extra={len(got-expected)}')
 smap={r['galaxy_id']:r for r in statuses}
 profiles=[]
 for p in profile_files: profiles += list(csv.DictReader(p.open()))
 pc=Counter(r['galaxy_id'] for r in profiles)
 pass_ids={r['galaxy_id'] for r in statuses if r['status']=='pass'}
 badprof=[g for g in pass_ids if pc[g]<8]
 extra_prof=set(pc)-pass_ids
 if badprof or extra_prof: raise SystemExit(f'profile consistency fail short={badprof[:8]} extra={list(extra_prof)[:8]}')
 statuses=sorted(statuses,key=lambda r:(r['field_id'],r['galaxy_id']))
 profiles=sorted(profiles,key=lambda r:(r['galaxy_id'],float(r['rad_hi_source_arcsec'])))
 with (OUT/'source_hi_status.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(statuses[0].keys()));w.writeheader();w.writerows(statuses)
 with (OUT/'hi_profile.csv').open('w',newline='') as f:
  fields=['galaxy_id','rad_hi_source_arcsec','sigma_hi_source_weight'];w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(profiles)
 meta=[]
 for _,r in final.iterrows():
  gid=str(r['name']); ok=eligible(r) and gid in smap and smap[gid]['status']=='pass'; field=RELEASES[str(r['team_release'])]
  def val(c):
   try:return float(r[c])
   except:return ''
  ev=val('ell_maj')
  meta.append({'galaxy_id':gid,'field_id':field,'release_phase':2,'distance_mpc':val('dist_h'),'logMstar':'','logMHI':val('log_m_hi_corr'),
               'ell_maj_beams':(ev/5.0 if ev!='' else ''),'robust_sample':bool(ok),'full_stellar':False,'team_release':str(r['team_release'])})
 meta=sorted(meta,key=lambda r:(r['field_id'],r['galaxy_id']))
 with (OUT/'source_meta_hi.csv').open('w',newline='') as f:
  fields=list(meta[0].keys());w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(meta)
 counts=defaultdict(Counter)
 for r in statuses: counts[r['field_id']][r['status']]+=1
 report={'schema':'WALLABY_STAGE1_SOURCE_HI_FULL_v1','scope':'all finite qflag=0 rows in frozen final Phase-2 source releases; source_data only; no response products',
         'n_final_release_rows':int(len(final)),'n_q0_finite_candidates':len(expected),'n_pass':len(pass_ids),'n_fail':len(expected)-len(pass_ids),
         'field_status_counts':{k:dict(v) for k,v in sorted(counts.items())},'n_profile_rows':len(profiles),
         'frozen_profiler_sha256':'eb290a85d5d4600a03602d5ab3a4c7c8ae1a7e40d18c74cfdee70f2b10804306',
         'source_engine_verification':'post-stage0/frozen-v4.99-source-engine/VERIFICATION.json','response_products_used':False,
         'next':'run optical/stellar adapter only for robust_sample=true galaxies; do not open validation kinematics'}
 (OUT/'SOURCE_HI_STAGE1_REPORT.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
 files=[OUT/'source_hi_status.csv',OUT/'hi_profile.csv',OUT/'source_meta_hi.csv',OUT/'SOURCE_HI_STAGE1_REPORT.json']
 with (OUT/'SHA256SUMS.txt').open('w') as f:
  for p in files:f.write(f'{sha(p)}  {p.relative_to(ROOT)}\n')
 print(json.dumps(report,indent=2,sort_keys=True))
if __name__=='__main__':main()
