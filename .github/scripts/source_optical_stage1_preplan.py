#!/usr/bin/env python3
"""Resolve DR10 primary bricks and spatial shards for Stage-1 H I survivors.

This is source-only planning. It reads the committed H I Stage-1 status and final
source catalogue, plus the public DR10 brick geometry table. It never reads a
WALLABY kinematic product or response column.
"""
from __future__ import annotations
import csv, hashlib, json, math
from pathlib import Path
from collections import Counter
from dl import authClient as ac
from dl import queryClient as qc

ROOT=Path.cwd(); OUT=ROOT/'post-stage0/source-optical-stage1-preplan'; OUT.mkdir(parents=True,exist_ok=True)
HI=ROOT/'post-stage0/source-hi-stage1/source_hi_status.csv'
CAT=ROOT/'post-stage0/source-only-census/phase2_source_ids.csv'
VELA=ROOT/'post-stage0/source-vela-dr10-field-coverage/VELA_DR10_FIELD_COVERAGE.json'
FIELDS={'NGC4808':'NGC 4808 TR1','NGC5044':'NGC 5044 TR3'}
SHARDS=24

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def shard(brick): return int.from_bytes(hashlib.sha256(brick.encode()).digest()[:8],'big')%SHARDS

def main():
    hs=list(csv.DictReader(HI.open())); src=list(csv.DictReader(CAT.open()))
    sm={(r['name'],r['team_release']):r for r in src}
    passes=[r for r in hs if r['status']=='pass']
    if not passes: raise SystemExit('no committed H I survivors')
    vela=json.loads(VELA.read_text())
    if not vela.get('all_vela_positions_outside_dr10_catalog_footprint'): raise SystemExit('Vela field coverage receipt not in expected locked state')
    token=ac.login('anonymous')
    plan=[]; field_bricks={}
    for field,release in FIELDS.items():
        rr=[r for r in passes if r['field_id']==field]
        if not rr: continue
        coords=[]
        for r in rr:
            z=sm.get((r['galaxy_id'],release))
            if z is None: raise SystemExit(f'missing source row {r["galaxy_id"]} {release}')
            coords.append((r,z,float(z['ra']),float(z['dec'])))
        ramin=min(x[2] for x in coords)-0.5; ramax=max(x[2] for x in coords)+0.5
        dmin=min(x[3] for x in coords)-0.5; dmax=max(x[3] for x in coords)+0.5
        if ramax-ramin>180: raise SystemExit(f'RA wrap unsupported for {field}')
        sql=("SELECT brickname,ra,dec,ra1,ra2,dec1,dec2 FROM ls_dr10.bricks "
             f"WHERE ra2>{ramin} AND ra1<{ramax} AND dec2>{dmin} AND dec1<{dmax} ORDER BY brickname")
        text=qc.query(token=token,sql=sql,fmt='csv'); bricks=list(csv.DictReader(text.splitlines()))
        if not bricks: raise SystemExit(f'no DR10 bricks returned for {field}')
        for b in bricks:
            for k in ('ra','dec','ra1','ra2','dec1','dec2'): b[k]=float(b[k])
        field_bricks[field]={'query':sql,'n_bricks':len(bricks),'bounds':[ramin,ramax,dmin,dmax]}
        for h,z,ra,dec in coords:
            hit=[b for b in bricks if ra>=b['ra1'] and ra<b['ra2'] and dec>=b['dec1'] and dec<b['dec2']]
            if len(hit)!=1: raise SystemExit(f'{field} {h["galaxy_id"]}: containing brick count {len(hit)}')
            b=hit[0]; brick=b['brickname']
            plan.append({'galaxy_id':h['galaxy_id'],'field_id':field,'team_release':release,'ra':ra,'dec':dec,
                         'distance_mpc':z['dist_h'],'logMHI':z['log_m_hi_corr'],'ell_maj_beams':z['ell_maj_beams_30arcsec'],
                         'primary_brick':brick,'primary_brick_subdir':brick[:3],'spatial_shard':shard(brick),
                         'hi_status':'pass','dr10_optical_field_status':'candidate'})
    # Carry every Vela H I survivor as a deterministic optical-false row; no brick query or substitution.
    for h in passes:
        if h['field_id']!='Vela': continue
        z=sm.get((h['galaxy_id'],'Vela TR1'))
        if z is None: raise SystemExit(f'missing Vela source row {h["galaxy_id"]}')
        plan.append({'galaxy_id':h['galaxy_id'],'field_id':'Vela','team_release':'Vela TR1','ra':z['ra'],'dec':z['dec'],
                     'distance_mpc':z['dist_h'],'logMHI':z['log_m_hi_corr'],'ell_maj_beams':z['ell_maj_beams_30arcsec'],
                     'primary_brick':'','primary_brick_subdir':'','spatial_shard':'','hi_status':'pass',
                     'dr10_optical_field_status':'no_dr10_field_coverage'})
    plan.sort(key=lambda r:(r['field_id'],str(r['spatial_shard']),r['primary_brick'],r['galaxy_id']))
    fields=list(plan[0].keys())
    with (OUT/'optical_stage1_plan.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(plan)
    c=Counter(r['field_id'] for r in plan); cs=Counter(r['dr10_optical_field_status'] for r in plan)
    report={'schema':'WALLABY_STAGE1_OPTICAL_PREPLAN_v1','scope':'committed H I survivors plus public DR10 brick geometry only; no response products',
            'n_hi_survivors':len(passes),'n_plan_rows':len(plan),'by_field':dict(sorted(c.items())),'optical_field_status':dict(sorted(cs.items())),
            'n_spatial_shards':SHARDS,'shard_rule':'SHA256(primary_brick) first 64 bits mod 24; all sources in a primary brick stay together',
            'field_brick_queries':field_bricks,'vela_coverage_receipt':str(VELA.relative_to(ROOT)),
            'frame_lock':'post-stage0/SOURCE_OPTICAL_PRODUCTION_FRAME_LOCK.json','response_products_used':False}
    (OUT/'OPTICAL_STAGE1_PREPLAN_REPORT.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    with (OUT/'SHA256SUMS.txt').open('w') as f:
        for p in (OUT/'optical_stage1_plan.csv',OUT/'OPTICAL_STAGE1_PREPLAN_REPORT.json'): f.write(f'{sha(p)}  {p.relative_to(ROOT)}\n')
    print(json.dumps(report,indent=2,sort_keys=True))
if __name__=='__main__': main()
