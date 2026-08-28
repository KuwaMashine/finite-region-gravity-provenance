#!/usr/bin/env python3
"""Prospective source-only field-level DR10 coverage audit for final Vela TR1."""
import csv, json, math
from pathlib import Path
from dl import authClient as ac
from dl import queryClient as qc

ROOT=Path.cwd(); OUT=ROOT/'post-stage0/source-vela-dr10-field-coverage'; OUT.mkdir(parents=True,exist_ok=True)
rows=[r for r in csv.DictReader((ROOT/'post-stage0/source-only-census/phase2_source_ids.csv').open()) if r['team_release']=='Vela TR1']
if len(rows)!=203: raise SystemExit(f'Vela final-release row count changed: {len(rows)}')
ras=[float(r['ra']) for r in rows]; decs=[float(r['dec']) for r in rows]
if max(ras)-min(ras)>180: raise SystemExit('RA wrap not supported by this field audit')
pad=0.05; bounds={'ra_min':min(ras)-pad,'ra_max':max(ras)+pad,'dec_min':min(decs)-pad,'dec_max':max(decs)+pad}
token=ac.login('anonymous')
sql=("SELECT COUNT(*) AS n FROM ls_dr10.tractor_s WHERE "
     f"ra>={bounds['ra_min']} AND ra<={bounds['ra_max']} AND dec>={bounds['dec_min']} AND dec<={bounds['dec_max']}")
text=qc.query(token=token,sql=sql,fmt='csv')
rr=list(csv.DictReader(text.splitlines()))
if len(rr)!=1: raise SystemExit(f'unexpected Data Lab response {text[:500]!r}')
n=int(rr[0]['n'])
# If the bounding rectangle is nonempty, count every Vela source locally to distinguish partial footprint from full coverage.
per=[]
if n>0:
    for i,r in enumerate(rows,1):
        ra=float(r['ra']); dec=float(r['dec'])
        q=f"SELECT COUNT(*) AS n FROM ls_dr10.tractor_s WHERE q3c_radial_query(ra,dec,{ra},{dec},0.02)"
        t=qc.query(token=token,sql=q,fmt='csv'); z=list(csv.DictReader(t.splitlines()))
        if len(z)!=1: raise SystemExit(f'bad local count for {r["name"]}')
        per.append({'galaxy_id':r['name'],'team_release':'Vela TR1','ra':ra,'dec':dec,'tractor_count_r0p02deg':int(z[0]['n'])})
        if i%25==0: print(f'{i}/{len(rows)}',flush=True)
    with (OUT/'vela_source_local_counts.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(per[0].keys()));w.writeheader();w.writerows(per)
report={'schema':'WALLABY_VELA_DR10_FIELD_COVERAGE_v1','scope':'source-only final Vela TR1 sky positions; no WALLABY kinematic product or response column queried',
        'n_vela_final_release_sources':len(rows),'bounds_deg':bounds,'dr10_tractor_count_in_full_padded_bounding_rectangle':n,
        'all_vela_positions_outside_dr10_catalog_footprint':(n==0),
        'n_sources_with_local_dr10_catalog_coverage':(0 if n==0 else sum(x['tractor_count_r0p02deg']>0 for x in per)),
        'query':sql,'response_products_used':False,
        'interpretation':('The complete padded bounding rectangle containing all final Vela TR1 source positions has zero DR10 Tractor objects. Under the frozen DR10 all-three-band optical gate, Vela has no production optical catalogue footprint; no substitute survey is permitted.' if n==0 else 'DR10 coverage is not uniformly absent across the Vela rectangle; use the recorded per-source local counts for source-side optical eligibility.')}
(OUT/'VELA_DR10_FIELD_COVERAGE.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
print(json.dumps(report,indent=2,sort_keys=True))
