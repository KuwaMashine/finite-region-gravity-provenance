#!/usr/bin/env python3
"""Resolve final-release WALLABY CAOM2 source_data plane spellings without listing planes.

Only a closed list of plausible source_data product IDs is queried for one known
source observation in each frozen validation release. No kinematic or response
product name is queried or discovered.
"""
from __future__ import annotations
import csv, hashlib, json, urllib.parse, urllib.request
from pathlib import Path

ROOT=Path.cwd(); OUT=ROOT/'post-stage0/source-plane-resolution'; OUT.mkdir(parents=True,exist_ok=True)
TAP='https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/argus/sync'
CASES={
 'NGC4808':{
   'observationID':'WALLABY_J130150+041953',
   'candidates':['source_data_NGC4808_TR1','source_data_NGC_4808_TR1','source_data_NGC4808_DR1','source_data_NGC_4808_DR1','source_data_N4808_TR1','source_data_N4808_DR1']},
 'NGC5044':{
   'observationID':'WALLABY_J133209-245132',
   'candidates':['source_data_NGC5044_TR3','source_data_NGC_5044_TR3','source_data_NGC5044_DR3','source_data_NGC_5044_DR3','source_data_N5044_TR3','source_data_N5044_DR3']},
 'Vela':{
   'observationID':'WALLABY_J094821-452603',
   'candidates':['source_data_Vela_TR1']},
}
FORBIDDEN=('kinematic','kinmodel','rotcur','model_data','vrot','gobs','gbar','rar_residual','vmax','vflat')

def query(obs,plane):
    q=("SELECT o.observationID,p.productID,a.uri,a.contentType,a.contentLength,a.contentChecksum "
       "FROM caom2.Observation AS o JOIN caom2.Plane AS p ON o.obsID=p.obsID JOIN caom2.Artifact AS a ON p.planeID=a.planeID "
       f"WHERE o.collection='WALLABY' AND o.observationID='{obs}' AND p.productID='{plane}'")
    data=urllib.parse.urlencode({'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'csv','QUERY':q}).encode()
    req=urllib.request.Request(TAP,data=data,headers={'User-Agent':'finite-region-gravity-source-plane-resolver/1'})
    with urllib.request.urlopen(req,timeout=180) as r: text=r.read().decode('utf-8-sig')
    return list(csv.DictReader(text.splitlines())),q

def main():
    result={'schema':'WALLABY_SOURCE_PLANE_RESOLUTION_v1','scope':'closed source_data productID probes only; no validation kinematic plane queried','response_products_used':False,'fields':{}}
    for field,spec in CASES.items():
        attempts=[]; chosen=None
        for plane in spec['candidates']:
            rows,q=query(spec['observationID'],plane); attempts.append({'productID':plane,'n_artifacts':len(rows)})
            if not rows: continue
            if any(r.get('productID')!=plane for r in rows): raise RuntimeError('unrequested productID returned')
            if any(any(x in (r.get('uri','')+' '+r.get('productID','')).lower() for x in FORBIDDEN) for r in rows): raise RuntimeError('source-plane firewall violation')
            mom=[r for r in rows if r.get('uri','').endswith('_mom0.fits')]; mask=[r for r in rows if r.get('uri','').endswith('_mask.fits')]
            if len(mom)!=1 or len(mask)!=1: raise RuntimeError(f'{field} {plane}: lacks unique mom0/mask')
            chosen={'productID':plane,'observationID':spec['observationID'],'n_artifacts':len(rows),'mom0_uri':mom[0]['uri'],'mask_uri':mask[0]['uri'],'query':q}
            break
        if chosen is None: raise RuntimeError(f'no closed source_data candidate resolved for {field}: {attempts}')
        result['fields'][field]={'chosen':chosen,'attempts':attempts}
    p=OUT/'SOURCE_PLANE_RESOLUTION.json'; p.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    h=hashlib.sha256(p.read_bytes()).hexdigest(); (OUT/'SHA256SUMS.txt').write_text(f'{h}  {p.relative_to(ROOT)}\n')
    print(json.dumps(result,indent=2,sort_keys=True))
if __name__=='__main__': main()
