#!/usr/bin/env python3
"""Repair layer for the TARGET source-stage preflight.

The first preflight assumed one CAOM2 source-plane spelling. This wrapper probes
only a closed list of plausible *source_data* product IDs for the already frozen
NGC 5044 DR3/TR3 release. It never lists all planes for the observation and thus
cannot discover or inspect a kinematic product by accident.
"""
from __future__ import annotations
import csv, importlib.util, urllib.parse, urllib.request
from pathlib import Path

BASE_PATH=Path(__file__).with_name('source_stage_target_preflight.py')
spec=importlib.util.spec_from_file_location('source_stage_target_preflight_base',BASE_PATH)
base=importlib.util.module_from_spec(spec); spec.loader.exec_module(base)

CANDIDATES=(
    'source_data_NGC5044_TR3',
    'source_data_NGC_5044_TR3',
    'source_data_NGC5044_DR3',
    'source_data_NGC_5044_DR3',
    'source_data_N5044_TR3',
    'source_data_N5044_DR3',
)

def source_artifacts_closed_probe():
    attempts=[]
    for plane in CANDIDATES:
        q=("SELECT o.observationID,p.productID,a.uri,a.productType,a.contentType,a.contentLength,a.contentChecksum "
           "FROM caom2.Observation AS o JOIN caom2.Plane AS p ON o.obsID=p.obsID "
           "JOIN caom2.Artifact AS a ON p.planeID=a.planeID "
           f"WHERE o.collection='WALLABY' AND o.observationID='{base.OBS}' AND p.productID='{plane}'")
        data=urllib.parse.urlencode({'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'csv','QUERY':q}).encode()
        req=urllib.request.Request(base.TAP,data=data,headers={'User-Agent':'finite-region-gravity-source-preflight/2'})
        with urllib.request.urlopen(req,timeout=180) as r: text=r.read().decode('utf-8-sig')
        rows=list(csv.DictReader(text.splitlines())); attempts.append({'productID':plane,'n_artifacts':len(rows)})
        if not rows: continue
        bad=[r for r in rows if any(x in (r.get('uri','')+' '+r.get('productID','')).lower() for x in base.FORBIDDEN)]
        if bad: raise RuntimeError('source-plane artifact firewall violation')
        if any(r.get('productID')!=plane for r in rows): raise RuntimeError('CADC returned an unrequested productID')
        mom=[r for r in rows if r.get('uri','').endswith('_mom0.fits')]
        mask=[r for r in rows if r.get('uri','').endswith('_mask.fits')]
        if len(mom)!=1 or len(mask)!=1: raise RuntimeError(f'{plane}: source plane lacks unique mom0/mask artifacts')
        prefix=mom[0]['uri'].split('/')[-1][:-len('_mom0.fits')]
        if mask[0]['uri'].split('/')[-1] != prefix+'_mask.fits': raise RuntimeError('mom0/mask source prefixes differ')
        if not prefix.startswith(base.OBS+'_'): raise RuntimeError(f'unexpected source artifact prefix {prefix}')
        base.PLANE=plane; base.OBJ=prefix
        # Preserve the exact successful query plus a source-only naming audit.
        return rows,text,q
    raise RuntimeError('no artifacts under closed source_data productID candidates: '+repr(attempts))

base.tap_source_artifacts=source_artifacts_closed_probe
base.main()
