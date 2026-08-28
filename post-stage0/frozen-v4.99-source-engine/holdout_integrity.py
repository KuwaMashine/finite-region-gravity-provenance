#!/usr/bin/env python3
"""Integrity/timestamp checks shared by the staged WALLABY holdout workflow."""
from __future__ import annotations
import hashlib, json
from pathlib import Path


def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):
            h.update(b)
    return h.hexdigest()


def require_hash(path: Path, expected: str, label: str) -> str:
    p=Path(path)
    if not p.exists():
        raise SystemExit(f"TARGET LOCK: missing {label}: {p}")
    got=sha256(p)
    if not expected or got.lower()!=str(expected).lower():
        raise SystemExit(f"TARGET LOCK: {label} SHA256 mismatch: {got} != {expected}")
    return got


def _load_receipt(path: Path, label: str) -> dict:
    p=Path(path)
    if not p.exists(): raise SystemExit(f"TARGET LOCK: missing {label}: {p}")
    try: obj=json.loads(p.read_text())
    except Exception as e: raise SystemExit(f"TARGET LOCK: invalid {label}: {e}")
    for k in ('timestamp_service','timestamp_utc','proof_reference'):
        if not str(obj.get(k,'')).strip(): raise SystemExit(f"TARGET LOCK: {label} missing {k}")
    return obj


def verify_initial_receipt(receipt_path: Path, manifest_path: Path) -> dict:
    r=_load_receipt(receipt_path,'initial timestamp receipt')
    m=Path(manifest_path)
    if not m.exists(): raise SystemExit(f"TARGET LOCK: pretarget manifest missing: {m}")
    mh=sha256(m)
    if str(r.get('frozen_manifest_file_sha256','')).lower()!=mh.lower():
        raise SystemExit('TARGET LOCK: initial timestamp receipt does not bind this pretarget manifest')
    return r


def verify_strata_receipt(receipt_path: Path, strata_path: Path, manifest_path: Path, beam_report_path: Path|None=None) -> dict:
    r=_load_receipt(receipt_path,'proxy-strata timestamp receipt')
    s=Path(strata_path); m=Path(manifest_path)
    sh=sha256(s); mh=sha256(m)
    if str(r.get('proxy_strata_map_sha256','')).lower()!=sh.lower():
        raise SystemExit('TARGET LOCK: proxy-strata receipt does not bind this strata map')
    if str(r.get('frozen_manifest_file_sha256','')).lower()!=mh.lower():
        raise SystemExit('TARGET LOCK: proxy-strata receipt was made against a different pretarget manifest')
    if beam_report_path is not None:
        bh=sha256(Path(beam_report_path))
        if str(r.get('source_beam_report_sha256','')).lower()!=bh.lower():
            raise SystemExit('TARGET LOCK: proxy-strata receipt does not bind the source-only beam report')
    return r


def verify_adapter_lock(lock_path: Path, base_dir: Path) -> dict:
    p=Path(lock_path)
    try: lock=json.loads(p.read_text())
    except Exception as e: raise SystemExit(f"TARGET LOCK: invalid adapter lock: {e}")
    comps=lock.get('components',{})
    if not comps: raise SystemExit('TARGET LOCK: adapter lock has no components')
    base=Path(base_dir)
    for rel,expected in sorted(comps.items()):
        require_hash((base/rel).resolve(),expected,f'adapter component {rel}')
    return lock


def verify_export_manifest(manifest_path: Path, actual_files: dict[str, Path], stage: str) -> dict:
    """Verify the per-run public-product -> normalized-export provenance manifest.

    This does not inspect scientific values. It binds the frozen adapter to exact
    product identifiers, column mappings and normalized file hashes for the run.
    """
    p=Path(manifest_path)
    try: obj=json.loads(p.read_text())
    except Exception as e: raise SystemExit(f"TARGET LOCK: invalid {stage} export manifest: {e}")
    if str(obj.get('stage','')).upper()!=stage.upper():
        raise SystemExit(f"TARGET LOCK: export manifest stage {obj.get('stage')} != {stage}")
    if str(obj.get('join_key',''))!='galaxy_id':
        raise SystemExit('TARGET LOCK: export manifest join_key must be galaxy_id')
    if not isinstance(obj.get('product_ids'),dict) or not obj['product_ids']:
        raise SystemExit('TARGET LOCK: export manifest has no public product identifiers')
    if not isinstance(obj.get('column_map'),dict) or not obj['column_map']:
        raise SystemExit('TARGET LOCK: export manifest has no raw-to-normalized column map')
    nf=obj.get('normalized_files',{})
    for key,path in actual_files.items():
        item=nf.get(key)
        if not isinstance(item,dict) or not item.get('sha256'):
            raise SystemExit(f'TARGET LOCK: export manifest missing normalized_files.{key}.sha256')
        require_hash(Path(path),str(item['sha256']),f'{stage} normalized file {key}')
    if stage.upper()=='SOURCE_ONLY':
        if obj.get('response_products_loaded') is not False:
            raise SystemExit('TARGET LOCK: source export manifest must state response_products_loaded=false')
        if int(obj.get('validation_release_phase',-1))!=2:
            raise SystemExit('TARGET LOCK: source export manifest validation_release_phase must be 2')
        vf=set(obj.get('validation_fields',[]))
        if vf!={'NGC5044','NGC4808','Vela'}:
            raise SystemExit('TARGET LOCK: source export manifest validation_fields mismatch')
    elif stage.upper()=='RESPONSE':
        if obj.get('source_stage_timestamp_verified') is not True:
            raise SystemExit('TARGET LOCK: response export manifest must state source_stage_timestamp_verified=true')
    return obj
