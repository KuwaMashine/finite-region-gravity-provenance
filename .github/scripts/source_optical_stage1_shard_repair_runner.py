#!/usr/bin/env python3
"""Mechanical repair runner for the frozen Stage-1 optical shard processor.

This wrapper does not change source selection, frame sizes, convergence thresholds,
photometric rules, MLCR coefficients, or any response-side firewall. It applies two
execution-only guards already used by source_optical_crossfield_canary.py: when
AutoProf returns success but does not emit the expected .prof/.aux pair, that frame
is classified as QualityFail so the registered frame ladder may continue.
"""
from pathlib import Path

SOURCE = Path('.github/scripts/source_optical_stage1_shard.py')
src = SOURCE.read_text(encoding='utf-8')

old_r = "base.run_autoprof(cfg,work); rp=work/f'{slug}_r_{size}.prof'; ra=work/f'{slug}_r_{size}.aux'; shutil.copy2(rp,work/'r.prof'); shutil.copy2(ra,work/'r.aux')"
new_r = "base.run_autoprof(cfg,work); rp=work/f'{slug}_r_{size}.prof'; ra=work/f'{slug}_r_{size}.aux';\n    if not rp.exists() or not ra.exists(): raise base.QualityFail('r-band AutoProf outputs missing');\n    shutil.copy2(rp,work/'r.prof'); shutil.copy2(ra,work/'r.aux')"

old_forced = "base.run_autoprof(cfg,work); pp=work/f'{slug}_{band}_{size}.prof'; aa=work/f'{slug}_{band}_{size}.aux'; shutil.copy2(pp,work/f'{band}.prof'); shutil.copy2(aa,work/f'{band}.aux')"
new_forced = "base.run_autoprof(cfg,work); pp=work/f'{slug}_{band}_{size}.prof'; aa=work/f'{slug}_{band}_{size}.aux';\n        if not pp.exists() or not aa.exists(): raise base.QualityFail(f'{band}-band forced AutoProf outputs missing');\n        shutil.copy2(pp,work/f'{band}.prof'); shutil.copy2(aa,work/f'{band}.aux')"

if src.count(old_r) != 1:
    raise SystemExit(f'mechanical repair refused: r-band target count={src.count(old_r)}')
if src.count(old_forced) != 1:
    raise SystemExit(f'mechanical repair refused: forced-band target count={src.count(old_forced)}')

src = src.replace(old_r, new_r).replace(old_forced, new_forced)
ns = {'__name__': '__main__', '__file__': str(SOURCE)}
exec(compile(src, str(SOURCE), 'exec'), ns, ns)
