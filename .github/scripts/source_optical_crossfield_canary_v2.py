#!/usr/bin/env python3
import csv
import importlib.util
import json
import math
import re
import shutil
from pathlib import Path

import numpy as np
from astropy.io import fits
from dl import authClient as ac

ROOT=Path.cwd()
BASE=ROOT/'.github/scripts/source_optical_crossfield_canary.py'
spec=importlib.util.spec_from_file_location('base_canary',BASE)
base=importlib.util.module_from_spec(spec); spec.loader.exec_module(base)

SIZE=1024
OUT=ROOT/'post-stage0/source-optical-crossfield-canary-v2'
TMP=Path('/tmp/wallaby_source_optical_crossfield_canary_v2')
MARGIN_PIX=2.0


def fitted_center(aux):
    m=re.search(r'center x:\s*([0-9.+-]+) pix, y:\s*([0-9.+-]+) pix',Path(aux).read_text(errors='replace'))
    if not m: raise base.QualityFail('unable to read fitted r-band center')
    return tuple(map(float,m.groups()))


def process(source,token):
    name=source['name']; release=source['team_release']; ra=float(source['ra']); dec=float(source['dec'])
    slug=base.slugify(name); work=TMP/slug; sout=OUT/slug
    work.mkdir(parents=True,exist_ok=True); sout.mkdir(parents=True,exist_ok=True)
    result={
        'scope':'source-only cross-field optical canary v2; no validation kinematics queried',
        'source_name':name,'team_release':release,'ra':ra,'dec':dec,
        'distance_mpc':float(source['dist_h']),'log_m_hi_corr':float(source['log_m_hi_corr']),
        'ell_maj_beams_30arcsec':float(source['ell_maj_beams_30arcsec']),'qflag':float(source['qflag']),
        'cutout_size_pix':SIZE,'cutout_side_arcsec':SIZE*base.PIX_SCALE,
        'containment_rule':'common aperture semi-major radius must be <= orientation-independent inscribed radius from fitted r-band center minus 2 pixels',
        'retry_policy':base.RETRY_POLICY,'full_stellar':None,
    }
    try:
        count,body=base.datalab_count(token,ra,dec)
        (sout/'dr10_tractor_count.csv').write_text(body if body.endswith('\n') else body+'\n')
        result['dr10_tractor_count_r0p02deg']=count
    except Exception as e:
        result['classification_status']='unresolved_catalog_query_failure'; result['error']=str(e); return result
    if count==0:
        result['classification_status']='no_dr10_catalog_coverage'; result['full_stellar']=False; return result

    cut=f'https://www.legacysurvey.org/viewer/cutout.fits?ra={ra}&dec={dec}&layer=ls-dr10&pixscale={base.PIX_SCALE}&bands=grz&size={SIZE}'
    dust=f'https://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust?locstr={ra}%20{dec}'
    try:
        base.curl_retry(cut,work/'dr10_grz.fits'); base.curl_retry(dust,work/'dust.xml')
        result['dr10_fits_sha256']=base.sha256(work/'dr10_grz.fits'); result['dust_xml_sha256']=base.sha256(work/'dust.xml')
    except base.NetworkUnresolved as e:
        result['classification_status']='unresolved_network_after_retries'; result['error']=str(e); return result

    try:
        with fits.open(work/'dr10_grz.fits',memmap=False) as h:
            d=np.asarray(h[0].data); bands=str(h[0].header.get('BANDS',''))
        if d.shape!=(3,SIZE,SIZE) or bands!='grz': raise base.QualityFail(f'unexpected DR10 cube shape={d.shape} BANDS={bands!r}')
        if not np.isfinite(d).all(): raise base.QualityFail('nonfinite DR10 pixels')
        for i,b in enumerate('grz'): fits.PrimaryHDU(d[i].astype('float32')).writeto(work/f'{b}.fits',overwrite=True)

        cfg=work/'r_config.py'
        cfg.write_text(
            "ap_process_mode='image'\n"
            f"ap_image_file=r'{work/'r.fits'}'\n"
            f"ap_name='{slug}_r'\n"
            f"ap_pixscale={base.PIX_SCALE}\n"
            f"ap_zeropoint={base.ZEROPOINT}\n"
            "ap_doplot=False\n"
            "ap_isoclip=True\n"
            f"ap_guess_center={{'x': {SIZE/2:.1f}, 'y': {SIZE/2:.1f}}}\n"
        )
        base.run_autoprof(cfg,work)
        rp=work/f'{slug}_r.prof'; ra_=work/f'{slug}_r.aux'
        if not rp.exists() or not ra_.exists(): raise base.QualityFail('missing r-band AutoProf outputs')
        shutil.copy2(rp,work/'r.prof'); shutil.copy2(ra_,work/'r.aux')
        rl=work/f'{slug}_r.log'
        if rl.exists(): shutil.copy2(rl,work/'r.log')

        for b in ('g','z'):
            cfg=work/f'{b}_config.py'
            cfg.write_text(
                "ap_process_mode='forced image'\n"
                f"ap_image_file=r'{work/f'{b}.fits'}'\n"
                f"ap_name='{slug}_{b}'\n"
                f"ap_pixscale={base.PIX_SCALE}\n"
                f"ap_zeropoint={base.ZEROPOINT}\n"
                "ap_doplot=False\n"
                "ap_isoclip=True\n"
                f"ap_forcing_profile=r'{work/'r.prof'}'\n"
            )
            base.run_autoprof(cfg,work)
            pp=work/f'{slug}_{b}.prof'; aa=work/f'{slug}_{b}.aux'
            if not pp.exists() or not aa.exists(): raise base.QualityFail(f'missing forced {b} outputs')
            shutil.copy2(pp,work/f'{b}.prof'); shutil.copy2(aa,work/f'{b}.aux')
            ll=work/f'{slug}_{b}.log'
            if ll.exists(): shutil.copy2(ll,work/f'{b}.log')

        ebv=base.dust_ebv(work/'dust.xml')
        details=base.full_stellar_from_profiles(source,work,sout,ebv)
        cx,cy=fitted_center(work/'r.aux')
        safe_pix=min(cx,(SIZE-1)-cx,cy,(SIZE-1)-cy)-MARGIN_PIX
        safe_arcsec=safe_pix*base.PIX_SCALE
        result['fitted_center_pix']={'x':cx,'y':cy}
        result['orientation_independent_safe_radius_arcsec']=safe_arcsec
        result.update(details)
        if details['common_aperture_radius_arcsec']>safe_arcsec:
            result['classification_status']='unresolved_frame_limited'
            result['full_stellar']=None
            result['error']=f"common aperture {details['common_aperture_radius_arcsec']:.6f} arcsec exceeds safe radius {safe_arcsec:.6f} arcsec"
        else:
            result['classification_status']='processed_full_stellar_true'
            result['full_stellar']=True
        for fn in ('r.prof','r.aux','r.log','g.prof','g.aux','g.log','z.prof','z.aux','z.log'):
            p=work/fn
            if p.exists(): shutil.copy2(p,sout/fn)
        return result
    except base.QualityFail as e:
        result['classification_status']='processed_full_stellar_false'; result['full_stellar']=False; result['error']=str(e); return result
    except Exception as e:
        result['classification_status']='unresolved_processing_exception'; result['full_stellar']=None; result['error']=f'{type(e).__name__}: {e}'; return result


def main():
    if TMP.exists(): shutil.rmtree(TMP)
    TMP.mkdir(parents=True); OUT.mkdir(parents=True,exist_ok=True)
    selected=base.select_canaries()
    with (OUT/'selected_canaries.csv').open('w',newline='') as f:
        fields=['name','team_release','ra','dec','dist_h','log_m_hi_corr','ell_maj_beams_30arcsec','rel','qflag','selection_distance_from_2p5_beams']
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader()
        for r in selected: w.writerow({k:r[k] for k in fields})
    token=ac.login('anonymous')
    results=[]
    for s in selected:
        print(f"processing {s['team_release']} :: {s['name']} :: beams={s['ell_maj_beams_30arcsec']}")
        r=process(s,token); results.append(r)
        so=OUT/base.slugify(s['name']); (so/'result.json').write_text(json.dumps(r,indent=2,sort_keys=True)+'\n')
        print(json.dumps(r,indent=2,sort_keys=True))
    old=json.loads((ROOT/'post-stage0/source-optical-crossfield-canary/crossfield_canary_summary.json').read_text())
    oldmap={r['source_name']:r for r in old['results']}
    comparisons={}
    for r in results:
        o=oldmap.get(r['source_name'])
        if o and r.get('full_stellar') is True and o.get('full_stellar') is True:
            comparisons[r['source_name']]={
                'delta_R50_arcsec':r['R50_r_arcsec']-o['R50_r_arcsec'],
                'delta_log10_Mstar':r['log10_Mstar_adopted_median']-o['log10_Mstar_adopted_median'],
            }
    summary={
        'scope':'source-only deterministic three-field 1024-pixel canary v2; no validation kinematics queried',
        'selection_rule':'within each frozen final release: qflag==0 and ell_maj_beams>=2; choose source closest to 2.5 beams, tie by source name',
        'cutout_size_pix':SIZE,'pixscale_arcsec':base.PIX_SCALE,'containment_margin_pix':MARGIN_PIX,
        'n_canaries':len(results),'n_full_stellar_true':sum(r['full_stellar'] is True for r in results),
        'n_full_stellar_false':sum(r['full_stellar'] is False for r in results),'n_unresolved':sum(r['full_stellar'] is None for r in results),
        'v1_to_v2_comparisons':comparisons,'results':results,
        'note':'AutoProf checkfit flags remain diagnostics only. Frame-limited and network-exhausted cases are unresolved, not astrophysical false classifications.'
    }
    (OUT/'crossfield_canary_v2_summary.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
    files=sorted(p for p in OUT.rglob('*') if p.is_file() and p.name!='SHA256SUMS.txt')
    with (OUT/'SHA256SUMS.txt').open('w') as f:
        for p in files: f.write(f'{base.sha256(p)}  {p.relative_to(ROOT)}\n')
    print(json.dumps(summary,indent=2,sort_keys=True))

if __name__=='__main__': main()
