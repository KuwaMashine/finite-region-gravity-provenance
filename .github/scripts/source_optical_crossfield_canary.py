#!/usr/bin/env python3
import csv
import hashlib
import json
import math
import re
import shutil
import statistics
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
from astropy.io import fits
from dl import authClient as ac
from dl import queryClient as qc

ROOT = Path.cwd()
OUT = ROOT / 'post-stage0' / 'source-optical-crossfield-canary'
TMP = Path('/tmp/wallaby_source_optical_crossfield_canary')
CENSUS = ROOT / 'post-stage0' / 'source-only-census' / 'phase2_source_ids.csv'
COEFF = ROOT / 'post-stage0' / 'STAGE0A_MLCR_COEFFICIENTS.csv'
RELEASES = ('NGC 4808 TR1', 'NGC 5044 TR3', 'Vela TR1')
TARGET_BEAMS = 2.5
PIX_SCALE = 0.262
ZEROPOINT = 22.5
COMMON_ERR = 0.22
RETRY_POLICY = {
    'attempts': 4,
    'delay_seconds': 5,
    'connect_timeout_seconds': 20,
    'max_time_seconds': 180,
}

class NetworkUnresolved(RuntimeError):
    pass

class QualityFail(RuntimeError):
    pass


def finite_float(x):
    try:
        v = float(x)
        return v if math.isfinite(v) else None
    except Exception:
        return None


def slugify(name):
    return re.sub(r'[^A-Za-z0-9]+', '_', name).strip('_').lower()


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def curl_retry(url, dest):
    dest = Path(dest)
    cmd = [
        'curl', '-fSsL', '--retry', str(RETRY_POLICY['attempts'] - 1),
        '--retry-delay', str(RETRY_POLICY['delay_seconds']), '--retry-all-errors',
        '--connect-timeout', str(RETRY_POLICY['connect_timeout_seconds']),
        '--max-time', str(RETRY_POLICY['max_time_seconds']), '-o', str(dest), url,
    ]
    cp = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if cp.returncode != 0 or not dest.exists() or dest.stat().st_size == 0:
        raise NetworkUnresolved(f'curl failed rc={cp.returncode}: {cp.stderr[-1000:]}')
    return cp


def select_canaries():
    rows = list(csv.DictReader(CENSUS.open(encoding='utf-8')))
    chosen = []
    for release in RELEASES:
        eligible = []
        for r in rows:
            if r['team_release'] != release:
                continue
            q = finite_float(r['qflag'])
            b = finite_float(r['ell_maj_beams_30arcsec'])
            D = finite_float(r['dist_h'])
            mhi = finite_float(r['log_m_hi_corr'])
            ra = finite_float(r['ra'])
            dec = finite_float(r['dec'])
            if q != 0.0 or b is None or b < 2.0:
                continue
            if None in (D, mhi, ra, dec):
                continue
            eligible.append((abs(b - TARGET_BEAMS), r['name'], r))
        if not eligible:
            raise SystemExit(f'no eligible canary in {release}')
        eligible.sort(key=lambda x: (x[0], x[1]))
        r = dict(eligible[0][2])
        r['selection_distance_from_2p5_beams'] = eligible[0][0]
        chosen.append(r)
    return chosen


def datalab_count(token, ra, dec):
    sql = f'SELECT COUNT(*) AS n FROM ls_dr10.tractor_s WHERE q3c_radial_query(ra,dec,{ra},{dec},0.02)'
    body = qc.query(token=token, sql=sql, fmt='csv')
    rr = list(csv.DictReader(body.splitlines()))
    if len(rr) != 1 or 'n' not in rr[0] or not re.fullmatch(r'\d+', str(rr[0]['n']).strip()):
        raise NetworkUnresolved(f'unexpected Data Lab response: {body[:500]!r}')
    return int(rr[0]['n']), body


def run_autoprof(config_path, cwd):
    cp = subprocess.run(['autoprof', str(config_path)], cwd=cwd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if cp.returncode != 0:
        raise QualityFail(f'AutoProf failed rc={cp.returncode}: {cp.stdout[-2000:]}')
    return cp.stdout


def read_prof(path):
    with Path(path).open() as f:
        f.readline()
        return list(csv.DictReader(f))


def parse_aux_checks(path):
    checks = {}
    if not Path(path).exists():
        return checks
    for line in Path(path).read_text(errors='replace').splitlines():
        if line.startswith('checkfit ') and ': ' in line:
            key, val = line[len('checkfit '):].rsplit(': ', 1)
            checks[key] = (val.strip().lower() == 'pass')
    return checks


def dust_ebv(path):
    root = ET.parse(path).getroot()
    vals = []
    for elem in root.iter():
        if elem.tag.endswith('refPixelValueSFD') and elem.text:
            m = re.search(r'[-+0-9.eE]+', elem.text.strip())
            if m:
                vals.append(float(m.group(0)))
    if not vals or not math.isfinite(vals[0]):
        raise QualityFail('no finite SFD E(B-V)')
    return vals[0]


def full_stellar_from_profiles(source, work, source_out, ebv):
    P = {b: read_prof(work / f'{b}.prof') for b in 'grz'}
    n = min(map(len, P.values()))
    if n < 5:
        raise QualityFail(f'fewer than five common profile rows: {n}')
    first_bad = None
    first_bad_err = None
    for i in range(n):
        rs = [float(P[b][i]['R']) for b in 'grz']
        if max(rs) - min(rs) > 1e-8:
            raise QualityFail(f'forced radius mismatch at row {i}: {rs}')
        errs = {b: float(P[b][i]['SB_e']) for b in 'grz'}
        if any((not math.isfinite(v)) or v >= COMMON_ERR for v in errs.values()):
            first_bad = i
            first_bad_err = errs
            break
    retained = n if first_bad is None else first_bad
    if retained < 5:
        raise QualityFail(f'common uncertainty aperture retains only {retained} rows')
    j = retained - 1
    Rcut = float(P['r'][j]['R'])

    rvals = []
    for i in range(retained):
        R = float(P['r'][i]['R'])
        m = float(P['r'][i]['totmag'])
        if not (math.isfinite(R) and math.isfinite(m) and m < 99):
            raise QualityFail(f'invalid retained r cumulative row {i}')
        rvals.append((R, 10 ** (-0.4 * m)))
    Ftot = rvals[-1][1]
    target = 0.5 * Ftot
    k = next((i for i, (_, F) in enumerate(rvals) if F >= target), None)
    if k is None:
        raise QualityFail('half-light crossing absent')
    if k == 0:
        R50 = rvals[0][0]
    else:
        R0, F0 = rvals[k - 1]
        R1, F1 = rvals[k]
        if not F1 > F0:
            raise QualityFail('non-increasing cumulative light at half-light bracket')
        R50 = R0 + (target - F0) / (F1 - F0) * (R1 - R0)
    if not math.isfinite(R50):
        raise QualityFail('nonfinite R50')

    D = float(source['dist_h'])
    Rd_arcsec = R50 / 1.6783469900166605
    arcsec_to_kpc = D * 1000.0 * math.pi / (180.0 * 3600.0)
    R50_kpc = R50 * arcsec_to_kpc
    Rd_kpc = Rd_arcsec * arcsec_to_kpc

    extcoeff = {'g': 3.214, 'r': 2.165, 'z': 1.211}
    A = {b: extcoeff[b] * ebv for b in extcoeff}
    mobs = {b: float(P[b][j]['totmag']) for b in 'grz'}
    if any(not math.isfinite(v) for v in mobs.values()):
        raise QualityFail('nonfinite aperture magnitude')
    mcorr = {b: mobs[b] - A[b] for b in 'grz'}
    DM = 5.0 * math.log10(D) + 25.0
    Mabs = {b: mcorr[b] - DM for b in 'grz'}
    Msun = {'g': 5.05, 'r': 4.61, 'z': 4.50}
    logL = {b: -0.4 * (Mabs[b] - Msun[b]) for b in 'grz'}
    colours = {
        'g-r': mcorr['g'] - mcorr['r'],
        'g-z': mcorr['g'] - mcorr['z'],
        'r-z': mcorr['r'] - mcorr['z'],
    }

    coeff = list(csv.DictReader(COEFF.open()))
    if len(coeff) != 45:
        raise SystemExit(f'frozen MLCR coefficient count changed: {len(coeff)}')
    est = []
    for row in coeff:
        c = row['color']
        b = row['luminosity_band']
        lm = logL[b] + float(row['intercept_a']) + float(row['slope_b']) * colours[c]
        if not math.isfinite(lm):
            raise QualityFail('nonfinite MLCR estimate')
        est.append({**row, 'color_value_mag': colours[c], 'log10_L_band_Lsun': logL[b], 'log10_Mstar_Msun': lm})
    masses = [x['log10_Mstar_Msun'] for x in est]
    with (source_out / 'mlcr_45.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(est[0].keys()))
        w.writeheader()
        w.writerows(est)

    checks = {b: parse_aux_checks(work / f'{b}.aux') for b in 'grz'}
    return {
        'n_common_profile_rows': n,
        'n_retained': retained,
        'common_aperture_radius_arcsec': Rcut,
        'first_bad_radius_arcsec': None if first_bad is None else float(P['r'][first_bad]['R']),
        'first_bad_SB_e': first_bad_err,
        'R50_r_arcsec': R50,
        'R50_r_kpc': R50_kpc,
        'Rd_star_arcsec': Rd_arcsec,
        'Rd_star_kpc': Rd_kpc,
        'SFD_EBV': ebv,
        'A_mag': A,
        'observed_common_aperture_mag': mobs,
        'foreground_corrected_common_aperture_mag': mcorr,
        'colours_mag': colours,
        'log10_Mstar_adopted_median': statistics.median(masses),
        'log10_Mstar_method_population_sd': statistics.pstdev(masses),
        'min_log10_Mstar': min(masses),
        'max_log10_Mstar': max(masses),
        'n_mlcr_estimates': len(masses),
        'autoprof_checkfit': checks,
    }


def process_source(source, token):
    name = source['name']
    release = source['team_release']
    ra = float(source['ra'])
    dec = float(source['dec'])
    slug = slugify(name)
    work = TMP / slug
    source_out = OUT / slug
    work.mkdir(parents=True, exist_ok=True)
    source_out.mkdir(parents=True, exist_ok=True)

    result = {
        'scope': 'source-only cross-field optical canary; no validation kinematics queried',
        'source_name': name,
        'team_release': release,
        'ra': ra,
        'dec': dec,
        'distance_mpc': float(source['dist_h']),
        'log_m_hi_corr': float(source['log_m_hi_corr']),
        'ell_maj_beams_30arcsec': float(source['ell_maj_beams_30arcsec']),
        'qflag': float(source['qflag']),
        'selection_rule': 'within release: qflag==0 and ell_maj_beams>=2; minimize abs(ell_maj_beams-2.5), tie by source name',
        'retry_policy': RETRY_POLICY,
        'full_stellar': None,
    }

    try:
        count, body = datalab_count(token, ra, dec)
        (source_out / 'dr10_tractor_count.csv').write_text(body if body.endswith('\n') else body + '\n')
        result['dr10_tractor_count_r0p02deg'] = count
    except Exception as e:
        result['classification_status'] = 'unresolved_catalog_query_failure'
        result['error'] = str(e)
        return result

    if count == 0:
        result['classification_status'] = 'no_dr10_catalog_coverage'
        result['full_stellar'] = False
        return result

    cut = f'https://www.legacysurvey.org/viewer/cutout.fits?ra={ra}&dec={dec}&layer=ls-dr10&pixscale={PIX_SCALE}&bands=grz&size=512'
    dust = f'https://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust?locstr={ra}%20{dec}'
    try:
        curl_retry(cut, work / 'dr10_grz.fits')
        curl_retry(dust, work / 'dust.xml')
        result['dr10_fits_sha256'] = sha256(work / 'dr10_grz.fits')
        result['dust_xml_sha256'] = sha256(work / 'dust.xml')
    except NetworkUnresolved as e:
        result['classification_status'] = 'unresolved_network_after_retries'
        result['error'] = str(e)
        return result

    try:
        with fits.open(work / 'dr10_grz.fits', memmap=False) as h:
            d = np.asarray(h[0].data)
            bands = str(h[0].header.get('BANDS', ''))
        if d.ndim != 3 or d.shape[0] != 3 or bands != 'grz':
            raise QualityFail(f'unexpected DR10 cube shape={d.shape} BANDS={bands!r}')
        if not np.isfinite(d).all():
            raise QualityFail('nonfinite DR10 pixels')
        for i, b in enumerate('grz'):
            fits.PrimaryHDU(d[i].astype('float32')).writeto(work / f'{b}.fits', overwrite=True)

        r_config = work / 'r_config.py'
        r_config.write_text(
            "ap_process_mode = 'image'\n"
            f"ap_image_file = r'{work / 'r.fits'}'\n"
            f"ap_name = '{slug}_r'\n"
            f"ap_pixscale = {PIX_SCALE}\n"
            f"ap_zeropoint = {ZEROPOINT}\n"
            "ap_doplot = False\n"
            "ap_isoclip = True\n"
            "ap_guess_center = {'x': 256.0, 'y': 256.0}\n"
        )
        run_autoprof(r_config, work)
        rprof = work / f'{slug}_r.prof'
        raux = work / f'{slug}_r.aux'
        if not rprof.exists() or not raux.exists():
            raise QualityFail('r-band AutoProf outputs missing')
        shutil.copy2(rprof, work / 'r.prof')
        shutil.copy2(raux, work / 'r.aux')
        rlog = work / f'{slug}_r.log'
        if rlog.exists():
            shutil.copy2(rlog, work / 'r.log')

        for b in ('g', 'z'):
            cfg = work / f'{b}_config.py'
            cfg.write_text(
                "ap_process_mode = 'forced image'\n"
                f"ap_image_file = r'{work / f'{b}.fits'}'\n"
                f"ap_name = '{slug}_{b}'\n"
                f"ap_pixscale = {PIX_SCALE}\n"
                f"ap_zeropoint = {ZEROPOINT}\n"
                "ap_doplot = False\n"
                "ap_isoclip = True\n"
                f"ap_forcing_profile = r'{work / 'r.prof'}'\n"
            )
            run_autoprof(cfg, work)
            prof = work / f'{slug}_{b}.prof'
            aux = work / f'{slug}_{b}.aux'
            if not prof.exists() or not aux.exists():
                raise QualityFail(f'{b}-band forced AutoProf outputs missing')
            shutil.copy2(prof, work / f'{b}.prof')
            shutil.copy2(aux, work / f'{b}.aux')
            log = work / f'{slug}_{b}.log'
            if log.exists():
                shutil.copy2(log, work / f'{b}.log')

        ebv = dust_ebv(work / 'dust.xml')
        details = full_stellar_from_profiles(source, work, source_out, ebv)
        result.update(details)
        result['full_stellar'] = True
        result['classification_status'] = 'processed_full_stellar_true'

        for fn in ('r.prof', 'r.aux', 'r.log', 'g.prof', 'g.aux', 'g.log', 'z.prof', 'z.aux', 'z.log'):
            p = work / fn
            if p.exists():
                shutil.copy2(p, source_out / fn)
        return result
    except QualityFail as e:
        result['classification_status'] = 'processed_full_stellar_false'
        result['full_stellar'] = False
        result['error'] = str(e)
        return result
    except Exception as e:
        result['classification_status'] = 'unresolved_processing_exception'
        result['error'] = f'{type(e).__name__}: {e}'
        return result


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    if TMP.exists():
        shutil.rmtree(TMP)
    TMP.mkdir(parents=True)

    selected = select_canaries()
    with (OUT / 'selected_canaries.csv').open('w', newline='') as f:
        fields = ['name', 'team_release', 'ra', 'dec', 'dist_h', 'log_m_hi_corr', 'ell_maj_beams_30arcsec', 'rel', 'qflag', 'selection_distance_from_2p5_beams']
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in selected:
            w.writerow({k: r[k] for k in fields})

    token = ac.login('anonymous')
    results = []
    for source in selected:
        print(f"processing {source['team_release']} :: {source['name']} :: beams={source['ell_maj_beams_30arcsec']}")
        res = process_source(source, token)
        results.append(res)
        source_out = OUT / slugify(source['name'])
        (source_out / 'result.json').write_text(json.dumps(res, indent=2, sort_keys=True) + '\n')
        print(json.dumps(res, indent=2, sort_keys=True))

    summary = {
        'scope': 'source-only deterministic three-field canary; no validation kinematics queried',
        'selection_rule': 'within each frozen final release: qflag==0 and ell_maj_beams>=2; choose source closest to 2.5 beams, tie by source name',
        'target_beams': TARGET_BEAMS,
        'releases': list(RELEASES),
        'retry_policy': RETRY_POLICY,
        'n_canaries': len(results),
        'n_full_stellar_true': sum(r['full_stellar'] is True for r in results),
        'n_full_stellar_false': sum(r['full_stellar'] is False for r in results),
        'n_unresolved': sum(r['full_stellar'] is None for r in results),
        'results': results,
        'note': 'AutoProf checkfit flags are diagnostics only; Light symmetry is not a hard veto. Network exhaustion is unresolved, not full_stellar=false.',
    }
    (OUT / 'crossfield_canary_summary.json').write_text(json.dumps(summary, indent=2, sort_keys=True) + '\n')

    files = sorted(p for p in OUT.rglob('*') if p.is_file() and p.name != 'SHA256SUMS.txt')
    with (OUT / 'SHA256SUMS.txt').open('w') as f:
        for p in files:
            f.write(f'{sha256(p)}  {p.relative_to(ROOT)}\n')
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
