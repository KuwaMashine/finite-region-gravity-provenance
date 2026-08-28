#!/usr/bin/env python3
import csv, hashlib, json, math
from pathlib import Path

ROOT=Path.cwd()
OUT=ROOT/'post-stage0/source-optical-stage1-hookup'
OUT.mkdir(parents=True,exist_ok=True)
ERR=0.22
XHALF=1.6783469900166605
R50_TOL=0.5

CASES={
 'frame_1024': ROOT/'post-stage0/source-optical-nersc-frame-convergence/n1024',
 'frame_1280': ROOT/'post-stage0/source-optical-nersc-frame-convergence/n1280',
 'bilinear_multibrick': ROOT/'post-stage0/source-optical-nersc-multibrick-probe',
 'native_multibrick': ROOT/'post-stage0/source-optical-nersc-native-multibrick-probe',
}

def read_prof(p):
    with p.open() as f:
        f.readline()
        return list(csv.DictReader(f))

def common_rows(d):
    P={b:read_prof(d/f'{b}.prof') for b in 'grz'}
    n=min(len(P[b]) for b in 'grz')
    if n<5: raise RuntimeError(f'{d}: fewer than five common rows')
    keep=n
    for i in range(n):
        rs=[float(P[b][i]['R']) for b in 'grz']
        if max(rs)-min(rs)>1e-8: raise RuntimeError(f'{d}: radius mismatch row {i}')
        es=[float(P[b][i]['SB_e']) for b in 'grz']
        if any((not math.isfinite(x)) or x>=ERR for x in es):
            keep=i; break
    if keep<5: raise RuntimeError(f'{d}: common aperture leaves {keep} rows')
    return {b:P[b][:keep] for b in 'grz'},keep

def edges(r):
    mid=[0.5*(r[i]+r[i+1]) for i in range(len(r)-1)]
    rin=[0.0]*len(r); rout=[0.0]*len(r)
    rin[1:]=mid; rout[:-1]=mid
    rin[0]=max(0.0,r[0]-0.5*(r[1]-r[0]))
    rout[-1]=r[-1]+0.5*(r[-1]-r[-2])
    return rin,rout

def stage1_r50(rows):
    vals=[]
    for z in rows:
        r=float(z['R']); mu=float(z['SB']); emu=float(z['SB_e']); ell=float(z['ellip'])
        if all(math.isfinite(x) for x in (r,mu,emu,ell)):
            vals.append((r,mu,emu,ell))
    vals.sort()
    if len(vals)<5 or vals[0][0]<=0 or any(vals[i+1][0]<=vals[i][0] for i in range(len(vals)-1)):
        raise RuntimeError('invalid normalized optical profile')
    r=[x[0] for x in vals]; mu=[x[1] for x in vals]; ell=[x[3] for x in vals]
    q=[1.0-x for x in ell]
    if any(x<=0 or x>1 for x in q): raise RuntimeError('invalid ellipticity')
    rin,rout=edges(r); m0=min(mu)
    flux=[10**(-0.4*(m-m0))*math.pi*qq*(ro*ro-ri*ri) for m,qq,ri,ro in zip(mu,q,rin,rout)]
    total=sum(flux); target=0.5*total; cum=0.0
    for j,f in enumerate(flux):
        if cum+f>=target:
            frac=(target-cum)/f
            return math.sqrt(rin[j]**2+frac*(rout[j]**2-rin[j]**2))
        cum+=f
    raise RuntimeError('half-light crossing absent')

def autoprof_totmag_r50(rows):
    rf=[]
    for z in rows:
        r=float(z['R']); m=float(z['totmag'])
        if not (math.isfinite(r) and math.isfinite(m) and m<99): raise RuntimeError('bad totmag row')
        rf.append((r,10**(-0.4*m)))
    target=0.5*rf[-1][1]
    k=next(i for i,(_,f) in enumerate(rf) if f>=target)
    if k==0:return rf[0][0]
    r0,f0=rf[k-1]; r1,f1=rf[k]
    if not f1>f0: raise RuntimeError('non-increasing totmag curve')
    return r0+(target-f0)/(f1-f0)*(r1-r0)

results={}
for name,d in CASES.items():
    P,n=common_rows(d)
    r50=stage1_r50(P['r']); auto=autoprof_totmag_r50(P['r'])
    results[name]={
      'directory':str(d.relative_to(ROOT)), 'n_common_retained':n,
      'common_aperture_radius_arcsec':float(P['r'][-1]['R']),
      'stage1_engine_R50_r_arcsec':r50,
      'stage1_engine_Rd_star_arcsec':r50/XHALF,
      'autoprof_totmag_R50_r_arcsec':auto,
      'stage1_minus_autoprof_R50_arcsec':r50-auto,
    }

comparisons={
 'frame_1024_to_1280':{
   'a':'frame_1024','b':'frame_1280',
   'abs_delta_stage1_R50_arcsec':abs(results['frame_1280']['stage1_engine_R50_r_arcsec']-results['frame_1024']['stage1_engine_R50_r_arcsec'])},
 'bilinear_to_native_multibrick':{
   'a':'bilinear_multibrick','b':'native_multibrick',
   'abs_delta_stage1_R50_arcsec':abs(results['native_multibrick']['stage1_engine_R50_r_arcsec']-results['bilinear_multibrick']['stage1_engine_R50_r_arcsec'])},
}
for x in comparisons.values(): x['passes_0p5_arcsec']=x['abs_delta_stage1_R50_arcsec']<=R50_TOL
status='pass' if all(x['passes_0p5_arcsec'] for x in comparisons.values()) else 'fail'
report={
 'scope':'source-only optical-to-Stage-1 estimator hookup; no validation kinematics queried',
 'normalization_rule':'retain grz profile rows strictly interior to first row where any band SB_e is nonfinite or >=0.22; export retained r R,SB,SB_e,ellip as normalized optical profile',
 'stage1_estimator':'frozen v4.99 optical_exponential_equivalent_scale: relative r-band SB integrated over local elliptical annuli; Rd=R50/1.6783469900166605',
 'comparison_tolerance_R50_arcsec':R50_TOL,
 'results':results,'comparisons':comparisons,'status':status,
 'decision':('The normalized common-aperture r profile is a stable input to the frozen Stage-1 estimator across the tested frame and resampling perturbations.' if status=='pass' else 'Do not launch production optical processing: frozen Stage-1 R50 is not stable under an already registered source-side perturbation.'),
}
p=OUT/'stage1_optical_hookup.json'; p.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
with (OUT/'SHA256SUMS.txt').open('w') as f: f.write(f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(ROOT)}\n")
print(json.dumps(report,indent=2,sort_keys=True))
if status!='pass': raise SystemExit('Stage-1 optical hookup stability failed')
