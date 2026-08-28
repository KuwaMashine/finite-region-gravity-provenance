#!/usr/bin/env python3
"""Process one spatial shard of Stage-1 source-only DR10 optical photometry.

Selection comes only from the committed source-only H I Stage-1 survivor plan.
The production method is frozen in SOURCE_OPTICAL_PRODUCTION_FRAME_LOCK.json.
No WALLABY kinematic product or response column is read.
"""
from __future__ import annotations
import argparse, copy, csv, hashlib, importlib.util, json, math, re, shutil, subprocess, time
from pathlib import Path
import numpy as np
import requests
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from dl import authClient as ac
from dl import queryClient as qc
from reproject import reproject_interp

ROOT=Path.cwd()
PLAN=ROOT/'post-stage0/source-optical-stage1-preplan/optical_stage1_plan.csv'
LOCK=ROOT/'post-stage0/SOURCE_OPTICAL_PRODUCTION_FRAME_LOCK.json'
BASE_SCRIPT=ROOT/'.github/scripts/source_optical_crossfield_canary.py'
spec=importlib.util.spec_from_file_location('base_canary',BASE_SCRIPT); base=importlib.util.module_from_spec(spec); spec.loader.exec_module(base)
NERSC='https://portal.nersc.gov/cfs/cosmo/data/legacysurvey/dr10/south/coadd'

class NoCoverage(RuntimeError): pass
class NetworkUnresolved(RuntimeError): pass


def ffloat(x):
    v=float(x)
    if not math.isfinite(v): raise ValueError(f'nonfinite {x!r}')
    return v

def sha256(p:Path):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''): h.update(b)
    return h.hexdigest()

def download(url: str, dest: Path):
    if dest.exists() and dest.stat().st_size>0: return
    last=None
    for attempt in range(4):
        try:
            with requests.get(url,stream=True,timeout=(20,240),headers={'User-Agent':'finite-region-gravity-stage1-optical/1'}) as r:
                if r.status_code==404: raise NoCoverage(f'HTTP 404 {url}')
                if r.status_code>=400: raise NetworkUnresolved(f'HTTP {r.status_code} {url}')
                with dest.open('wb') as f:
                    for chunk in r.iter_content(1<<20):
                        if chunk: f.write(chunk)
            if dest.stat().st_size<=0: raise NetworkUnresolved(f'zero-byte {url}')
            return
        except NoCoverage: raise
        except Exception as e:
            last=e
            try: dest.unlink()
            except FileNotFoundError: pass
            if attempt<3: time.sleep(2**attempt)
    raise NetworkUnresolved(str(last))

def load_brick(cache:Path, brick:str, band:str):
    p=cache/f'{brick}_{band}.fits.fz'; sub=brick[:3]
    url=f'{NERSC}/{sub}/{brick}/legacysurvey-{brick}-image-{band}.fits.fz'
    download(url,p)
    with fits.open(p,memmap=False) as h:
        hh=h[1] if len(h)>1 and h[1].data is not None else h[0]
        if hh.data is None: raise NoCoverage(f'no image data {url}')
        return np.asarray(hh.data).copy(), WCS(hh.header).celestial, url, sha256(p)

def query_bricks(token, tw, size):
    pts=np.array([[0,0],[size-1,0],[0,size-1],[size-1,size-1],[(size-1)/2,0],[(size-1)/2,size-1],[0,(size-1)/2],[size-1,(size-1)/2]],float)
    ra,dec=tw.pixel_to_world_values(pts[:,0],pts[:,1]); ra=np.asarray(ra,float); dec=np.asarray(dec,float)
    if np.nanmax(ra)-np.nanmin(ra)>180: raise RuntimeError('RA wrap unsupported')
    ramin=float(np.nanmin(ra))-1e-4; ramax=float(np.nanmax(ra))+1e-4; dmin=float(np.nanmin(dec))-1e-4; dmax=float(np.nanmax(dec))+1e-4
    sql=("SELECT brickname,ra,dec,ra1,ra2,dec1,dec2 FROM ls_dr10.bricks "
         f"WHERE ra2>{ramin} AND ra1<{ramax} AND dec2>{dmin} AND dec1<{dmax} ORDER BY brickname")
    body=qc.query(token=token,sql=sql,fmt='csv'); rows=list(csv.DictReader(body.splitlines()))
    if not rows: raise NetworkUnresolved(f'no brick geometry response for footprint; query={sql}')
    for b in rows:
        for k in ('ra','dec','ra1','ra2','dec1','dec2'): b[k]=float(b[k])
    return rows,sql

def target_wcs(primary_wcs, coord, size):
    sx,sy=primary_wcs.world_to_pixel(coord); sx=float(np.asarray(sx).reshape(-1)[0]); sy=float(np.asarray(sy).reshape(-1)[0])
    x0=int(round(sx))-size//2; y0=int(round(sy))-size//2
    tw=copy.deepcopy(primary_wcs); tw.wcs.crpix[0]-=x0; tw.wcs.crpix[1]-=y0
    tx,ty=tw.world_to_pixel(coord)
    return tw,x0,y0,float(np.asarray(tx).reshape(-1)[0]),float(np.asarray(ty).reshape(-1)[0])

def build_mosaics(cache, work, primary_brick, primary_data, primary_wcs, coord, size, token):
    tw,x0,y0,tx,ty=target_wcs(primary_wcs,coord,size)
    yy,xx=np.indices((size,size),dtype=float); tra,tdec=tw.pixel_to_world_values(xx,yy)
    bricks,sql=query_bricks(token,tw,size)
    owners=[]; count=np.zeros((size,size),np.int16)
    for b in bricks:
        own=(tra>=b['ra1'])&(tra<b['ra2'])&(tdec>=b['dec1'])&(tdec<b['dec2']); owners.append(own); count+=own.astype(np.int16)
    if not np.all(count==1):
        vals,cnts=np.unique(count,return_counts=True); raise RuntimeError(f'nonunique brick owner map {dict(zip(vals.tolist(),cnts.tolist()))}')
    brecord={}
    for band in 'grz':
        out=np.full((size,size),np.nan,np.float32); rec=[]
        for b,own in zip(bricks,owners):
            if not np.any(own): continue
            brick=b['brickname']
            if brick==primary_brick:
                data,sw,url,hsh=primary_data[band]
                srcx=xx.astype(np.int64)+x0; srcy=yy.astype(np.int64)+y0
                valid=own&(srcx>=0)&(srcx<data.shape[1])&(srcy>=0)&(srcy<data.shape[0])
                if int(valid.sum())!=int(own.sum()): raise RuntimeError(f'{band}: primary owner outside primary array')
                out[own]=data[srcy[own],srcx[own]].astype(np.float32); mode='exact_native_integer_pixels'
            else:
                data,sw,url,hsh=load_brick(cache,brick,band)
                repro,foot=reproject_interp((data,sw),tw,shape_out=(size,size),order='bilinear',return_footprint=True)
                bad=own&((foot<=0)|(~np.isfinite(repro)))
                if np.any(bad): raise NoCoverage(f'{brick} {band}: {int(bad.sum())} geometrically owned pixels have no finite DR10 image')
                out[own]=repro[own].astype(np.float32); mode='bilinear_neighbor_only'
            rec.append({'brick':brick,'band':band,'mode':mode,'owned_pixels':int(own.sum()),'url':url,'sha256':hsh})
        if not np.isfinite(out).all(): raise NoCoverage(f'{band}: nonfinite final mosaic')
        fits.PrimaryHDU(out,header=tw.to_header()).writeto(work/f'{band}.fits',overwrite=True)
        brecord[band]=rec
    return tw,tx,ty,bricks,sql,brecord

def parse_aux(path):
    txt=Path(path).read_text(errors='replace')
    def one(pat):
        m=re.search(pat,txt); return None if not m else float(m.group(1))
    cm=re.search(r'center x:\s*([-+0-9.eE]+) pix, y:\s*([-+0-9.eE]+) pix',txt)
    return {'center_x_pix':None if not cm else float(cm.group(1)),'center_y_pix':None if not cm else float(cm.group(2)),
            'background_flux_per_pix':one(r'background:\s*([-+0-9.eE]+)'),
            'background_noise_flux_per_pix':one(r'noise:\s*([-+0-9.eE]+) flux/pix')}

def common_rows(work):
    P={b:base.read_prof(work/f'{b}.prof') for b in 'grz'}; n=min(len(P[b]) for b in 'grz')
    if n<5: raise base.QualityFail(f'fewer than five common rows: {n}')
    keep=n
    for i in range(n):
        rs=[float(P[b][i]['R']) for b in 'grz']
        if max(rs)-min(rs)>1e-8: raise base.QualityFail(f'forced radius mismatch row {i}')
        es=[float(P[b][i]['SB_e']) for b in 'grz']
        if any((not math.isfinite(e)) or e>=base.COMMON_ERR for e in es): keep=i; break
    if keep<5: raise base.QualityFail(f'common uncertainty aperture retains only {keep} rows')
    return {b:P[b][:keep] for b in 'grz'}

def stage1_r50(rows):
    vals=[]
    for z in rows:
        r=float(z['R']); mu=float(z['SB']); emu=float(z['SB_e']); ell=float(z['ellip'])
        if all(math.isfinite(x) for x in (r,mu,emu,ell)): vals.append((r,mu,emu,ell))
    vals.sort()
    if len(vals)<5 or vals[0][0]<=0 or any(vals[i+1][0]<=vals[i][0] for i in range(len(vals)-1)): raise base.QualityFail('invalid normalized optical profile')
    r=[x[0] for x in vals]; mu=[x[1] for x in vals]; q=[1.0-x[3] for x in vals]
    if any(x<=0 or x>1 for x in q): raise base.QualityFail('invalid optical axis ratio')
    mid=[0.5*(r[i]+r[i+1]) for i in range(len(r)-1)]; rin=[0.0]*len(r); rout=[0.0]*len(r); rin[1:]=mid; rout[:-1]=mid
    rin[0]=max(0.0,r[0]-0.5*(r[1]-r[0])); rout[-1]=r[-1]+0.5*(r[-1]-r[-2]); m0=min(mu)
    flux=[10**(-0.4*(m-m0))*math.pi*qq*(ro*ro-ri*ri) for m,qq,ri,ro in zip(mu,q,rin,rout)]
    total=sum(flux); target=0.5*total; cum=0.0
    for j,f in enumerate(flux):
        if cum+f>=target:
            frac=(target-cum)/f; return math.sqrt(rin[j]**2+frac*(rout[j]**2-rin[j]**2))
        cum+=f
    raise base.QualityFail('Stage-1 half-light crossing absent')

def run_frame(row, cache, source_tmp, primary_brick, primary_data, primary_wcs, coord, size, token, ebv):
    work=source_tmp/f'n{size}'; work.mkdir(parents=True,exist_ok=True)
    tw,tx,ty,bricks,brick_query,brec=build_mosaics(cache,work,primary_brick,primary_data,primary_wcs,coord,size,token)
    slug=base.slugify(row['galaxy_id'])
    cfg=work/'r_config.py'; cfg.write_text("ap_process_mode='image'\n"+f"ap_image_file=r'{work/'r.fits'}'\n"+f"ap_name='{slug}_r_{size}'\n"+
        f"ap_pixscale={base.PIX_SCALE}\nap_zeropoint={base.ZEROPOINT}\nap_doplot=False\nap_isoclip=True\n"+f"ap_guess_center={{'x': {tx!r}, 'y': {ty!r}}}\n")
    base.run_autoprof(cfg,work); rp=work/f'{slug}_r_{size}.prof'; ra=work/f'{slug}_r_{size}.aux'; shutil.copy2(rp,work/'r.prof'); shutil.copy2(ra,work/'r.aux')
    rl=work/f'{slug}_r_{size}.log';
    if rl.exists(): shutil.copy2(rl,work/'r.log')
    for band in ('g','z'):
        cfg=work/f'{band}_config.py'; cfg.write_text("ap_process_mode='forced image'\n"+f"ap_image_file=r'{work/f'{band}.fits'}'\n"+f"ap_name='{slug}_{band}_{size}'\n"+
            f"ap_pixscale={base.PIX_SCALE}\nap_zeropoint={base.ZEROPOINT}\nap_doplot=False\nap_isoclip=True\n"+f"ap_forcing_profile=r'{work/'r.prof'}'\n")
        base.run_autoprof(cfg,work); pp=work/f'{slug}_{band}_{size}.prof'; aa=work/f'{slug}_{band}_{size}.aux'; shutil.copy2(pp,work/f'{band}.prof'); shutil.copy2(aa,work/f'{band}.aux')
        ll=work/f'{slug}_{band}_{size}.log';
        if ll.exists(): shutil.copy2(ll,work/f'{band}.log')
    sink=work/'derived'; sink.mkdir(exist_ok=True)
    src={'dist_h':row['distance_mpc']}; details=base.full_stellar_from_profiles(src,work,sink,ebv)
    P=common_rows(work); r50=stage1_r50(P['r']); aux=parse_aux(work/'r.aux'); cx=aux['center_x_pix']; cy=aux['center_y_pix']
    if None in (cx,cy,aux['background_flux_per_pix'],aux['background_noise_flux_per_pix']): raise base.QualityFail('cannot parse required r-band nuisance metadata')
    safe=(min(cx,(size-1)-cx,cy,(size-1)-cy)-2.0)*base.PIX_SCALE; contained=float(P['r'][-1]['R'])<=safe
    norm=[{'galaxy_id':row['galaxy_id'],'radius_arcsec':float(z['R']),'mu_r':float(z['SB']),'muerr_r':float(z['SB_e']),'ellipticity':float(z['ellip'])} for z in P['r']]
    return {'size_pix':size,'contained':bool(contained),'safe_radius_arcsec':safe,'common_aperture_arcsec':float(P['r'][-1]['R']),
            'stage1_R50_arcsec':r50,'stage1_Rd_arcsec':r50/1.6783469900166605,'logMstar':details['log10_Mstar_adopted_median'],
            'logMstar_method_sd':details['log10_Mstar_method_population_sd'],'n_retained':len(P['r']),'r_background':aux['background_flux_per_pix'],
            'r_background_rms':aux['background_noise_flux_per_pix'],'checkfit_json':json.dumps(details['autoprof_checkfit'],sort_keys=True),
            'n_intersecting_bricks':len(bricks),'brick_query':brick_query,'brick_records_json':json.dumps(brec,sort_keys=True),'norm':norm}

def pair_pass(a,b,lock):
    c=lock['convergence_criteria']; rms=0.5*(abs(a['r_background_rms'])+abs(b['r_background_rms']))
    bg=math.inf if not rms>0 else abs(b['r_background']-a['r_background'])/rms
    metrics={'abs_delta_stage1_R50_arcsec':abs(b['stage1_R50_arcsec']-a['stage1_R50_arcsec']),
             'abs_delta_log10_Mstar_dex':abs(b['logMstar']-a['logMstar']),
             'abs_delta_common_aperture_arcsec':abs(b['common_aperture_arcsec']-a['common_aperture_arcsec']),
             'abs_delta_background_over_mean_RMS':bg}
    ok=(a['contained'] and b['contained'] and metrics['abs_delta_stage1_R50_arcsec']<=c['max_abs_delta_stage1_R50_arcsec'] and
        metrics['abs_delta_log10_Mstar_dex']<=c['max_abs_delta_log10_Mstar_dex'] and metrics['abs_delta_common_aperture_arcsec']<=c['max_abs_delta_common_aperture_arcsec'] and
        metrics['abs_delta_background_over_mean_RMS']<=c['max_abs_delta_background_over_mean_RMS'])
    return ok,metrics

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--shard-index',type=int,required=True); ap.add_argument('--out-dir',type=Path,required=True); a=ap.parse_args()
    lock=json.loads(LOCK.read_text()); sizes=lock['sizes_pix']; a.out_dir.mkdir(parents=True,exist_ok=True)
    rows=[r for r in csv.DictReader(PLAN.open()) if r['dr10_optical_field_status']=='candidate' and int(r['spatial_shard'])==a.shard_index]
    cache=Path('/tmp')/f'wallaby_stage1_optical_cache_{a.shard_index:02d}'; cache.mkdir(exist_ok=True)
    token=ac.login('anonymous'); statuses=[]; profiles=[]; audits=[]; unresolved=[]
    for idx,row in enumerate(rows,1):
        gid=row['galaxy_id']; rec={'galaxy_id':gid,'field_id':row['field_id'],'team_release':row['team_release'],'spatial_shard':a.shard_index,'full_stellar':False,'status':'pending','response_products_used':False}
        st=Path('/tmp')/f'optical_{a.shard_index:02d}_{hashlib.sha256(gid.encode()).hexdigest()[:12]}'; shutil.rmtree(st,ignore_errors=True); st.mkdir()
        try:
            coord=SkyCoord(float(row['ra'])*u.deg,float(row['dec'])*u.deg); pb=row['primary_brick']; primary={}
            for band in 'grz': primary[band]=load_brick(cache,pb,band)
            primary_wcs=primary['r'][1]
            dust=f'https://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust?locstr={float(row["ra"])}%20{float(row["dec"])}'; dp=st/'dust.xml'
            try: base.curl_retry(dust,dp); ebv=base.dust_ebv(dp)
            except base.NetworkUnresolved as e: raise NetworkUnresolved(str(e))
            frames=[]; chosen=None; previous_success=None
            for size in sizes:
                try:
                    fr=run_frame(row,cache,st,pb,primary,primary_wcs,coord,int(size),token,ebv); frames.append(fr)
                    audits.append({k:v for k,v in fr.items() if k not in ('norm','brick_records_json','brick_query')} | {'galaxy_id':gid,'frame_status':'success','reason':''})
                    if previous_success is not None:
                        ok,metrics=pair_pass(previous_success,fr,lock)
                        for ar in reversed(audits):
                            if ar['galaxy_id']==gid and ar['size_pix']==previous_success['size_pix']:
                                ar.update({f'pair_to_next_{k}':v for k,v in metrics.items()}); ar['pair_to_next_pass']=ok; break
                        if ok: chosen=previous_success; break
                    previous_success=fr
                except NoCoverage as e:
                    rec.update({'status':'no_dr10_imaging_coverage','reason':str(e)}); chosen='terminal_false'; break
                except base.QualityFail as e:
                    audits.append({'galaxy_id':gid,'size_pix':size,'frame_status':'quality_fail','reason':str(e)})
                    previous_success=None
                except base.NetworkUnresolved as e: raise NetworkUnresolved(str(e))
            if chosen=='terminal_false': pass
            elif chosen is None:
                rec.update({'status':'frame_nonconverged_after_registered_ladder','reason':'no adjacent registered frame pair satisfied all frozen convergence criteria'})
            else:
                rec.update({'status':'pass','full_stellar':True,'logMstar':chosen['logMstar'],'logMstar_method_sd':chosen['logMstar_method_sd'],
                            'selected_size_pix':chosen['size_pix'],'stage1_R50_arcsec':chosen['stage1_R50_arcsec'],'stage1_Rd_arcsec':chosen['stage1_Rd_arcsec'],
                            'common_aperture_arcsec':chosen['common_aperture_arcsec'],'n_optical_profile_rows':len(chosen['norm']),'SFD_EBV':ebv,'reason':''})
                profiles.extend(chosen['norm'])
        except NoCoverage as e: rec.update({'status':'no_dr10_imaging_coverage','reason':str(e)})
        except NetworkUnresolved as e:
            rec.update({'status':'unresolved_network','reason':str(e)}); unresolved.append(gid)
        except base.QualityFail as e: rec.update({'status':'optical_quality_fail','reason':str(e)})
        except Exception as e: rec.update({'status':'implementation_failure','reason':f'{type(e).__name__}: {e}'})
        statuses.append(rec); shutil.rmtree(st,ignore_errors=True)
        print(f'shard {a.shard_index}: {idx}/{len(rows)} {gid} {rec["status"]}',flush=True)
    sf=['galaxy_id','field_id','team_release','spatial_shard','status','full_stellar','reason','logMstar','logMstar_method_sd','selected_size_pix','stage1_R50_arcsec','stage1_Rd_arcsec','common_aperture_arcsec','n_optical_profile_rows','SFD_EBV','response_products_used']
    with (a.out_dir/'status.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=sf,extrasaction='ignore');w.writeheader();w.writerows(statuses)
    pf=['galaxy_id','radius_arcsec','mu_r','muerr_r','ellipticity']
    with (a.out_dir/'optical_profile.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=pf);w.writeheader();w.writerows(profiles)
    af=sorted({k for r in audits for k in r.keys()}) if audits else ['galaxy_id','size_pix','frame_status','reason']
    with (a.out_dir/'frame_audit.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=af,extrasaction='ignore');w.writeheader();w.writerows(audits)
    counts={s:sum(r['status']==s for r in statuses) for s in sorted({r['status'] for r in statuses})}
    summary={'shard_index':a.shard_index,'n_candidates':len(rows),'status_counts':counts,'n_full_stellar':sum(bool(r['full_stellar']) for r in statuses),
             'n_unresolved_network':len(unresolved),'unresolved_network_ids':unresolved,'frame_lock_sha256':sha256(LOCK),'response_products_used':False}
    (a.out_dir/'summary.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n'); print(json.dumps(summary,sort_keys=True))
    if unresolved: raise SystemExit(f'unresolved network cases remain: {len(unresolved)}')
    impl=[r for r in statuses if r['status']=='implementation_failure']
    if impl: raise SystemExit(f'implementation failures remain: {len(impl)}')
if __name__=='__main__': main()
