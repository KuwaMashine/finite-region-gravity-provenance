#!/usr/bin/env python3
import json
from pathlib import Path

ROOT=Path.cwd()
src=ROOT/'post-stage0/source-optical-nersc-frame-convergence/nersc_frame_convergence_summary.json'
outdir=ROOT/'post-stage0/source-optical-nersc-frame-convergence-v2'
outdir.mkdir(parents=True,exist_ok=True)
d=json.loads(src.read_text())
res=d['results']
criteria={
  'max_abs_delta_log10_Mstar_dex':0.02,
  'max_abs_delta_R50_arcsec':0.5,
  'max_abs_delta_common_aperture_arcsec':0.5,
  'max_abs_delta_background_over_mean_RMS':0.01,
  'containment_required':True,
}
steps=[]; chosen=None
for a,b in zip(res[:-1],res[1:]):
    db=b['r_aux']['background_flux_per_pix']-a['r_aux']['background_flux_per_pix']
    mean_noise=0.5*(abs(a['r_aux']['background_noise_flux_per_pix'])+abs(b['r_aux']['background_noise_flux_per_pix']))
    row={
      'from_size_pix':a['cutout_size_pix'],'to_size_pix':b['cutout_size_pix'],
      'from_contained':bool(a['common_aperture_fully_contained']),
      'abs_delta_log10_Mstar_dex':abs(b['log10_Mstar_adopted_median']-a['log10_Mstar_adopted_median']),
      'abs_delta_R50_arcsec':abs(b['R50_r_arcsec']-a['R50_r_arcsec']),
      'abs_delta_common_aperture_arcsec':abs(b['common_aperture_radius_arcsec']-a['common_aperture_radius_arcsec']),
      'abs_delta_background_flux_per_pix':abs(db),
      'mean_background_RMS_flux_per_pix':mean_noise,
      'abs_delta_background_over_mean_RMS':abs(db)/mean_noise,
    }
    row['passes']=(row['from_contained'] and
                   row['abs_delta_log10_Mstar_dex']<=criteria['max_abs_delta_log10_Mstar_dex'] and
                   row['abs_delta_R50_arcsec']<=criteria['max_abs_delta_R50_arcsec'] and
                   row['abs_delta_common_aperture_arcsec']<=criteria['max_abs_delta_common_aperture_arcsec'] and
                   row['abs_delta_background_over_mean_RMS']<=criteria['max_abs_delta_background_over_mean_RMS'])
    steps.append(row)
    if chosen is None and row['passes']: chosen=a['cutout_size_pix']
report={
  'scope':'source-only repair of finite-frame convergence decision; no validation kinematics queried',
  'source_name':d['source_name'],'team_release':d['team_release'],'input_summary':str(src.relative_to(ROOT)),
  'previous_decision_rule':d['decision_rule'],'previous_smallest_converged_size_pix':d['smallest_converged_size_pix'],
  'repair_reason':'The previous scalar-only rule admitted 512 pixels even though the fitted background and common uncertainty aperture changed strongly at 768 pixels. Scalar output agreement can therefore occur by cancellation while nuisance/background estimation is not converged. Add explicit aperture and background-RMS stability requirements before response opening.',
  'criteria_v2':criteria,'steps':steps,'smallest_converged_size_pix_v2':chosen,
  'status':'freeze_candidate' if chosen is not None else 'no_converged_size_in_tested_ladder',
  'note':'This is a mechanical source-side criterion repair triggered by an internal nuisance-convergence failure. It does not use or inspect validation kinematics or Gate outcomes.'
}
(outdir/'convergence_decision_v2.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
print(json.dumps(report,indent=2,sort_keys=True))
