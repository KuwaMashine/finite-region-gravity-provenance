# WALLABY Gate-A/B Stage-0A source/optical implementation lock

Status: prospective supplemental lock, created after the original Stage-0 target seal and before any validation optical pixel values or any validation kinematic catalogue values are opened.

This document does **not** rewrite the Stage-0 hypothesis or target-side constants. It resolves source/optical adapter details that the published WALLABY description and the original Stage-0 payload did not specify sufficiently for byte-for-byte reproduction.

## 1. Firewall

Until the Stage-1 source-only receipt is externally timestamped, no validation kinematic catalogue/product may be queried or read. In particular the source-side process must not read `kflag`, `kin_pa`, `Vrot_model`, `QFlag_model`, `gobs`, `gbar`, `delta_ref`, `rar_residual`, `Vmax`, or `Vflat`.

The final Phase-2 source population is defined by the final source releases only:

- NGC 4808 TR1
- NGC 5044 TR3
- Vela TR1

NGC 5044 TR1 and TR2 are superseded and excluded.

## 2. Source quality

`robust_sample` is reconstructed without kinematic quality information.

Primary source-quality entry requires:

1. final Phase-2 source release above;
2. source catalogue `qflag == 0`;
3. finite source-only distance, H I mass and source ellipse metadata;
4. successful source moment-0/mask extraction;
5. at least eight H I radial samples after deterministic source-profile construction;
6. signed source moment-0 pixels are retained; no pixel-level clipping is permitted;
7. a galaxy fails source-profile construction if any integrated radial annular source weight used by the frozen predictor is negative or non-finite.

`full_stellar` is a separate optical completeness flag and is not inferred from the historical kinematic robust sample.

The primary resolution lane remains `ell_maj_beams >= 2.0`; `1.3 <= ell_maj_beams < 2.0` is diagnostic only.

Machine-readable contract: `STAGE0A_SOURCE_QUALITY_LOCK.json`.
SHA-256 of its exact UTF-8 contents: `f86549a599dad5e7951c8d5ec066a0a29e76ee7f22d3e2520d4caa0bcc7517f6`.

## 3. H I source-map geometry

For released WALLABY source products, the source moment-0 map is 2-D and the source mask is a 3-D spectral cube. The 2-D source footprint is frozen as

`mask2d(x,y) = any(mask3d(nu,x,y) > 0, spectral axis)`.

The released moment-0 map is used with its signed pixel values inside this projected footprint. No response or kinematic map is used.

The source-catalogue `ell_pa` is treated as the sky-plane major-axis position angle, north through east, modulo 180 degrees. This convention was verified prospectively against the WCS-transformed principal axis of a source-only Vela moment-0 map; no kinematic PA is substituted.

## 4. Optical image/profile extraction

Optical imaging is DESI Legacy Imaging Surveys DR10 `g,r,z` imaging. AutoProf is pinned to public release **v1.3.4**.

The r band is the geometric fit band. The g and z profiles are forced photometry using the r-band isophotal solution, matching the published WALLABY procedure.

The source position is the H I source-catalogue sky centroid and is the only target-side centering input. No kinematic center, PA, inclination or model is permitted.

The photometric zero point is 22.5 for nanomaggy images. Galactic foreground extinction is corrected multiplicatively in flux using the local Legacy Survey/Schlafly-Finkbeiner transmission for each band. No internal/inclination attenuation correction is applied.

No k-correction is applied in the primary Stage-1 optical lane. This is an explicit prospective operational convention because the WALLABY paper specifies the 0.22-mag profile truncation and MLCR construction but does not specify a k-correction in that procedure. A k-corrected result, if ever examined, must be registered separately as a source-only sensitivity and may not replace the primary after response opening.

For a common colour aperture, retain profile radii strictly interior to the first forced-profile radius at which **any** of `muerr_g`, `muerr_r`, or `muerr_z` reaches or exceeds `0.22 mag arcsec^-2`. This removes an ambiguity in the published phrase while ensuring all colours refer to one aperture.

`full_stellar` requires all three bands, at least five valid r-band profile points before the threshold, a finite r-band half-light radius, and all 45 MLCR mass estimates finite.

## 5. R50 and stellar scale length

Within the frozen optical aperture, integrate the r-band profile cumulatively using the local AutoProf ellipticity for annular area. `R50,r` is the radius enclosing 50% of that truncated r-band light, with linear interpolation between bracketing cumulative samples.

The stellar exponential-equivalent scale is fixed as

`R_d,* = R50,r / 1.6783469900166605`.

No galaxy-specific exponential fitting interval is introduced.

## 6. Stellar mass architecture

The WALLABY publication explicitly specifies five MLCR families, three colours (`g-r`, `g-z`, `r-z`) and three luminosity bands (`g`,`r`,`z`), yielding **45** stellar-mass estimates. The adopted stellar mass is their median and the quoted method uncertainty is their standard deviation.

The publication does not identify all hidden table/model-row choices needed to reconstruct those five families exactly. Therefore the following choices are frozen prospectively before optical target values are opened; they are not represented as recovered hidden WALLABY source code:

1. Roediger & Courteau 2015 BC03;
2. Roediger & Courteau 2015 FSPS;
3. Zhang, Puzia & Weisz 2017 FSPS — **40 Local-Group dwarf galaxies** table;
4. Zhang, Puzia & Weisz 2017 BC03 — **40 Local-Group dwarf galaxies** table;
5. García-Benito et al. 2019 — CBe base, Chabrier IMF, HUB3 `Late` (Sc+Sd), `SynR` reddened synthetic spectra.

Every relation is normalized to

`log10(M/L_band) = a + b * colour`.

The exact 45 coefficients are frozen in `STAGE0A_MLCR_COEFFICIENTS.csv`.
SHA-256 of its exact UTF-8 contents: `93aa2e7d1cdc4e619fd6942fffad199ff42218d376b23839ccd5bfb9dbb9cae5`.

The Zhang 40-LG choice is made because it is the directly reconstructed real Local-Group dwarf calibration rather than the synthetically expanded full sample. The García-Benito `Late/SynR` choice is made because the public PyCASSO database defines `Late` as Sc+Sd and recommends `SynR` for observational use.

## 7. Solar magnitudes and luminosities

For Legacy Survey southern/DECam-equivalent AB photometry, freeze the Willmer (2018) solar absolute magnitudes:

- `M_sun,g = 5.05`
- `M_sun,r = 4.61`
- `M_sun,z = 4.50`

For each band, after foreground correction,

`L_band/L_sun = 10^[-0.4 (M_band - M_sun,band)]`.

For each of the 45 coefficient rows,

`log10 Mstar = log10 L_band + a + b*colour`.

The source-side stellar mass is

`log10 Mstar_adopted = median(45 log10 Mstar estimates)`

and the method spread is the standard deviation of those 45 log masses.

## 8. No outcome-driven repair

Once this Stage-0A lock is committed, the coefficient families, quality rule, common colour aperture, extinction-only primary correction, AutoProf version, `R50 -> Rd` conversion and resolution lanes may not be changed in response to validation kinematics or Gate-A/B outcomes.

Any later source-only implementation bug must be documented as a mechanical correction with before/after hashes and may not use response values to choose between alternatives.

## 9. Stage-1 condition

Optical/H I source processing may begin only after this lock. Stage-1 is earned only when the completed source-only predictor table, baryonic H I fold table, optical/stellar metadata and beam/recovery report are hashed and externally timestamped together with the original Stage-0 receipt and this Stage-0A receipt.

Only after that external Stage-1 receipt may the validation kinematic catalogue be opened.

## Literature basis

- Deg et al., WALLABY Pilot Survey scaling-relations paper: r-band AutoProf geometry, forced g/z photometry, 0.22 mag arcsec^-2 uncertainty truncation, five MLCR families, 3 colours x 3 luminosity bands = 45 masses, median adopted and standard deviation uncertainty.
- Roediger & Courteau (2015), MNRAS 452, 3209: BC03 and FSPS MLCR coefficient tables.
- Zhang, Puzia & Weisz (2017), MNRAS 466, 3217: FSPS/BC03 Local-Group-dwarf MLCR tables.
- García-Benito et al. (2019), A&A 621, A120 and public PyCASSO MLCR database: morphology/mode-resolved MLCRs.
- Willmer (2018), ApJS 236, 47: solar absolute magnitudes.
- Legacy Survey DR10 documentation: AB nanomaggy fluxes and `MW_TRANSMISSION` foreground-extinction metadata.
- AutoProf public release v1.3.4.
