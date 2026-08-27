# Stage-0A AutoProf checkfit interpretation

Status: source-only implementation interpretation recorded before any validation kinematic catalogue/product is opened and before the production optical pass.

This note does **not** alter the Stage-0A `full_stellar` gate. It records how AutoProf v1.3.4 `checkfit` diagnostics are to be carried.

## Decision

AutoProf `checkfit` booleans are recorded as diagnostic metadata and are **not** additional hard sample-exclusion cuts.

The Stage-0A hard optical completeness gate remains exactly the already frozen rule: all three bands; at least five valid r-band profile points before the common 0.22 mag arcsec^-2 uncertainty threshold; finite r-band half-light radius; and all 45 MLCR stellar-mass estimates finite.

A pipeline execution failure, malformed/non-finite required profile, failed common-aperture construction, non-finite R50, or non-finite 45-way MLCR result still fails closed through that existing gate.

## Rationale

AutoProf v1.3.4 documents `Check_Fit` as a detector of potentially failed isophote solutions, but explicitly states that failure of one or more checks can instead reflect strong non-axisymmetric galaxy structure while the fit itself remains acceptable. In particular its `Light symmetry` test uses the first Fourier coefficient to detect either a failed centre or a lopsided/disturbed galaxy.

Making `Light symmetry == pass` a new hard cut after Stage-0A would therefore do two undesirable things:

1. change the prospectively frozen `full_stellar` definition; and
2. preferentially remove lopsided/disturbed optical systems, creating a source-side morphology selection that could be physically entangled with the structural/history variables under test.

Therefore all four AutoProf check flags (`isophote variability`, `FFT coefficients`, `initial fit compare`, `Light symmetry`) are retained in the source-only metadata and later reported as diagnostics/sensitivity strata, not used to select the primary sample.

## Pilot observation

The first source-only r-band pilot (`WALLABY J124947+035042`) completed under AutoProf v1.3.4 with FFT coefficients, initial-fit comparison and isophote-variability checks passing, while Light symmetry failed. The profile itself remains measured with small surface-brightness uncertainty through approximately 59.94 arcsec and first crosses the frozen 0.22-mag uncertainty threshold at approximately 65.93 arcsec. No validation kinematic quantity was consulted in making this interpretation.

## Firewall

This note contains no validation kinematic values and does not authorize opening any validation kinematic catalogue/product. Stage-1 remains unearned until the complete source-only predictor/optical tables and required receipts are sealed externally.
