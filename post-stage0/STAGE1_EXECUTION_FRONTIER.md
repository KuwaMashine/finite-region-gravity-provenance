# WALLABY Stage-1 execution frontier

Status: **source-only execution in progress; validation kinematics remain closed**.

This note is operational provenance only. It does not change any frozen scientific rule, source-quality gate, optical rule, statistic, threshold, or holdout adapter.

## Paid source-side evidence

- Public Stage-0 timestamp proof: commit `15abeac9b6d285ce42caf26498115c0686d9bcd3`; bound pretarget manifest file SHA-256 `61fa4a6f69d93df4616cc0e0ad664797e8cf89c0cd5d6f65a31c0dd695c82a0c`.
- External Stage-0 receipt: `post-stage0/STAGE0_EXTERNAL_TIMESTAMP_RECEIPT.json`.
- Exact v4.99 source-engine copy verification: `post-stage0/frozen-v4.99-source-engine/VERIFICATION.json`, status `pass`.
- Full H I source-only freeze: commit `229a4780c7ce7d1f3dc327c4392473a7d0e776a5`.
  - final release rows: 1760
  - finite qflag=0 source candidates: 1525
  - H I source-profile pass: 747
  - H I source-profile reject: 778
  - pass by field: NGC4808 110; NGC5044 539; Vela 98
- Vela DR10 source-side coverage audit: commit `66522f8639ffa530fe05b07080a128d7d903cc75`.
  - complete padded final Vela TR1 rectangle has zero DR10 Tractor objects
  - no substitute survey is permitted under the frozen all-g/r/z DR10 gate
  - therefore the 98 Vela H I survivors are source-side `full_stellar=false`
- Optical spatial preplan: commit `200f353f6122470d082a59d05142280505e0f99c`.
  - H I survivors: 747
  - NGC4808/NGC5044 optical candidates: 649
  - Vela no-DR10 classifications: 98
  - spatial shards: 24
- Production optical frame lock: `post-stage0/SOURCE_OPTICAL_PRODUCTION_FRAME_LOCK.json`, SHA-256 `0a0d028d771d799cb71c933acbd4604befb728a6c187ae7d39a8444845d15c56`.
- Exact production optical preflight pass: commit `c7b19380bb39fc53c555643e724f31b5b87225eb`; receipt SHA-256 `f75a9275fa4e32055f3ad09701d4587ff03e2d24427e5cad46e1a9245c90c350`.
  - real H I survivor `WALLABY J124947+035042`
  - `full_stellar=true`
  - selected registered frame: 2560 pixels
  - response products used: false

## Live optical production

The GitHub recursion guard prevented the bot-authored preflight receipt from triggering a second push workflow. The workflow was changed only to accept an explicit orchestration trigger path; scientific job content was not changed.

- orchestration-only workflow change: commit `bbafabf4b2e6261ca9dab172fecbbf5fe1da6fbc`
- explicit trigger receipt: commit `8cf5fda5d79939d1d42379c310158429729bd8fc`
- live GitHub Actions run: **`33139915501`**, `WALLABY Stage-1 full source-only optical pass`
- 24 spatial shards; maximum six parallel workers
- pinned `autoprof==1.3.4`; exact registered 1024--3072 frame ladder
- exact primary-brick pixels, neighbour-only bilinear reprojection, per-shard brick cache
- consolidation refuses to freeze with unresolved network or implementation failures

Expected successful consolidation location:
`post-stage0/source-optical-stage1/`
with `meta.csv`, `optical_profile.csv`, `optical_status.csv`, `frame_audit.csv`, `SOURCE_OPTICAL_STAGE1_REPORT.json`, and hashes.

## Armed next step

Workflow `source_stage1_target_full.yml` was added at commit `0da31db9d97071f303b91b68db285e69fea06ec0` but has **not** been triggered.

After the optical population is successfully committed:

1. create an explicit source-stage trigger receipt that binds that optical freeze commit;
2. build `post-stage0/SOURCE_EXPORT_MANIFEST.json` from the committed normalized source files;
3. run byte-verified v4.99 `run_source_stage.py --mode TARGET` using `STAGE0_EXTERNAL_TIMESTAMP_RECEIPT.json`;
4. commit the resulting `source_predictors.csv`, `baryonic_hi_fold.csv`, `proxy_strata.csv`, `proxy_strata_metadata.json`, `source_beam_report.json`, `source_prepare_report.json`, `STAGE1_TIMESTAMP_PAYLOAD.json`, and hashes;
5. use that public source-stage commit as the external Stage-1 timestamp proof and create a receipt binding the frozen manifest, proxy-strata map and source beam report;
6. only after that receipt verifies may validation kinematics or Gate-A/B outcomes be opened.

## Roadmap state

Do **not** tick the Phase-7 target/timestamp gate yet. Stage-0 is paid, but Stage-1 source-map/strata timestamping remains outstanding. No new roadmap checkbox has been claimed by this execution block, so the three-item ZIP checkpoint counter has not advanced.

## Firewall

No validation kinematic catalogue/product, `kflag`, `kin_pa`, `Vrot_model`, `QFlag_model`, `gobs`, `gbar`, `delta_ref`, `rar_residual`, `Vmax`, `Vflat`, or Gate-A/B outcome has been opened in this Stage-1 source-side work.
