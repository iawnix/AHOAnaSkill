---
name: ahoa-na-skill
description: Analyze asymmetric hydrogenation experimental datasets for SAR and AAR patterns. Use when working with AHO/AAR/asymmetric hydrogenation CSV plus SDF data, enantioselectivity trends, catalyst/reactant/solvent effects, RDKit descriptors, structure-activity relationship reports, or ad-hoc structural hypotheses that need traceable scripts and figures.
---

# AHOAnaSkill

## Overview

Use this skill to ingest asymmetric hydrogenation data, compute a stable RDKit baseline, inspect structures, test agent-authored structural hypotheses, and write traceable SAR reports. Keep the Python modules responsible for deterministic data plumbing and objective computation; keep interpretation, hypothesis generation, and final conclusions agent-led.

## Required Data

Use a CSV with exactly these columns in this order:

```text
DATA_ID,CAT_NAME,SOL_NAME,PRO_R_NAME,PRO_S_NAME,REA_NAME,TEMP,PRESSURE,EE
```

Expect `TEMP` in K, `PRESSURE` in bar, and signed decimal `EE` in `[-1, 1]` where positive means R-major and negative means S-major. `PRO_S_NAME` may be empty.

Use a flat SDF directory. Each CSV molecule name must match `<NAME>.sdf` exactly, case-sensitive. Each SDF should contain a `SMILES` property; the tools fall back to `MolToSmiles` when it is absent.

## Command Rule

Call the project through `aho`:

```bash
aho check-env
aho ingest --csv data.csv --sdf-dir sdf/
aho feat compute --family all
aho db query-with-features --json
```

Do not call `python scripts/*.py` in normal use, because that can bypass the conda environment. For local development before the conda env exists, `AHO_NO_CONDA=1 bin/aho ...` is acceptable.

Run `./install.sh` first in the skill root when setting up a fresh copy. The installer creates or updates the `AHOAnaSkill` Conda environment and writes `.aho-runtime.env`; `bin/aho` reads that file and uses the installed environment by default.

## Workflow

1. Check the environment:

   ```bash
   aho check-env
   aho stats
   ```

2. Ingest the dataset:

   ```bash
   aho ingest --csv X.csv --sdf-dir SDF_DIR
   ```

   Treat missing SDF files as blocking data-quality issues. Do not silently skip rows.

3. Compute or refresh baseline features:

   ```bash
   aho feat compute --family all
   ```

4. Export a wide analysis table:

   ```bash
   aho db query-with-features --json > reports/working_table.json
   ```

5. Inspect structures before making structural claims. Render at least one representative molecule for every role present and at least one `PRO_R`/`PRO_S` pair when both are available:

   ```bash
   aho render --sdf sdf/CAT-53.sdf --out reports/figures/cat-53.svg
   aho overlay --sdf-a sdf/PRO-R.sdf --sdf-b sdf/PRO-S.sdf --out reports/figures/product_pair.svg
   ```

6. Form hypotheses from the data:

   - Same `REA` and same `CAT` with changed `TEMP`, `PRESSURE`, or `SOL`: condition effects.
   - Same `REA` with changed `CAT`: catalyst screening.
   - Same `CAT` with changed `REA`: substrate scope.
   - EE sign reversal: analyze as a dedicated case.
   - Steric, ortho, coordination, or conformational claims: confirm with structure queries and figures before reporting.

7. Quantify each ad-hoc structural hypothesis with a script archived under `reports/scripts/<timestamp>_<topic>.py`. Use `StructQuery` and `DB`; then run:

   ```bash
   aho run-custom --script reports/scripts/<timestamp>_<topic>.py
   ```

   Persist new features under `features.family='custom_<topic>'` and register `custom_features_meta` with a verdict of `supported`, `rejected`, or `inconclusive`.

8. Use basic statistics: correlations, grouped means, Mann-Whitney U where appropriate, and explicit outlier checks. Keep small-sample limitations visible.

9. Write the report to `reports/<timestamp>_<topic>.md`. For every conclusion, include `n`, evidence, counterexamples, and figure paths. Include a `Reasoning Trace` table listing every custom feature family, hypothesis, verdict, molecule count, and script path.

## Resource Map

- `bin/aho`: only supported command entry point.
- `scripts/aar_db.py`: SQLite schema, import, query, custom script execution.
- `scripts/aar_feat.py`: RDKit descriptor, count, functional group, and Morgan fingerprint features.
- `scripts/aar_mol.py`: SDF reading, structure queries, and rendering.
- `scripts/ingest.py`: end-to-end CSV/SDF ingest plus optional feature computation.
- `docs/AHOAnaSkill_PLAN.md`: original implementation plan and design rationale.
