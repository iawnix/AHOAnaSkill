# AHOAnaSkill

AHOAnaSkill is a standalone skill/project for asymmetric hydrogenation (AHO) SAR analysis. It ingests fixed-schema CSV plus flat SDF directories, stores reactions and molecules in SQLite, computes RDKit baseline features, provides structure-query primitives, and supports traceable ad-hoc custom features.

The project is not installed into `~/.codex/skills` or `~/.claude/skills` by default. The repository root is the skill root and can be copied or linked into a runtime skill directory later.

## Install

```bash
./install.sh
export PATH="$PWD/bin:$PATH"
```

`install.sh` creates or updates the Conda environment named `AHOAnaSkill`,
then writes `.aho-runtime.env` in this skill root. `bin/aho` reads that file and
uses the recorded Conda executable and environment by default, so normal use
should call `aho ...` without `AHO_NO_CONDA=1`.

For development without the conda environment:

```bash
AHO_NO_CONDA=1 bin/aho check-env --allow-missing-renderer
```

## Data Contract

CSV columns must be exactly:

```text
DATA_ID,CAT_NAME,SOL_NAME,PRO_R_NAME,PRO_S_NAME,REA_NAME,TEMP,PRESSURE,EE
```

SDF files must live in one flat directory and use exact molecule names:

```text
CAT-53.sdf
SOL-01.sdf
REA-01.sdf
PRO-R-01.sdf
```

Each SDF should contain a `SMILES` property. The reader falls back to RDKit canonical SMILES if the property is absent.

### aho commands

```bash
aho check-env
aho ingest --csv X.csv --sdf-dir sdf/
aho stats
aho db query --catalyst CAT-53 --json
aho db query-with-features --json
aho feat compute --family all
aho feat compute --family rdkit_morgan --role CAT
aho mol smarts-count --smi "COc1ccccc1" --pattern c-OMe
aho mol atom-report --smi "C[C@H](O)c1ccccc1" --json
aho render --sdf sdf/CAT-53.sdf --out reports/figures/cat-53.svg
aho overlay --sdf-a sdf/PRO-R.sdf --sdf-b sdf/PRO-S.sdf --out reports/figures/product_pair.svg
aho smarts --sdf sdf/REA-01.sdf --pattern c-OMe --out reports/figures/rea_ome.svg
aho run-custom --script reports/scripts/20260622_ortho_ome.py
```

### Implementation Notes

- `scripts/aar_db.py`: schema, CSV import, reaction queries, custom feature execution.
- `scripts/aar_feat.py`: RDKit descriptors, counts, functional groups, Morgan bits.
- `scripts/aar_mol.py`: SDF reader, structure primitives, RDKit fallback renderer.
- `scripts/ingest.py`: end-to-end import plus optional baseline feature computation.
- `docs/AHOAnaSkill_PLAN.md`: original Chinese implementation plan.

The rendering layer uses RDKit drawing as a reliable fallback. `xyzrender` is declared in `environment.yml`; richer rendering can be added behind the same `Renderer` API without changing the `aho` command surface.

## Custom Feature Script Contract

Archive custom scripts under `reports/scripts/`. A script may expose:

```python
def run(db_path: str, skill_dir: str) -> dict:
    return {
        "features": [
            {
                "molecule_name": "REA-01",
                "role": "REA",
                "family": "custom_ortho_ome",
                "feature_name": "n_ortho_ome",
                "feature_value": 1,
            }
        ],
        "meta": {
            "family_topic": "custom_ortho_ome",
            "description": "Counts aromatic methoxy motif",
            "hypothesis": "Ortho OMe correlates with lower |ee|",
            "verdict": "inconclusive",
            "n_molecules": 1,
        },
    }
```

Run it with:

```bash
aho run-custom --script reports/scripts/20260622_ortho_ome.py
```

The runner sets `AHO_DB_PATH` and `AHO_SKILL_DIR`, adds `scripts/` to `sys.path`, archives external scripts into `reports/scripts/`, inserts returned features, and updates `custom_features_meta`.

## Validation

```bash
python3 -m pytest -q tests
python3 /home/iaw/.codex/skills/.system/skill-creator/scripts/quick_validate.py .
```
