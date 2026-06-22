from __future__ import annotations

import csv
import sys
from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

from aar_db import DB
from aar_feat import RDKitFeaturizer
from aar_mol import MolReader, Renderer, StructQuery


def write_sdf(path: Path, smiles: str) -> None:
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=42) == 0
    AllChem.UFFOptimizeMolecule(mol, maxIters=100)
    mol.SetProp("SMILES", smiles)
    writer = Chem.SDWriter(str(path))
    writer.write(mol)
    writer.close()


@pytest.fixture()
def sample_dataset(tmp_path: Path):
    sdf_dir = tmp_path / "sdf"
    sdf_dir.mkdir()
    molecules = {
        "CAT-1": "P(c1ccccc1)(c1ccccc1)c1ccccc1",
        "CAT-2": "CC(C)(C)P(c1ccccc1)c1ccccc1",
        "SOL-1": "CCO",
        "REA-1": "COc1ccccc1C=C",
        "PRO-R-1": "COc1ccccc1[C@H](C)O",
        "PRO-S-1": "COc1ccccc1[C@@H](C)O",
    }
    for name, smiles in molecules.items():
        write_sdf(sdf_dir / f"{name}.sdf", smiles)

    csv_path = tmp_path / "aho.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "DATA_ID",
                "CAT_NAME",
                "SOL_NAME",
                "PRO_R_NAME",
                "PRO_S_NAME",
                "REA_NAME",
                "TEMP",
                "PRESSURE",
                "EE",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "DATA_ID": "D1",
                "CAT_NAME": "CAT-1",
                "SOL_NAME": "SOL-1",
                "PRO_R_NAME": "PRO-R-1",
                "PRO_S_NAME": "PRO-S-1",
                "REA_NAME": "REA-1",
                "TEMP": "298.15",
                "PRESSURE": "20",
                "EE": "0.75",
            }
        )
        writer.writerow(
            {
                "DATA_ID": "D2",
                "CAT_NAME": "CAT-2",
                "SOL_NAME": "SOL-1",
                "PRO_R_NAME": "PRO-R-1",
                "PRO_S_NAME": "",
                "REA_NAME": "REA-1",
                "TEMP": "308.15",
                "PRESSURE": "30",
                "EE": "-0.20",
            }
        )
    return csv_path, sdf_dir


def test_import_features_query_and_render(sample_dataset, tmp_path: Path):
    csv_path, sdf_dir = sample_dataset
    db_path = tmp_path / "aho.sqlite"
    summary = DB.import_csv(csv_path, sdf_dir, db_path)
    assert summary["n_rows"] == 2
    assert Path(summary["sdf_index"]) == tmp_path / "sdf_index.json"
    assert (tmp_path / "sdf_index.json").exists()

    with DB(db_path) as db:
        stats = db.stats()
        assert stats["n_reactions"] == 2
        assert stats["molecules_by_role"]["PRO_S"] == 1

    featurizer = RDKitFeaturizer(db_path)
    try:
        feature_summary = featurizer.compute_all_new(family="all")
        assert feature_summary["n_molecules"] == 6
    finally:
        featurizer.close()

    with DB(db_path) as db:
        rows = db.query_with_features(families=["rdkit_count"])
        assert len(rows) == 2
        assert any(key.endswith("__n_atoms") for key in rows[0])

    assert StructQuery.smarts_count("COc1ccccc1", "c-OMe") == 1
    assert MolReader.read_sdf(sdf_dir / "REA-1.sdf")["smiles"] == "COc1ccccc1C=C"

    out = tmp_path / "rea.svg"
    Renderer.render(sdf_dir / "REA-1.sdf", out)
    assert out.exists()
    assert "<svg" in out.read_text(encoding="utf-8")[:300]


def test_missing_sdf_is_blocking(sample_dataset, tmp_path: Path):
    csv_path, sdf_dir = sample_dataset
    (sdf_dir / "CAT-2.sdf").unlink()
    with pytest.raises(FileNotFoundError, match="Missing required SDF"):
        DB.import_csv(csv_path, sdf_dir, tmp_path / "aho.sqlite")
