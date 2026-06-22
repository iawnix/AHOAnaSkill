#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem import rdFingerprintGenerator

from aar_db import DB
from aar_mol import MolReader, StructQuery
from constants import DEFAULT_DB

FAMILIES = ("rdkit_desc", "rdkit_count", "rdkit_fg", "rdkit_morgan")


class RDKitFeaturizer:
    def __init__(self, db_path: str | Path = DEFAULT_DB):
        self.db = DB(db_path)

    def close(self) -> None:
        self.db.close()

    def compute_for_molecule(self, molecule_id: int, families: list[str] | None = None) -> int:
        row = self.db.conn.execute("SELECT * FROM molecules WHERE id = ?", (molecule_id,)).fetchone()
        if row is None:
            raise ValueError(f"Unknown molecule_id: {molecule_id}")
        mol = MolReader.load_sdf(row["sdf_path"])
        mol_no_h = Chem.RemoveHs(mol, sanitize=False)
        smiles = row["smiles"] or Chem.MolToSmiles(mol_no_h)
        selected = families or list(FAMILIES)
        features: dict[str, dict[str, float]] = {}
        if "rdkit_desc" in selected:
            features["rdkit_desc"] = self._standard_desc(mol_no_h)
        if "rdkit_count" in selected:
            features["rdkit_count"] = self._counts(mol_no_h, smiles)
        if "rdkit_fg" in selected:
            features["rdkit_fg"] = {k: float(v) for k, v in StructQuery.functional_groups(smiles).items()}
        if "rdkit_morgan" in selected:
            features["rdkit_morgan"] = self._morgan_fp(mol_no_h)
        inserted = 0
        for family, family_features in features.items():
            for name, value in family_features.items():
                if value is None or not math.isfinite(float(value)):
                    continue
                self.db.insert_feature(molecule_id, family, name, float(value))
                inserted += 1
        self.db.conn.commit()
        return inserted

    def compute_all_new(self, family: str = "all", role: str | None = None) -> dict:
        if family == "all":
            families = list(FAMILIES)
        elif family in FAMILIES:
            families = [family]
        else:
            raise ValueError(f"Unknown feature family: {family}")
        sql = "SELECT id, name, role FROM molecules"
        params: list[str] = []
        if role:
            sql += " WHERE role = ?"
            params.append(role)
        rows = self.db.conn.execute(sql, params).fetchall()
        total = 0
        for row in rows:
            total += self.compute_for_molecule(int(row["id"]), families=families)
        return {
            "family": family,
            "role": role,
            "n_molecules": len(rows),
            "n_feature_values_upserted": total,
        }

    def stats(self) -> dict[str, int]:
        return {
            row["family"]: int(row["n"])
            for row in self.db.conn.execute("SELECT family, COUNT(*) AS n FROM features GROUP BY family ORDER BY family")
        }

    @staticmethod
    def _standard_desc(mol: Chem.Mol) -> dict[str, float]:
        values = {}
        for name, fn in Descriptors._descList:
            try:
                value = fn(mol)
            except Exception:
                continue
            if value is not None and math.isfinite(float(value)):
                values[name] = float(value)
        return values

    @staticmethod
    def _counts(mol: Chem.Mol, smiles: str) -> dict[str, float]:
        centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
        return {
            "n_atoms": float(mol.GetNumAtoms()),
            "n_heavy_atoms": float(mol.GetNumHeavyAtoms()),
            "n_rings": float(rdMolDescriptors.CalcNumRings(mol)),
            "n_aromatic_rings": float(rdMolDescriptors.CalcNumAromaticRings(mol)),
            "n_aliphatic_rings": float(rdMolDescriptors.CalcNumAliphaticRings(mol)),
            "n_stereo_centers": float(len(centers)),
            "fraction_csp3": float(Lipinski.FractionCSP3(mol)),
            "logp": float(Crippen.MolLogP(mol)),
            "tpsa": float(rdMolDescriptors.CalcTPSA(mol)),
            "n_ring_systems": float(StructQuery.count_ring_systems(smiles)),
            "n_aromatic_atoms": float(StructQuery.count_aromatic_atoms(smiles)),
            "longest_aliphatic_chain": float(StructQuery.longest_aliphatic_chain(smiles)),
        }

    @staticmethod
    def _morgan_fp(mol: Chem.Mol, radius: int = 2, n_bits: int = 1024) -> dict[str, float]:
        generator = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
        fp = generator.GetFingerprint(mol)
        return {f"bit_{bit}": 1.0 for bit in fp.GetOnBits()}


def _json_print(obj: object) -> None:
    print(json.dumps(obj, ensure_ascii=False, indent=2, sort_keys=True))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="AHO RDKit baseline featurizer")
    parser.add_argument("--db", dest="root_db", default=str(DEFAULT_DB))
    sub = parser.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("compute")
    p.add_argument("--db")
    p.add_argument("--family", default="all", choices=["all", *FAMILIES])
    p.add_argument("--role", choices=["REA", "CAT", "SOL", "PRO_R", "PRO_S"])
    p.add_argument("--json", action="store_true")

    p = sub.add_parser("stats")
    p.add_argument("--db")
    p.add_argument("--json", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    db_path = args.db or args.root_db
    featurizer = RDKitFeaturizer(db_path)
    try:
        if args.cmd == "compute":
            result = featurizer.compute_all_new(family=args.family, role=args.role)
            _json_print(result) if args.json else print(result)
        elif args.cmd == "stats":
            result = featurizer.stats()
            _json_print(result) if args.json else print(result)
    finally:
        featurizer.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
