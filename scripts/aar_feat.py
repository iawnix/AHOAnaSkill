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
        mol, smiles = self._analysis_mol(row)
        selected = families or list(FAMILIES)
        features: dict[str, dict[str, float]] = {}
        if "rdkit_desc" in selected:
            features["rdkit_desc"] = self._standard_desc(mol)
        if "rdkit_count" in selected:
            features["rdkit_count"] = self._counts(mol, smiles)
        if "rdkit_fg" in selected:
            features["rdkit_fg"] = self._functional_groups(smiles)
        if "rdkit_morgan" in selected:
            features["rdkit_morgan"] = self._morgan_fp(mol)
        inserted = 0
        for family, family_features in features.items():
            for name, value in family_features.items():
                if value is None or not math.isfinite(float(value)):
                    continue
                self.db.insert_feature(molecule_id, family, name, float(value))
                inserted += 1
        self.db.conn.commit()
        return inserted

    @staticmethod
    def _analysis_mol(row) -> tuple[Chem.Mol, str]:
        smiles = (row["smiles"] or "").strip()
        if smiles:
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                return mol, smiles

        mol = MolReader.load_sdf(row["sdf_path"])
        mol_no_h = Chem.RemoveHs(mol, sanitize=False)
        try:
            Chem.SanitizeMol(mol_no_h)
        except Exception:
            pass
        if not smiles:
            try:
                smiles = Chem.MolToSmiles(mol_no_h)
            except Exception:
                smiles = ""
        return mol_no_h, smiles

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
        values = {
            "n_atoms": float(mol.GetNumAtoms()),
            "n_heavy_atoms": float(mol.GetNumHeavyAtoms()),
        }
        for name, fn in {
            "n_rings": lambda: rdMolDescriptors.CalcNumRings(mol),
            "n_aromatic_rings": lambda: rdMolDescriptors.CalcNumAromaticRings(mol),
            "n_aliphatic_rings": lambda: rdMolDescriptors.CalcNumAliphaticRings(mol),
            "fraction_csp3": lambda: Lipinski.FractionCSP3(mol),
            "logp": lambda: Crippen.MolLogP(mol),
            "tpsa": lambda: rdMolDescriptors.CalcTPSA(mol),
        }.items():
            value = RDKitFeaturizer._safe_float(fn)
            if value is not None:
                values[name] = value

        centers = RDKitFeaturizer._safe_chiral_centers(mol)
        if centers is not None:
            values["n_stereo_centers"] = float(len(centers))

        if smiles:
            for name, fn in {
                "n_ring_systems": lambda: StructQuery.count_ring_systems(smiles),
                "n_aromatic_atoms": lambda: StructQuery.count_aromatic_atoms(smiles),
                "longest_aliphatic_chain": lambda: StructQuery.longest_aliphatic_chain(smiles),
            }.items():
                value = RDKitFeaturizer._safe_float(fn)
                if value is not None:
                    values[name] = value
        return values

    @staticmethod
    def _safe_chiral_centers(mol: Chem.Mol):
        try:
            return Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
        except Exception:
            try:
                return Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=True)
            except Exception:
                return None

    @staticmethod
    def _functional_groups(smiles: str) -> dict[str, float]:
        if not smiles:
            return {}
        try:
            return {k: float(v) for k, v in StructQuery.functional_groups(smiles).items()}
        except Exception:
            return {}

    @staticmethod
    def _safe_float(fn) -> float | None:
        try:
            value = fn()
        except Exception:
            return None
        if value is None or not math.isfinite(float(value)):
            return None
        return float(value)

    @staticmethod
    def _morgan_fp(mol: Chem.Mol, radius: int = 2, n_bits: int = 1024) -> dict[str, float]:
        try:
            generator = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
            fp = generator.GetFingerprint(mol)
        except Exception:
            return {}
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
