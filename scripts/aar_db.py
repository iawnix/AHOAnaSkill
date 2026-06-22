#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import os
import runpy
import shutil
import sqlite3
import sys
from pathlib import Path
from typing import Any

from aar_mol import MolReader
from constants import (
    CSV_COLUMNS,
    DEFAULT_CUSTOM_SCRIPTS_DIR,
    DEFAULT_DB,
    ROLE_TO_COLUMN,
    ROLES,
    SDF_INDEX_PATH,
    SKILL_DIR,
    utc_now,
)

SCHEMA = """
PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS molecules (
    id            INTEGER PRIMARY KEY AUTOINCREMENT,
    name          TEXT NOT NULL,
    role          TEXT NOT NULL CHECK(role IN ('REA','CAT','SOL','PRO_R','PRO_S')),
    smiles        TEXT,
    inchi_key     TEXT,
    sdf_path      TEXT,
    n_atoms       INTEGER,
    mw            REAL,
    created_at    TEXT,
    UNIQUE(name, role)
);

CREATE TABLE IF NOT EXISTS reactions (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    data_id         TEXT UNIQUE NOT NULL,
    reactant_id     INTEGER REFERENCES molecules(id),
    catalyst_id     INTEGER REFERENCES molecules(id),
    solvent_id      INTEGER REFERENCES molecules(id),
    product_r_id    INTEGER REFERENCES molecules(id),
    product_s_id    INTEGER REFERENCES molecules(id),
    temperature_k   REAL,
    pressure_bar    REAL,
    ee              REAL CHECK(ee >= -1 AND ee <= 1),
    source          TEXT,
    notes           TEXT,
    created_at      TEXT
);

CREATE TABLE IF NOT EXISTS features (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    molecule_id     INTEGER NOT NULL REFERENCES molecules(id),
    family          TEXT NOT NULL,
    feature_name    TEXT NOT NULL,
    feature_value   REAL,
    computed_at     TEXT,
    UNIQUE(molecule_id, family, feature_name)
);

CREATE TABLE IF NOT EXISTS custom_features_meta (
    family_topic    TEXT PRIMARY KEY,
    description     TEXT,
    hypothesis      TEXT,
    script_path     TEXT,
    verdict         TEXT CHECK(verdict IN ('supported','rejected','inconclusive') OR verdict IS NULL),
    n_molecules     INTEGER,
    created_at      TEXT
);

CREATE INDEX IF NOT EXISTS idx_reactions_data_id ON reactions(data_id);
CREATE INDEX IF NOT EXISTS idx_features_lookup ON features(molecule_id, family);
"""


class DB:
    def __init__(self, db_path: str | Path = DEFAULT_DB):
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self.conn = sqlite3.connect(self.db_path)
        self.conn.row_factory = sqlite3.Row
        self.conn.execute("PRAGMA foreign_keys = ON")
        self.init_schema()

    def close(self) -> None:
        self.conn.close()

    def __enter__(self) -> "DB":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()

    @classmethod
    def init(cls, db_path: str | Path = DEFAULT_DB) -> Path:
        with cls(db_path) as db:
            db.init_schema()
        return Path(db_path)

    def init_schema(self) -> None:
        self.conn.executescript(SCHEMA)
        self.conn.commit()

    def upsert_molecule(self, name: str, role: str, sdf_path: str | Path) -> int:
        if role not in ROLES:
            raise ValueError(f"Invalid role: {role}")
        info = MolReader.read_sdf(sdf_path)
        now = utc_now()
        self.conn.execute(
            """
            INSERT INTO molecules
                (name, role, smiles, inchi_key, sdf_path, n_atoms, mw, created_at)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(name, role) DO UPDATE SET
                smiles=excluded.smiles,
                inchi_key=excluded.inchi_key,
                sdf_path=excluded.sdf_path,
                n_atoms=excluded.n_atoms,
                mw=excluded.mw
            """,
            (
                name,
                role,
                info["smiles"],
                info["inchi_key"],
                str(Path(sdf_path).resolve()),
                info["n_atoms"],
                info["mw"],
                now,
            ),
        )
        row = self.conn.execute(
            "SELECT id FROM molecules WHERE name = ? AND role = ?",
            (name, role),
        ).fetchone()
        if row is None:
            raise RuntimeError(f"Failed to upsert molecule: {name}/{role}")
        return int(row["id"])

    def molecule_id(self, name: str, role: str) -> int | None:
        row = self.conn.execute(
            "SELECT id FROM molecules WHERE name = ? AND role = ?",
            (name, role),
        ).fetchone()
        return int(row["id"]) if row else None

    def insert_feature(self, molecule_id: int, family: str, feature_name: str, feature_value: float) -> None:
        self.conn.execute(
            """
            INSERT INTO features (molecule_id, family, feature_name, feature_value, computed_at)
            VALUES (?, ?, ?, ?, ?)
            ON CONFLICT(molecule_id, family, feature_name) DO UPDATE SET
                feature_value=excluded.feature_value,
                computed_at=excluded.computed_at
            """,
            (molecule_id, family, feature_name, float(feature_value), utc_now()),
        )

    def upsert_reaction(self, row: dict[str, str], sdf_index: dict[str, str], source: str | None = None) -> int:
        mol_ids: dict[str, int | None] = {}
        for role, column in ROLE_TO_COLUMN.items():
            name = (row.get(column) or "").strip()
            if not name and role == "PRO_S":
                mol_ids[role] = None
                continue
            if not name:
                raise ValueError(f"Missing required molecule name for role {role} in row {row.get('DATA_ID')}")
            sdf_path = MolReader.validate(name, sdf_index)
            if sdf_path is None:
                raise FileNotFoundError(f"Missing SDF for {role} {name}")
            mol_ids[role] = self.upsert_molecule(name, role, sdf_path)
        ee = float(row["EE"])
        if ee < -1.0 or ee > 1.0:
            raise ValueError(f"EE out of range [-1, 1] for {row['DATA_ID']}: {ee}")
        self.conn.execute(
            """
            INSERT INTO reactions
                (data_id, reactant_id, catalyst_id, solvent_id, product_r_id, product_s_id,
                 temperature_k, pressure_bar, ee, source, created_at)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(data_id) DO UPDATE SET
                reactant_id=excluded.reactant_id,
                catalyst_id=excluded.catalyst_id,
                solvent_id=excluded.solvent_id,
                product_r_id=excluded.product_r_id,
                product_s_id=excluded.product_s_id,
                temperature_k=excluded.temperature_k,
                pressure_bar=excluded.pressure_bar,
                ee=excluded.ee,
                source=excluded.source
            """,
            (
                row["DATA_ID"],
                mol_ids["REA"],
                mol_ids["CAT"],
                mol_ids["SOL"],
                mol_ids["PRO_R"],
                mol_ids["PRO_S"],
                float(row["TEMP"]),
                float(row["PRESSURE"]),
                ee,
                source,
                utc_now(),
            ),
        )
        out = self.conn.execute("SELECT id FROM reactions WHERE data_id = ?", (row["DATA_ID"],)).fetchone()
        if out is None:
            raise RuntimeError(f"Failed to upsert reaction: {row['DATA_ID']}")
        return int(out["id"])

    @classmethod
    def import_csv(cls, csv_path: str | Path, sdf_dir: str | Path, db_path: str | Path = DEFAULT_DB) -> dict[str, Any]:
        csv_path = Path(csv_path)
        sdf_index = MolReader.build_index(sdf_dir)
        rows = _read_csv_rows(csv_path)
        missing = _find_missing_sdfs(rows, sdf_index)
        if missing:
            sample = "\n".join(f"- {role}: {name}" for role, name in missing[:30])
            raise FileNotFoundError(f"Missing required SDF files ({len(missing)}):\n{sample}")
        cls.init(db_path)
        with cls(db_path) as db:
            for row in rows:
                db.upsert_reaction(row, sdf_index, source=str(csv_path.resolve()))
            db.conn.commit()
        SDF_INDEX_PATH.parent.mkdir(parents=True, exist_ok=True)
        SDF_INDEX_PATH.write_text(json.dumps(sdf_index, indent=2, sort_keys=True), encoding="utf-8")
        return {
            "csv": str(csv_path.resolve()),
            "sdf_dir": str(Path(sdf_dir).resolve()),
            "db": str(Path(db_path).resolve()),
            "n_rows": len(rows),
            "n_sdf_indexed": len(sdf_index),
        }

    def query_reactions(self, filters: dict[str, str] | None = None) -> list[dict[str, Any]]:
        filters = filters or {}
        role_aliases = {
            "reactant": "rea",
            "catalyst": "cat",
            "solvent": "sol",
            "product_r": "pr",
            "product_s": "ps",
        }
        sql = """
        SELECT
            r.id, r.data_id, r.temperature_k, r.pressure_bar, r.ee, r.source,
            r.reactant_id, r.catalyst_id, r.solvent_id, r.product_r_id, r.product_s_id,
            rea.name AS rea_name, cat.name AS cat_name, sol.name AS sol_name,
            pr.name AS pro_r_name, ps.name AS pro_s_name,
            rea.smiles AS rea_smiles, cat.smiles AS cat_smiles, sol.smiles AS sol_smiles,
            pr.smiles AS pro_r_smiles, ps.smiles AS pro_s_smiles
        FROM reactions r
        LEFT JOIN molecules rea ON rea.id = r.reactant_id
        LEFT JOIN molecules cat ON cat.id = r.catalyst_id
        LEFT JOIN molecules sol ON sol.id = r.solvent_id
        LEFT JOIN molecules pr ON pr.id = r.product_r_id
        LEFT JOIN molecules ps ON ps.id = r.product_s_id
        """
        clauses = []
        params: list[str] = []
        if filters.get("data_id"):
            clauses.append("r.data_id = ?")
            params.append(filters["data_id"])
        for key, alias in role_aliases.items():
            if filters.get(key):
                clauses.append(f"{alias}.name = ?")
                params.append(filters[key])
        if clauses:
            sql += " WHERE " + " AND ".join(clauses)
        sql += " ORDER BY r.data_id"
        return [dict(row) for row in self.conn.execute(sql, params).fetchall()]

    def query_with_features(
        self,
        filters: dict[str, str] | None = None,
        families: list[str] | None = None,
    ) -> list[dict[str, Any]]:
        rows = self.query_reactions(filters)
        family_clause = ""
        params: list[Any] = []
        if families:
            placeholders = ",".join("?" for _ in families)
            family_clause = f" AND family IN ({placeholders})"
            params.extend(families)
        id_fields = {
            "rea": "reactant_id",
            "cat": "catalyst_id",
            "sol": "solvent_id",
            "pro_r": "product_r_id",
            "pro_s": "product_s_id",
        }
        for row in rows:
            for role_prefix, id_field in id_fields.items():
                molecule_id = row.get(id_field)
                if molecule_id is None:
                    continue
                f_rows = self.conn.execute(
                    f"""
                    SELECT family, feature_name, feature_value
                    FROM features
                    WHERE molecule_id = ? {family_clause}
                    """,
                    [molecule_id, *params],
                ).fetchall()
                for feature in f_rows:
                    key = f"{role_prefix}__{feature['family']}__{feature['feature_name']}"
                    row[key] = feature["feature_value"]
        return rows

    def stats(self) -> dict[str, Any]:
        role_counts = {
            row["role"]: int(row["n"])
            for row in self.conn.execute("SELECT role, COUNT(*) AS n FROM molecules GROUP BY role").fetchall()
        }
        n_reactions = self.conn.execute("SELECT COUNT(*) AS n FROM reactions").fetchone()["n"]
        feature_counts = {
            row["family"]: int(row["n"])
            for row in self.conn.execute("SELECT family, COUNT(*) AS n FROM features GROUP BY family ORDER BY family").fetchall()
        }
        ee_values = [float(row["ee"]) for row in self.conn.execute("SELECT ee FROM reactions WHERE ee IS NOT NULL")]
        bins = {"[-1,-0.5)": 0, "[-0.5,0)": 0, "[0,0.5)": 0, "[0.5,1]": 0}
        for ee in ee_values:
            if ee < -0.5:
                bins["[-1,-0.5)"] += 1
            elif ee < 0:
                bins["[-0.5,0)"] += 1
            elif ee < 0.5:
                bins["[0,0.5)"] += 1
            else:
                bins["[0.5,1]"] += 1
        return {
            "db": str(self.db_path.resolve()),
            "n_reactions": int(n_reactions),
            "molecules_by_role": role_counts,
            "features_by_family": feature_counts,
            "ee_histogram": bins,
        }

    def register_custom_meta(
        self,
        family_topic: str,
        description: str | None = None,
        hypothesis: str | None = None,
        script_path: str | None = None,
        verdict: str | None = None,
        n_molecules: int | None = None,
    ) -> None:
        if verdict not in {None, "supported", "rejected", "inconclusive"}:
            raise ValueError(f"Invalid verdict: {verdict}")
        self.conn.execute(
            """
            INSERT INTO custom_features_meta
                (family_topic, description, hypothesis, script_path, verdict, n_molecules, created_at)
            VALUES (?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(family_topic) DO UPDATE SET
                description=excluded.description,
                hypothesis=excluded.hypothesis,
                script_path=excluded.script_path,
                verdict=excluded.verdict,
                n_molecules=excluded.n_molecules,
                created_at=excluded.created_at
            """,
            (family_topic, description, hypothesis, script_path, verdict, n_molecules, utc_now()),
        )

    def run_custom_script(self, script_path: str | Path) -> dict[str, Any]:
        archive_path = _archive_custom_script(script_path)
        os.environ["AHO_DB_PATH"] = str(self.db_path.resolve())
        os.environ["AHO_SKILL_DIR"] = str(SKILL_DIR.resolve())
        sys.path.insert(0, str((SKILL_DIR / "scripts").resolve()))
        namespace = runpy.run_path(str(archive_path))
        result = None
        if callable(namespace.get("run")):
            result = namespace["run"](db_path=str(self.db_path), skill_dir=str(SKILL_DIR))
        elif callable(namespace.get("main")):
            result = namespace["main"](db_path=str(self.db_path), skill_dir=str(SKILL_DIR))
        applied = self._apply_custom_result(result, archive_path) if isinstance(result, dict) else {}
        self.conn.commit()
        return {
            "script": str(archive_path),
            "returned": isinstance(result, dict),
            **applied,
        }

    def _apply_custom_result(self, result: dict[str, Any], archive_path: Path) -> dict[str, Any]:
        n_features = 0
        for feature in result.get("features", []):
            molecule_id = feature.get("molecule_id")
            if molecule_id is None:
                molecule_id = self.molecule_id(feature["molecule_name"], feature["role"])
            if molecule_id is None:
                raise ValueError(f"Could not resolve molecule for custom feature: {feature}")
            self.insert_feature(
                int(molecule_id),
                feature["family"],
                feature["feature_name"],
                float(feature["feature_value"]),
            )
            n_features += 1
        meta = result.get("meta")
        if meta:
            self.register_custom_meta(
                family_topic=meta["family_topic"],
                description=meta.get("description"),
                hypothesis=meta.get("hypothesis"),
                script_path=str(archive_path),
                verdict=meta.get("verdict"),
                n_molecules=meta.get("n_molecules"),
            )
        return {"n_features": n_features, "registered_meta": bool(meta)}


def _read_csv_rows(csv_path: Path) -> list[dict[str, str]]:
    if not csv_path.exists():
        raise FileNotFoundError(f"CSV not found: {csv_path}")
    with csv_path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != CSV_COLUMNS:
            raise ValueError(f"CSV columns must be exactly {CSV_COLUMNS}; got {reader.fieldnames}")
        rows = [dict(row) for row in reader]
    if not rows:
        raise ValueError("CSV has no data rows")
    return rows


def _find_missing_sdfs(rows: list[dict[str, str]], sdf_index: dict[str, str]) -> list[tuple[str, str]]:
    missing = []
    seen = set()
    for row in rows:
        for role, column in ROLE_TO_COLUMN.items():
            name = (row.get(column) or "").strip()
            if not name and role == "PRO_S":
                continue
            key = (role, name)
            if name and key not in seen and name not in sdf_index:
                missing.append(key)
                seen.add(key)
    return missing


def _archive_custom_script(script_path: str | Path) -> Path:
    src = Path(script_path).resolve()
    if not src.exists():
        raise FileNotFoundError(f"Custom script not found: {src}")
    DEFAULT_CUSTOM_SCRIPTS_DIR.mkdir(parents=True, exist_ok=True)
    if src.parent.resolve() == DEFAULT_CUSTOM_SCRIPTS_DIR.resolve():
        return src
    dst = DEFAULT_CUSTOM_SCRIPTS_DIR / f"{utc_now().replace(':', '').replace('+', 'Z')}_{src.name}"
    shutil.copy2(src, dst)
    return dst


def _json_print(obj: object) -> None:
    print(json.dumps(obj, ensure_ascii=False, indent=2, sort_keys=True))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="AHO SQLite utilities")
    parser.add_argument("--db", dest="root_db", default=str(DEFAULT_DB))
    sub = parser.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("init")
    p.add_argument("--db")

    p = sub.add_parser("import-csv")
    p.add_argument("--db")
    p.add_argument("--csv", required=True)
    p.add_argument("--sdf-dir", required=True)

    p = sub.add_parser("query")
    p.add_argument("--db")
    p.add_argument("--data-id")
    p.add_argument("--reactant")
    p.add_argument("--catalyst")
    p.add_argument("--solvent")
    p.add_argument("--product-r")
    p.add_argument("--product-s")
    p.add_argument("--json", action="store_true")

    p = sub.add_parser("query-with-features")
    p.add_argument("--db")
    p.add_argument("--family", action="append", default=[])
    p.add_argument("--json", action="store_true")

    p = sub.add_parser("stats")
    p.add_argument("--db")
    p.add_argument("--json", action="store_true")

    p = sub.add_parser("run-custom-script")
    p.add_argument("--db")
    p.add_argument("--script", required=True)

    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    db_path = args.db or args.root_db
    if args.cmd == "init":
        path = DB.init(db_path)
        print(path)
    elif args.cmd == "import-csv":
        _json_print(DB.import_csv(args.csv, args.sdf_dir, db_path))
    elif args.cmd == "query":
        filters = {
            key: value
            for key, value in {
                "data_id": args.data_id,
                "reactant": args.reactant,
                "catalyst": args.catalyst,
                "solvent": args.solvent,
                "product_r": args.product_r,
                "product_s": args.product_s,
            }.items()
            if value
        }
        with DB(db_path) as db:
            rows = db.query_reactions(filters)
        _json_print(rows) if args.json else print(rows)
    elif args.cmd == "query-with-features":
        with DB(db_path) as db:
            rows = db.query_with_features(families=args.family or None)
        _json_print(rows) if args.json else print(rows)
    elif args.cmd == "stats":
        with DB(db_path) as db:
            stats = db.stats()
        _json_print(stats) if args.json else print(json.dumps(stats, indent=2, sort_keys=True))
    elif args.cmd == "run-custom-script":
        with DB(db_path) as db:
            result = db.run_custom_script(args.script)
        _json_print(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
