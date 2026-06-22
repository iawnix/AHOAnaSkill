#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json

from aar_db import DB
from aar_feat import RDKitFeaturizer
from constants import DEFAULT_DB


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Ingest AHO CSV/SDF data into SQLite")
    parser.add_argument("--csv", required=True)
    parser.add_argument("--sdf-dir", required=True)
    parser.add_argument("--db", default=str(DEFAULT_DB))
    parser.add_argument("--no-features", action="store_true")
    parser.add_argument("--feature-family", default="all")
    parser.add_argument("--json", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    result = {"import": DB.import_csv(args.csv, args.sdf_dir, args.db)}
    if not args.no_features:
        featurizer = RDKitFeaturizer(args.db)
        try:
            result["features"] = featurizer.compute_all_new(family=args.feature_family)
        finally:
            featurizer.close()
    if args.json:
        print(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True))
    else:
        print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
