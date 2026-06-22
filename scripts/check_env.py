#!/usr/bin/env python3
from __future__ import annotations

import argparse
import sqlite3
import sys


def check(name: str, fn, required: bool = True) -> bool:
    try:
        value = fn()
        print(f"  [ok] {name}: {value}")
        return True
    except Exception as exc:
        marker = "[x]" if required else "[warn]"
        print(f"  {marker} {name}: {exc}")
        return not required


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Check AHOAnaSkill runtime dependencies")
    parser.add_argument("--allow-missing-renderer", action="store_true")
    args = parser.parse_args(argv)

    print(f"Python: {sys.version.split()[0]}")
    ok = all(
        [
            check("rdkit", lambda: __import__("rdkit").__version__),
            check("pandas", lambda: __import__("pandas").__version__),
            check("numpy", lambda: __import__("numpy").__version__),
            check("scipy", lambda: __import__("scipy").__version__),
            check("sqlite3", lambda: sqlite3.sqlite_version),
            check(
                "xyzrender",
                lambda: getattr(__import__("xyzrender"), "__version__", "installed"),
                required=not args.allow_missing_renderer,
            ),
        ]
    )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
