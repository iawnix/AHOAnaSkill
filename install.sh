#!/usr/bin/env bash
set -euo pipefail

SKILL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_NAME="AHOAnaSkill"

command -v conda >/dev/null 2>&1 || {
    echo "[x] conda is not installed or not on PATH."
    exit 1
}

if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    echo "[i] Updating env '$ENV_NAME'..."
    conda env update -n "$ENV_NAME" -f "$SKILL_DIR/environment.yml" --prune
else
    echo "[i] Creating env '$ENV_NAME'..."
    conda env create -f "$SKILL_DIR/environment.yml"
fi

chmod +x "$SKILL_DIR/bin/aho"
mkdir -p "$SKILL_DIR/data" "$SKILL_DIR/reports/figures" "$SKILL_DIR/reports/scripts"

conda run -n "$ENV_NAME" python "$SKILL_DIR/scripts/aar_db.py" init \
    --db "$SKILL_DIR/data/aho.sqlite"

echo "[i] Checking environment..."
conda run -n "$ENV_NAME" python "$SKILL_DIR/scripts/check_env.py"

echo
echo "[ok] Installed AHOAnaSkill."
echo "Add this to PATH if desired:"
echo "    export PATH=\"$SKILL_DIR/bin:\$PATH\""
