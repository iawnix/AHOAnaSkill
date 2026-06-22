#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="AHOAnaSkill"
SKILL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

read -r -p "Remove conda env '$ENV_NAME'? [y/N] " ans
[[ "$ans" =~ ^[Yy]$ ]] || exit 0

conda env remove -n "$ENV_NAME" -y || true
echo "[i] Env removed. Project directory remains: $SKILL_DIR"
