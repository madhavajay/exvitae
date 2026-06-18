#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

BIOSCRIPT_BIN="${BIOSCRIPT_BIN:-${SCRIPT_DIR}/bioscript/bs}"
FORMATTER="${FORMATTER:-${SCRIPT_DIR}/format-yaml.sh}"

SELF_CONTAINED_PACKAGES=(
  assays/risk/APOL1/manifest.yaml
  assays/pgx/pgx-1/manifest.yaml
  assays/pgx/glp1/manifest.yaml
  assays/traits/blood-type/manifest.yaml
  assays/risk/longevity/manifest.yaml
  assays/risk/pcsk9-ldl/panel.yaml
  assays/risk/prostate-cancer-prs/panel.yaml
)

usage() {
  cat >&2 <<'EOF'
usage: ./build-assay-zips.sh [--self-contained] [PACKAGE_OR_DIR...]

Build BioScript package zips.

Default builds all publishable packages under assays/.
--self-contained builds only packages that do not require repo-external build hooks.

Environment:
  BIOSCRIPT_BIN   bioscript binary/wrapper to use
  FORMATTER       YAML formatter command to use
EOF
}

self_contained=0
paths=()
while (($#)); do
  case "$1" in
    --self-contained)
      self_contained=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      paths+=("$1")
      shift
      ;;
  esac
done

if ((self_contained)); then
  if ((${#paths[@]})); then
    echo "--self-contained does not accept explicit package paths" >&2
    exit 2
  fi
  paths=("${SELF_CONTAINED_PACKAGES[@]}")
elif ((${#paths[@]} == 0)); then
  paths=(assays)
fi

exec uv run --with pyyaml python tools/package_bioscript.py \
  --bioscript "${BIOSCRIPT_BIN}" \
  --formatter "${FORMATTER}" \
  "${paths[@]}"
