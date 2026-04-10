#!/usr/bin/env bash
set -euo pipefail

# Fetch test data defined in test-data/sources.yaml.
# Only downloads files that are not already present locally.
#
# Usage:
#   ./fetch_test_data.sh                         # download all missing files
#   ./fetch_test_data.sh --check                 # report what's present/missing
#   ./fetch_test_data.sh --dataset 1k-genomes    # only this dataset
#   ./fetch_test_data.sh --dataset 23andme       # only this dataset
#   ./fetch_test_data.sh --list                  # show full manifest
#   ./fetch_test_data.sh --only "*.crai,*.fai"   # only matching files
#   ./fetch_test_data.sh --exclude "*.fa,*.cram"  # skip matching files

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)/test-data"
SOURCES="${SCRIPT_DIR}/sources.yaml"

# --- Helpers -----------------------------------------------------------------

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
CYAN='\033[0;36m'
NC='\033[0m'

info()  { printf "${GREEN}[OK]${NC}    %s\n" "$*"; }
warn()  { printf "${YELLOW}[MISS]${NC}  %s\n" "$*"; }
fetching() { printf "${CYAN}[FETCH]${NC} %s\n" "$*"; }
error() { printf "${RED}[ERR]${NC}   %s\n" "$*"; }

human_size() {
  local bytes=$1
  if   [ "$bytes" -ge 1073741824 ] 2>/dev/null; then echo "$(( bytes / 1073741824 )) GB"
  elif [ "$bytes" -ge 1048576 ]    2>/dev/null; then echo "$(( bytes / 1048576 )) MB"
  elif [ "$bytes" -ge 1024 ]       2>/dev/null; then echo "$(( bytes / 1024 )) KB"
  else echo "${bytes} B"
  fi
}

file_size() {
  stat -f%z "$1" 2>/dev/null || stat --printf="%s" "$1" 2>/dev/null || echo "0"
}

# --- Parse sources.yaml via python and emit tab-separated lines --------------
# Output: dataset<TAB>subdir<TAB>filename<TAB>url

parse_sources() {
  python3 -c "
import yaml, sys
with open(sys.argv[1]) as f:
    cfg = yaml.safe_load(f)
for ds_name, ds in cfg.get('datasets', {}).items():
    for sd_name, sd in ds.get('subdirs', {}).items():
        for f in sd.get('files', []):
            print(f\"{ds_name}\t{sd_name}\t{f['name']}\t{f['url']}\")
" "$SOURCES"
}

# --- CLI ---------------------------------------------------------------------

MODE="download"
FILTER_DATASET=""
ONLY_PATTERNS=""
EXCLUDE_PATTERNS=""

# Check if a filename matches any pattern in a comma-separated list
# Supports: exact name, glob (*), partial substring
matches_any() {
  local filename="$1"
  local patterns="$2"
  [ -z "$patterns" ] && return 1
  local IFS=','
  for pat in $patterns; do
    pat="$(echo "$pat" | sed 's/^ *//;s/ *$//')"
    case "$filename" in
      $pat) return 0 ;;
    esac
  done
  return 1
}

while [ $# -gt 0 ]; do
  case "$1" in
    --check)    MODE="check"; shift ;;
    --list)     MODE="list"; shift ;;
    --dataset)  FILTER_DATASET="$2"; shift 2 ;;
    --only)     ONLY_PATTERNS="$2"; shift 2 ;;
    --exclude)  EXCLUDE_PATTERNS="$2"; shift 2 ;;
    --help|-h)
      cat <<EOF
Usage: $0 [OPTIONS]

Downloads test data files defined in sources.yaml.
Only fetches files that are not already present locally.

Options:
  --check              Report what's present/missing without downloading
  --list               Show the full file manifest with URLs
  --dataset <name>     Only process this dataset (e.g. 1k-genomes, 23andme)
  --only <patterns>    Only include files matching patterns (comma-separated globs)
  --exclude <patterns> Skip files matching patterns (comma-separated globs)
  -h, --help           Show this help

Filter examples:
  --only "*.zip"                  Only 23andMe zip files
  --only "*.fai,*.crai"           Only index files
  --exclude "*.fa,*.cram"         Skip large files, get indexes only
  --dataset 23andme --only "*v5*" Only the v5 23andMe file

Datasets:
  1k-genomes           GRCh38 reference genome + NA06985 high-coverage CRAM
  23andme              23andMe SNP exports (v2-v5) from OpenMined biovault-data

Data directory: ${DATA_DIR}/
Config file:    ${SOURCES}
EOF
      exit 0
      ;;
    *) error "Unknown argument: $1"; exit 1 ;;
  esac
done

if [ ! -f "$SOURCES" ]; then
  error "Sources file not found: ${SOURCES}"
  exit 1
fi

for cmd in curl python3; do
  if ! command -v "$cmd" &>/dev/null; then
    error "$cmd is required but not found."
    exit 1
  fi
done

echo ""
echo "PGx Test Data Fetcher"
echo "====================="
echo "Config: ${SOURCES}"
echo "Data:   ${DATA_DIR}"
[ -n "$FILTER_DATASET" ]   && echo "Dataset: ${FILTER_DATASET}"
[ -n "$ONLY_PATTERNS" ]    && echo "Only:    ${ONLY_PATTERNS}"
[ -n "$EXCLUDE_PATTERNS" ] && echo "Exclude: ${EXCLUDE_PATTERNS}"
echo ""

missing=0
present=0
downloaded=0
failed=0

while IFS=$'\t' read -r dataset subdir filename url; do
  # apply dataset filter
  if [ -n "$FILTER_DATASET" ] && [ "$dataset" != "$FILTER_DATASET" ]; then
    continue
  fi

  # apply --only filter (match against filename or full path)
  if [ -n "$ONLY_PATTERNS" ]; then
    if ! matches_any "$filename" "$ONLY_PATTERNS" && ! matches_any "${dataset}/${subdir}/${filename}" "$ONLY_PATTERNS"; then
      continue
    fi
  fi

  # apply --exclude filter
  if [ -n "$EXCLUDE_PATTERNS" ]; then
    if matches_any "$filename" "$EXCLUDE_PATTERNS" || matches_any "${dataset}/${subdir}/${filename}" "$EXCLUDE_PATTERNS"; then
      continue
    fi
  fi

  label="${dataset}/${subdir}/${filename}"
  dir="${DATA_DIR}/${dataset}/${subdir}"
  full_path="${dir}/${filename}"

  if [ "$MODE" = "list" ]; then
    printf "  %-60s %s\n" "$label" "$url"
    continue
  fi

  if [ -f "$full_path" ]; then
    size=$(file_size "$full_path")
    info "${label} ($(human_size "$size"))"
    present=$((present + 1))
    continue
  fi

  if [ "$MODE" = "check" ]; then
    warn "${label}"
    missing=$((missing + 1))
    continue
  fi

  # download
  mkdir -p "$dir"
  fetching "${label}"
  echo "       ${url}"

  if curl -fSL --progress-bar --retry 3 --retry-delay 5 -o "$full_path" "$url"; then
    size=$(file_size "$full_path")
    info "${label} ($(human_size "$size"))"
    downloaded=$((downloaded + 1))
  else
    error "${label} — download failed"
    rm -f "$full_path"
    failed=$((failed + 1))
  fi

done < <(parse_sources)

echo ""
echo "Summary: ${present} present, ${downloaded} downloaded, ${missing} missing, ${failed} failed"

if [ "$failed" -gt 0 ]; then
  exit 1
fi
