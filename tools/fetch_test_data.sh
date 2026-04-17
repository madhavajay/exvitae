#!/usr/bin/env bash
set -euo pipefail

# Fetch test data defined in tools/sources.yaml.
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
DEFAULT_REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
REPO_ROOT="${BIOSCRIPT_TEST_DATA_REPO_ROOT:-${DEFAULT_REPO_ROOT}}"
DATA_DIR="${REPO_ROOT}/test-data"
DEFAULT_CACHE_ROOT="${HOME}/.bioscript/cache/test-data"
CACHE_ROOT="${BIOSCRIPT_TEST_DATA_CACHE_DIR:-${DEFAULT_CACHE_ROOT}}"
SOURCES="${BIOSCRIPT_TEST_DATA_SOURCES:-${SCRIPT_DIR}/sources.yaml}"

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

ensure_parent_dir() {
  mkdir -p "$(dirname "$1")"
}

is_split_archive_part() {
  case "$1" in
    *.tar.gz.??) return 0 ;;
    *) return 1 ;;
  esac
}

reconstruct_split_archive() {
  local cache_dir="$1"
  local repo_dir="$2"
  local filename="$3"

  is_split_archive_part "$filename" || return 0

  local archive_name="${filename%.[A-Za-z][A-Za-z]}"
  local final_name="${archive_name%.tar.gz}"
  local final_cache_path="${cache_dir}/${final_name}"
  local final_repo_path="${repo_dir}/${final_name}"
  local first_part="${cache_dir}/${archive_name}.aa"

  if [ -f "$final_cache_path" ]; then
    ensure_repo_symlink "$final_cache_path" "$final_repo_path"
    return 0
  fi

  if [ ! -f "$first_part" ]; then
    return 0
  fi

  local parts=()
  local candidate=""
  for candidate in "${cache_dir}/${archive_name}".??; do
    [ -e "$candidate" ] || continue
    parts+=("$candidate")
  done

  if [ "${#parts[@]}" -eq 0 ]; then
    return 0
  fi

  local tmpdir
  tmpdir="$(mktemp -d "${cache_dir}/reconstruct.XXXXXX")"

  if cat "${parts[@]}" | tar tzf - >/dev/null 2>&1; then
    if cat "${parts[@]}" | tar xzf - -C "$tmpdir"; then
      local extracted
      extracted="$(find "$tmpdir" -type f | head -n 1)"
      if [ -n "$extracted" ] && [ -f "$extracted" ]; then
        mv "$extracted" "$final_cache_path"
        ensure_repo_symlink "$final_cache_path" "$final_repo_path"
      fi
    fi
  fi

  rm -rf "$tmpdir"
}

ensure_repo_materialized() {
  local cache_path="$1"
  local repo_path="$2"
  local filename="$3"

  if is_split_archive_part "$filename"; then
    return 0
  fi

  ensure_repo_symlink "$cache_path" "$repo_path"
}

ensure_repo_symlink() {
  local cache_path="$1"
  local repo_path="$2"

  ensure_parent_dir "$repo_path"

  if [ -L "$repo_path" ] && [ "$(readlink "$repo_path")" = "$cache_path" ]; then
    return 0
  fi

  rm -f "$repo_path"
  ln -s "$cache_path" "$repo_path"
}

migrate_repo_file_into_cache() {
  local repo_path="$1"
  local cache_path="$2"

  ensure_parent_dir "$cache_path"
  mv "$repo_path" "$cache_path"
  ensure_repo_symlink "$cache_path" "$repo_path"
}

materialize_local_file() {
  local label="$1"
  local repo_path="$2"
  local cache_path="$3"

  if [ -f "$cache_path" ]; then
    ensure_repo_symlink "$cache_path" "$repo_path"
    return 0
  fi

  if [ -f "$repo_path" ] && [ ! -L "$repo_path" ]; then
    if migrate_repo_file_into_cache "$repo_path" "$cache_path"; then
      info "${label} (migrated into cache)"
      return 0
    fi
    error "${label} — failed to migrate into cache"
    return 1
  fi

  return 1
}

# --- Parse sources.yaml via python and emit tab-separated lines --------------
# Output: dataset<TAB>subdir<TAB>filename<TAB>url

parse_sources() {
  python3 -c "
import os, sys
from urllib.parse import unquote, urlparse

import yaml

def basename(url: str) -> str:
    return os.path.basename(unquote(urlparse(url).path))

with open(sys.argv[1]) as f:
    cfg = yaml.safe_load(f)

for ds_name, ds in cfg.get('datasets', {}).items():
    for sd_name, sd in ds.get('subdirs', {}).items():
        for f in sd.get('files', []):
            print(f\"{ds_name}\t{sd_name}\t{f['name']}\t{f['url']}\")

for sample_name, sample in cfg.get('sample_data_urls', {}).items():
    for key, subdir in (
        ('ref', 'ref'),
        ('ref_index', 'ref'),
        ('aligned', 'aligned'),
        ('aligned_index', 'aligned'),
        ('vcf', 'vcf'),
        ('vcf_index', 'vcf'),
        ('snp', 'snp'),
    ):
        value = sample.get(key)
        if not value:
            continue
        urls = value if isinstance(value, list) else [value]
        for url in urls:
            print(f\"{sample_name}\t{subdir}\t{basename(url)}\t{url}\")
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
  1k-genomes           GRCh38 reference genome + NA06985 CRAM + cleaned NA06985 VCF
  23andme              23andMe SNP exports (v2-v5) from OpenMined biovault-data
  dynamicdna           Dynamic DNA GSAv3-DTC synthetic export on GRCh38

Data directory: ${DATA_DIR}/
Cache directory: ${CACHE_ROOT}/
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
echo "Cache:  ${CACHE_ROOT}"
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
  repo_dir="${DATA_DIR}/${dataset}/${subdir}"
  repo_path="${repo_dir}/${filename}"
  cache_dir="${CACHE_ROOT}/${dataset}/${subdir}"
  cache_path="${cache_dir}/${filename}"

  if [ "$MODE" = "list" ]; then
    printf "  %-60s %s\n" "$label" "$url"
    continue
  fi

  if materialize_local_file "$label" "$repo_path" "$cache_path"; then
    reconstruct_split_archive "$cache_dir" "$repo_dir" "$filename"
    size=$(file_size "$cache_path")
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
  mkdir -p "$cache_dir"
  fetching "${label}"
  echo "       ${url}"

  if curl -fSL --progress-bar --retry 3 --retry-delay 5 -o "$cache_path" "$url"; then
    ensure_repo_materialized "$cache_path" "$repo_path" "$filename"
    reconstruct_split_archive "$cache_dir" "$repo_dir" "$filename"
    size=$(file_size "$cache_path")
    info "${label} ($(human_size "$size"))"
    downloaded=$((downloaded + 1))
  else
    error "${label} — download failed"
    rm -f "$cache_path"
    failed=$((failed + 1))
  fi

done < <(parse_sources)

echo ""
echo "Summary: ${present} present, ${downloaded} downloaded, ${missing} missing, ${failed} failed"

if [ "$failed" -gt 0 ]; then
  exit 1
fi
