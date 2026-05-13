#!/usr/bin/env bash
set -euo pipefail

ROOT="/Users/madhavajay/dev"
REPO="/Users/madhavajay/dev/exvitae-data/exvitae"
PROJECTS="/Users/madhavajay/dev/exvitae-data/projects"
PRIVATE_DNA="/Users/madhavajay/dev/my_private_data/dna"

usage() {
  cat <<'USAGE'
usage:
  ./test-report.sh <data-alias> <assay-or-panel> [--output-dir DIR] [--no-open] [-- EXTRA_ARGS...]
  ./test-report.sh --list

data aliases:
  23andme
  23andme_v5
  carika
  sequencing
  sequencing-vcf
  sequencing-dv-vcf
  carigenetics
  carigenetics-vcf
  NA06985
  NA06985-vcf

assay/panel aliases:
  apol1
  apol1-zip
  pgx-1
  pgx-1-zip

examples:
  ./test-report.sh 23andme apol1
  ./test-report.sh sequencing apol1-zip
  ./test-report.sh carigenetics-vcf /Users/madhavajay/dev/exvitae-data/exvitae/assays/risk/APOL1/APOL1.zip
  ./test-report.sh carika /Users/madhavajay/dev/exvitae-data/projects/pgx-1/manifest.yaml -- --analysis-max-duration-ms 60000
USAGE
}

list_aliases() {
  usage
}

sanitize_name() {
  printf '%s' "$1" | tr '/: ' '---' | tr -cd 'A-Za-z0-9._-'
}

root_relative_path() {
  local path="$1"
  if [[ "$path" == "$ROOT/"* ]]; then
    printf '%s\n' "${path#"$ROOT/"}"
  elif [[ "$path" == /* ]]; then
    printf '%s\n' "$path"
  else
    printf '%s\n' "${REPO#"$ROOT/"}/$path"
  fi
}

resolve_manifest() {
  local value="$1"
  case "$value" in
    apol1)
      printf '%s\n' "$REPO/assays/risk/APOL1/manifest.yaml"
      ;;
    apol1-zip)
      printf '%s\n' "$REPO/assays/risk/APOL1/APOL1.zip"
      ;;
    pgx-1)
      printf '%s\n' "$PROJECTS/pgx-1/manifest.yaml"
      ;;
    pgx-1-zip)
      printf '%s\n' "$PROJECTS/pgx-1/pgx-1.zip"
      ;;
    *)
      printf '%s\n' "$value"
      ;;
  esac
}

configure_data_alias() {
  local alias="$1"
  INPUT_FILE=""
  INPUT_INDEX=""
  INPUT_FORMAT=""
  REFERENCE_FILE=""
  REFERENCE_INDEX=""
  SAMPLE_SEX=""
  NEED_VCF_INDEX="false"
  ALLOW_MD5_MISMATCH="false"
  DEFAULT_ANALYSIS_MAX_DURATION_MS=""

  case "$alias" in
    23andme)
      INPUT_FILE="$PRIVATE_DNA/23andme.com/genome_Madhava_Jay_v4_Full_20250611034825.zip"
      ;;
    23andme_v5)
      ROOT="/Users/madhavajay"
      INPUT_FILE="$REPO/test-data/23andme/v5/hu50B3F5/genome_hu50B3F5_v5_Full.zip"
      DEFAULT_ANALYSIS_MAX_DURATION_MS="30000"
      ;;
    carika)
      INPUT_FILE="$PRIVATE_DNA/dynamicdnalabs.com/carika.txt"
      DEFAULT_ANALYSIS_MAX_DURATION_MS="30000"
      ;;
    sequencing)
      local data_dir="$PRIVATE_DNA/sequencing.com"
      INPUT_FILE="$data_dir/MadhavaJay-SQA28U66-30x-WGS-Sequencing_com.recal.cram"
      INPUT_INDEX="$data_dir/MadhavaJay-SQA28U66-30x-WGS-Sequencing_com.recal.cram.crai"
      REFERENCE_FILE="$data_dir/Homo_sapiens_assembly38.fasta"
      REFERENCE_INDEX="$data_dir/Homo_sapiens_assembly38.fasta.fai"
      SAMPLE_SEX="male"
      ;;
    sequencing-vcf)
      local data_dir="$PRIVATE_DNA/sequencing.com"
      INPUT_FILE="$data_dir/MadhavaJay-SQA28U66-30x-WGS-Sequencing_com-05-20-25.snp-indel.genome.vcf.gz"
      INPUT_INDEX="$INPUT_FILE.tbi"
      INPUT_FORMAT="vcf"
      NEED_VCF_INDEX="true"
      DEFAULT_ANALYSIS_MAX_DURATION_MS="30000"
      ;;
    sequencing-dv-vcf)
      local data_dir="$PRIVATE_DNA/sequencing.com"
      INPUT_FILE="$data_dir/MadhavaJay-SQA28U66-30x-WGS-Sequencing_com.deepvariant.vcf.gz"
      INPUT_INDEX="$INPUT_FILE.tbi"
      INPUT_FORMAT="vcf"
      NEED_VCF_INDEX="true"
      DEFAULT_ANALYSIS_MAX_DURATION_MS="30000"
      ;;
    carigenetics)
      local data_dir="$PRIVATE_DNA/carigenetics.com"
      INPUT_FILE="$data_dir/MadhavaJay-WGS-carigenetics-PBE09234-PBE09980-ENT0001-ONT-20250802-all.cram"
      INPUT_INDEX="$data_dir/MadhavaJay-WGS-carigenetics-PBE09234-PBE09980-ENT0001-ONT-20250802-all.cram.crai"
      REFERENCE_FILE="$data_dir/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
      REFERENCE_INDEX="$data_dir/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai"
      SAMPLE_SEX="male"
      DEFAULT_ANALYSIS_MAX_DURATION_MS="30000"
      ;;
    carigenetics-vcf)
      local data_dir="$PRIVATE_DNA/carigenetics.com"
      INPUT_FILE="$data_dir/MadhavaJay-WGS-carigenetics-PBE09234-PBE09980-ENT0001-ONT-20250802.wf_snp.vcf.gz"
      INPUT_INDEX="$INPUT_FILE.tbi"
      INPUT_FORMAT="vcf"
      NEED_VCF_INDEX="true"
      DEFAULT_ANALYSIS_MAX_DURATION_MS="30000"
      ;;
    NA06985)
      ROOT="/Users/madhavajay"
      local data_dir="$REPO/test-data/1k-genomes"
      INPUT_FILE="$data_dir/aligned/NA06985.final.cram"
      INPUT_INDEX="$data_dir/aligned/NA06985.final.cram.crai"
      REFERENCE_FILE="$data_dir/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa"
      REFERENCE_INDEX="$data_dir/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai"
      ALLOW_MD5_MISMATCH="true"
      DEFAULT_ANALYSIS_MAX_DURATION_MS="30000"
      ;;
    NA06985-vcf)
      ROOT="/Users/madhavajay"
      local data_dir="$REPO/test-data/1k-genomes/vcf"
      INPUT_FILE="$data_dir/NA06985.clean.vcf.gz"
      INPUT_INDEX="$INPUT_FILE.tbi"
      INPUT_FORMAT="vcf"
      NEED_VCF_INDEX="true"
      DEFAULT_ANALYSIS_MAX_DURATION_MS="30000"
      ;;
    *)
      echo "unknown data alias: $alias" >&2
      echo >&2
      usage >&2
      exit 2
      ;;
  esac
}

if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
  usage
  exit 0
fi

if [[ "${1:-}" == "--list" ]]; then
  list_aliases
  exit 0
fi

if [[ $# -lt 2 ]]; then
  usage >&2
  exit 2
fi

DATA_ALIAS="$1"
MANIFEST_ARG="$2"
shift 2

OPEN_REPORT="true"
OUTPUT_DIR=""
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --output-dir)
      if [[ $# -lt 2 ]]; then
        echo "--output-dir requires a value" >&2
        exit 2
      fi
      OUTPUT_DIR="$2"
      shift 2
      ;;
    --no-open)
      OPEN_REPORT="false"
      shift
      ;;
    --)
      shift
      EXTRA_ARGS=("$@")
      break
      ;;
    *)
      echo "unexpected argument: $1" >&2
      echo >&2
      usage >&2
      exit 2
      ;;
  esac
done

configure_data_alias "$DATA_ALIAS"
MANIFEST="$(resolve_manifest "$MANIFEST_ARG")"

if [[ ! -e "$MANIFEST" ]]; then
  echo "assay/panel path does not exist: $MANIFEST" >&2
  exit 1
fi

if [[ ! -f "$INPUT_FILE" ]]; then
  echo "input file does not exist for alias '$DATA_ALIAS': $INPUT_FILE" >&2
  exit 1
fi

if [[ "$NEED_VCF_INDEX" == "true" && ! -f "$INPUT_INDEX" ]]; then
  bcftools index -t "$INPUT_FILE"
fi

if [[ -n "$INPUT_INDEX" && ! -f "$INPUT_INDEX" ]]; then
  echo "input index does not exist for alias '$DATA_ALIAS': $INPUT_INDEX" >&2
  exit 1
fi

if [[ -n "$REFERENCE_FILE" && ! -f "$REFERENCE_FILE" ]]; then
  echo "reference file does not exist for alias '$DATA_ALIAS': $REFERENCE_FILE" >&2
  exit 1
fi

if [[ -n "$REFERENCE_INDEX" && ! -f "$REFERENCE_INDEX" ]]; then
  echo "reference index does not exist for alias '$DATA_ALIAS': $REFERENCE_INDEX" >&2
  exit 1
fi

if [[ -z "$OUTPUT_DIR" ]]; then
  OUTPUT_DIR="${REPO#"$ROOT/"}/test-output/app-$(sanitize_name "$MANIFEST_ARG")-$(sanitize_name "$DATA_ALIAS")"
else
  OUTPUT_DIR="$(root_relative_path "$OUTPUT_DIR")"
fi

cmd=(
  "$REPO/bioscript/bs" report "$MANIFEST"
  --root "$ROOT"
  --input-file "$INPUT_FILE"
  --detect-sex
  --output-dir "$OUTPUT_DIR"
)

if [[ -n "$INPUT_INDEX" ]]; then
  cmd+=(--input-index "$INPUT_INDEX")
fi

if [[ -n "$INPUT_FORMAT" ]]; then
  cmd+=(--input-format "$INPUT_FORMAT")
fi

if [[ -n "$REFERENCE_FILE" ]]; then
  cmd+=(--reference-file "$REFERENCE_FILE")
fi

if [[ -n "$REFERENCE_INDEX" ]]; then
  cmd+=(--reference-index "$REFERENCE_INDEX")
fi

if [[ "$ALLOW_MD5_MISMATCH" == "true" ]]; then
  cmd+=(--allow-md5-mismatch)
fi

if [[ -n "$SAMPLE_SEX" ]]; then
  cmd+=(--sample-sex "$SAMPLE_SEX")
fi

if [[ -n "$DEFAULT_ANALYSIS_MAX_DURATION_MS" ]]; then
  cmd+=(--analysis-max-duration-ms "$DEFAULT_ANALYSIS_MAX_DURATION_MS")
fi

if [[ "$OPEN_REPORT" == "true" ]]; then
  cmd+=(--open)
fi

if [[ ${#EXTRA_ARGS[@]} -gt 0 ]]; then
  cmd+=("${EXTRA_ARGS[@]}")
fi

cd "$REPO"
printf 'running:'
printf ' %q' "${cmd[@]}"
printf '\n'
exec "${cmd[@]}"
