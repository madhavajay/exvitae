#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUST_DIR="${SCRIPT_DIR}/bioscript/rust"
ASSAYS_DIR="${SCRIPT_DIR}/assays"
VALIDATE_RUSTFLAGS="${RUSTFLAGS:-} -Aunused-assignments -Amissing-docs"

run_validator() {
  RUSTFLAGS="${VALIDATE_RUSTFLAGS}" cargo run -q -p bioscript-cli --manifest-path "${RUST_DIR}/Cargo.toml" -- "$@"
}

echo "==> Validating variant YAML files"
run_validator validate-variants "${ASSAYS_DIR}"

echo
echo "==> Validating panel YAML files"
run_validator validate-panels "${ASSAYS_DIR}"

echo
echo "==> Validating assay YAML files"
run_validator validate-assays "${ASSAYS_DIR}"

echo
echo "Validation passed."
