#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUST_DIR="${SCRIPT_DIR}/bioscript/rust"
ASSAYS_DIR="${SCRIPT_DIR}/assays"

echo "==> Validating variant YAML files"
cargo run -q -p bioscript-cli --manifest-path "${RUST_DIR}/Cargo.toml" -- validate-variants "${ASSAYS_DIR}"

echo
echo "==> Validating panel YAML files"
cargo run -q -p bioscript-cli --manifest-path "${RUST_DIR}/Cargo.toml" -- validate-panels "${ASSAYS_DIR}"

echo
echo "Validation passed."
