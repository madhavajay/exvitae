#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IMAGE="${EXVITAE_DOCKER_IMAGE:-exvitae:local}"
TARGET="${EXVITAE_DOCKER_TARGET:-runtime}"

if [[ ! -f "${SCRIPT_DIR}/bioscript/rust/Cargo.toml" ]]; then
  echo "bioscript submodule is missing. Run: git submodule update --init --recursive" >&2
  exit 1
fi

docker build \
  -f "${SCRIPT_DIR}/docker/Dockerfile" \
  --target "${TARGET}" \
  -t "${IMAGE}" \
  "${SCRIPT_DIR}" \
  "$@"

echo "Built ${IMAGE}"
