#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VERSION="$(tr -d '[:space:]' < "${SCRIPT_DIR}/VERSION")"
IMAGE="${EXVITAE_DOCKER_IMAGE:-ghcr.io/madhavajay/exvitae:${VERSION}}"
CACHE_DIR="${BIOSCRIPT_TEST_DATA_CACHE_DIR:-${HOME}/.bioscript/cache/test-data}"

if [[ "${SKIP_DOCKER_BUILD:-0}" != "1" ]]; then
  "${SCRIPT_DIR}/build-docker.sh"
fi

docker run --rm \
  -v "${SCRIPT_DIR}:/workspace" \
  -v "${CACHE_DIR}:${CACHE_DIR}:ro" \
  -w /workspace \
  -e BIOSCRIPT_BIN=/usr/local/bin/bioscript \
  -e BIOSCRIPT_TEST_ROOT=/ \
  -e BIOSCRIPT_TEST_DATA_CACHE_DIR="${CACHE_DIR}" \
  "${IMAGE}" \
  bash ./test.sh "$@"
