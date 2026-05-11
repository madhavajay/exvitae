#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLS_PROJECT="$(cd "${SCRIPT_DIR}/../tools" && pwd)"

exec uv run --project "${TOOLS_PROJECT}" python "${SCRIPT_DIR}/tools/test_reports.py" "$@"
