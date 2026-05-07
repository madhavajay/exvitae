#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

if command -v prettier >/dev/null 2>&1; then
  mode="--write"
  args=()
  for arg in "$@"; do
    case "${arg}" in
      --check)
        mode="--check"
        ;;
      --diff)
        ;;
      *)
        args+=("${arg}")
        ;;
    esac
  done
  if ((${#args[@]} == 0)); then
    args=("projects")
  fi

  yaml_files=()
  for path in "${args[@]}"; do
    if [[ -d "${path}" ]]; then
      while IFS= read -r -d '' file; do
        yaml_files+=("${file}")
      done < <(
        find "${path}" \
          \( -name .git -o -name .mypy_cache -o -name .pytest_cache -o -name __pycache__ -o -name target -o -name test-output \) -prune \
          -o -type f \( -iname '*.yaml' -o -iname '*.yml' \) -print0
      )
    else
      case "${path}" in
        *.yaml | *.yml | *.YAML | *.YML)
          yaml_files+=("${path}")
          ;;
      esac
    fi
  done

  if ((${#yaml_files[@]} == 0)); then
    exit 0
  fi

  exec prettier "${mode}" "${yaml_files[@]}"
fi

exec uv run --project "${REPO_ROOT}/tools" python "${SCRIPT_DIR}/tools/format_yaml.py" "$@"
