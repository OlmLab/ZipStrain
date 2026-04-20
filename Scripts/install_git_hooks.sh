#!/usr/bin/env bash
set -euo pipefail

repo_root="$(git rev-parse --show-toplevel)"
hooks_dir="${repo_root}/.githooks"

if [[ ! -d "${hooks_dir}" ]]; then
  echo "Hook directory not found: ${hooks_dir}" >&2
  exit 1
fi

git config core.hooksPath "${hooks_dir}"
echo "Configured Git hooks path: ${hooks_dir}"
echo "The pre-commit hook will now run the full ZipStrain test suite before each commit."
