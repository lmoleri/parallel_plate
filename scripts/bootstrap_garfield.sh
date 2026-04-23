#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"
workspace_root="$(cd "${repo_root}/../.." && pwd)"

garfield_source_dir="${GARFIELD_SOURCE_DIR:-${workspace_root}/vendor/garfieldpp}"
garfield_build_dir="${GARFIELD_BUILD_DIR:-${workspace_root}/vendor/garfieldpp-build}"
garfield_install_prefix="${GARFIELD_INSTALL_PREFIX:-${workspace_root}/local/garfield}"

if [[ ! -d "${garfield_source_dir}" ]]; then
  echo "Garfield++ source directory not found at: ${garfield_source_dir}" >&2
  echo "Set GARFIELD_SOURCE_DIR or place a shared checkout at ../../vendor/garfieldpp." >&2
  exit 1
fi

cmake -S "${garfield_source_dir}" \
  -B "${garfield_build_dir}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="${garfield_install_prefix}" \
  -DGARFIELD_WITH_CUDA=OFF \
  -DGARFIELD_WITH_DEGRADE=OFF \
  -DGARFIELD_WITH_EXAMPLES=OFF \
  -DGARFIELD_WITH_TESTS=OFF \
  -DGARFIELD_WITH_GSL=ON

cmake --build "${garfield_build_dir}" --target install --parallel "${JOBS:-4}"
