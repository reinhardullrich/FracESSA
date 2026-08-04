#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

if [[ "$#" -ne 0 ]]; then
    echo "Usage: $0" >&2
    exit 2
fi

./build.sh

CTEST_JOBS="$(nproc)"
ctest --test-dir cpp/build --output-on-failure -j "${CTEST_JOBS}"
PYTHONPATH=python python3 -m unittest discover -s python/tests -p "test_*.py"
