#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

cmake -B build -S cpp -DCMAKE_BUILD_TYPE=Release
cmake --build build -j"$(nproc)"
