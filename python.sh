#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

FULL_MODE=0
for arg in "$@"; do
    case "$arg" in
        --full)
            FULL_MODE=1
            ;;
        -h|--help)
            echo "Usage: $0 [--full]"
            echo "  default: run FAST speed benchmarks (all matrices except IDs 32 and 34)"
            echo "  --full : run speed benchmarks on ALL verification matrices"
            exit 0
            ;;
        *)
            echo "Unknown option: $arg" >&2
            echo "Usage: $0 [--full]" >&2
            exit 2
            ;;
    esac
done

if [[ "$FULL_MODE" -eq 1 ]]; then
    python3 python/run_matrices.py --full
else
    python3 python/run_matrices.py
fi
