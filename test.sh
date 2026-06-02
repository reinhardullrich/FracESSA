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
            echo "  default: build + CTest core suite + wrapper Python tests + FAST correctness verification matrices"
            echo "  --full : build + CTest core suite + wrapper Python tests + ALL correctness verification matrices"
            exit 0
            ;;
        *)
            echo "Unknown option: $arg" >&2
            echo "Usage: $0 [--full]" >&2
            exit 2
            ;;
    esac
done

./build.sh

CTEST_JOBS="$(nproc)"
ctest --test-dir build --output-on-failure -j "${CTEST_JOBS}" -E "^VerificationMatrix_"
PYTHONPATH=python python3 -m unittest discover -s python/wrapper_v1/tests -p "test_*.py"

if [[ "$FULL_MODE" -eq 1 ]]; then
    ctest --test-dir build --output-on-failure -j "${CTEST_JOBS}" -R "^VerificationMatrix_"
else
    FAST_IDS_RAW="$(python3 python/verification/matrix_selection.py --fast-ids --verification-dir python/verification)"
    if [[ -z "${FAST_IDS_RAW// }" ]]; then
        echo "No fast verification matrices selected." >&2
        exit 1
    fi

    read -r -a FAST_IDS <<< "${FAST_IDS_RAW}"
    FAST_REGEX="$(IFS='|'; echo "${FAST_IDS[*]}")"
    ctest --test-dir build --output-on-failure -j "${CTEST_JOBS}" -R "^VerificationMatrix_(${FAST_REGEX})$"
fi
