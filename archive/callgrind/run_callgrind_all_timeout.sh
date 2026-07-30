#!/bin/bash
# Run callgrind profiling on all matrices with a profiling-only 60s timeout.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
EXECUTABLE="$PROJECT_ROOT/cpp/build/fracessa"
MATRICES_FILE="$PROJECT_ROOT/python/verification/verification_matrices.json"
OUTPUT_DIR="$SCRIPT_DIR"
TIMEOUT=60

# Check if valgrind is available
if ! command -v valgrind &> /dev/null; then
    echo "Error: valgrind is not installed. Please install it with: sudo apt-get install valgrind"
    exit 1
fi

# Check if executable exists
if [ ! -f "$EXECUTABLE" ]; then
    echo "Error: Executable not found at $EXECUTABLE"
    exit 1
fi

# Check if matrices file exists
if [ ! -f "$MATRICES_FILE" ]; then
    echo "Error: Matrices file not found at $MATRICES_FILE"
    exit 1
fi

# Parse all verification matrices.
echo "Loading matrices from $MATRICES_FILE..."
MATRICES=$(python3 << EOF
import json
import sys

with open("$MATRICES_FILE", 'r') as f:
    data = json.load(f)

matrices = []
for m in data.get('matrices', []):
    matrices.append({
        'id': m['id'],
        'dimension': m['dimension'],
        'matrix': m['matrix'],
        'is_cs': m.get('is_cs', False),
        'number_ess': m.get('number_ess', 0)
    })

# Output as space-separated values for bash to parse
for m in matrices:
    print(f"{m['id']}|{m['dimension']}|{m['matrix']}|{m['is_cs']}|{m['number_ess']}")
EOF
)

if [ -z "$MATRICES" ]; then
    echo "Error: No matrices found"
    exit 1
fi

# Count matrices
MATRIX_COUNT=$(echo "$MATRICES" | wc -l)
echo "Found $MATRIX_COUNT matrices to profile"
echo "Timeout: ${TIMEOUT}s per matrix"
echo ""

# Process each matrix
SUCCESS=0
FAILED=0
TIMEOUT_COUNT=0
FAILED_IDS=()

while IFS='|' read -r ID DIMENSION MATRIX_STR IS_CS NUMBER_ESS; do
    OUTPUT_FILE="$OUTPUT_DIR/callgrind.out.$ID"
    
    # Construct CLI format: "<dimension>#<matrix_string>"
    CLI_MATRIX="${DIMENSION}#${MATRIX_STR}"
    
    # Run valgrind with callgrind and timeout
    timeout ${TIMEOUT}s valgrind --tool=callgrind \
             --callgrind-out-file="$OUTPUT_FILE" \
             --quiet \
             "$EXECUTABLE" "$CLI_MATRIX" > /dev/null 2>&1
    
    EXIT_CODE=$?
    
    if [ $EXIT_CODE -eq 124 ]; then
        # Timeout occurred
        if [ -f "$OUTPUT_FILE" ]; then
            FILE_SIZE=$(du -h "$OUTPUT_FILE" | cut -f1)
            echo "  ID $ID: Timeout (${TIMEOUT}s) - partial profile saved ($FILE_SIZE)"
            ((TIMEOUT_COUNT++))
            ((SUCCESS++))
        else
            echo "  ID $ID: Timeout (${TIMEOUT}s) - no profile generated"
            FAILED_IDS+=($ID)
            ((FAILED++))
        fi
    elif [ $EXIT_CODE -eq 0 ] && [ -f "$OUTPUT_FILE" ]; then
        FILE_SIZE=$(du -h "$OUTPUT_FILE" | cut -f1)
        echo "  ID $ID: Completed - profile saved ($FILE_SIZE)"
        ((SUCCESS++))
    else
        echo "  ID $ID: Failed (exit code: $EXIT_CODE)"
        FAILED_IDS+=($ID)
        ((FAILED++))
    fi
done <<< "$MATRICES"

echo ""
echo "=========================================="
echo "Profiling Summary:"
echo "  Total matrices: $MATRIX_COUNT"
echo "  Successful: $SUCCESS"
echo "  Timeouts (partial): $TIMEOUT_COUNT"
echo "  Failed: $FAILED"
if [ $FAILED -gt 0 ]; then
    echo "  Failed IDs: ${FAILED_IDS[*]}"
fi
echo "=========================================="
