#!/bin/bash
# Run callgrind profiling on all verification matrices.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
EXECUTABLE="$PROJECT_ROOT/cpp/build/fracessa"
MATRICES_FILE="$PROJECT_ROOT/python/verification/verification_matrices.json"
OUTPUT_DIR="$SCRIPT_DIR"

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

# Parse every verification matrix; the dataset no longer has an in_use field.
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
echo ""

# Process each matrix
SUCCESS=0
FAILED=0
FAILED_IDS=()

while IFS='|' read -r ID DIMENSION MATRIX_STR IS_CS NUMBER_ESS; do
    OUTPUT_FILE="$OUTPUT_DIR/callgrind.out.$ID"
    
    # Construct CLI format: "<dimension>#<matrix_string>"
    CLI_MATRIX="${DIMENSION}#${MATRIX_STR}"
    
    echo "[$((SUCCESS + FAILED + 1))/$MATRIX_COUNT] Processing matrix ID $ID (dim $DIMENSION, ESS: $NUMBER_ESS)..."
    echo "  Matrix: $CLI_MATRIX"
    echo "  Output: $OUTPUT_FILE"
    
    # Run valgrind with callgrind
    valgrind --tool=callgrind \
             --callgrind-out-file="$OUTPUT_FILE" \
             --dump-instr=yes \
             --collect-jumps=yes \
             --quiet \
             "$EXECUTABLE" -m "$ID" "$CLI_MATRIX" > /dev/null 2>&1
    
    if [ $? -eq 0 ] && [ -f "$OUTPUT_FILE" ]; then
        FILE_SIZE=$(du -h "$OUTPUT_FILE" | cut -f1)
        echo "  ✓ Successfully generated $OUTPUT_FILE ($FILE_SIZE)"
        ((SUCCESS++))
    else
        echo "  ✗ Failed to generate profile for ID $ID"
        FAILED_IDS+=($ID)
        ((FAILED++))
    fi
    echo ""
done <<< "$MATRICES"

echo "=========================================="
echo "Profiling complete!"
echo "  Successful: $SUCCESS"
echo "  Failed: $FAILED"
if [ $FAILED -gt 0 ]; then
    echo "  Failed IDs: ${FAILED_IDS[*]}"
fi
echo "=========================================="
