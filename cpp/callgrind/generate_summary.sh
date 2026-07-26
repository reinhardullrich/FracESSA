#!/bin/bash
# Generate comprehensive summary from callgrind output files

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
MATRICES_FILE="$PROJECT_ROOT/python/verification/verification_matrices.json"
OUTPUT_DIR="$SCRIPT_DIR"
SUMMARY_FILE="$PROJECT_ROOT/aidocs/profiling/SUMMARY.md"
mkdir -p "$(dirname "$SUMMARY_FILE")"

echo "Generating comprehensive summary..."

# Get all callgrind output files
CALLGRIND_FILES=("$OUTPUT_DIR"/callgrind.out.*)

if [ ${#CALLGRIND_FILES[@]} -eq 0 ] || [ ! -f "${CALLGRIND_FILES[0]}" ]; then
    echo "Error: No callgrind output files found in $OUTPUT_DIR"
    exit 1
fi

# Load matrix metadata
MATRIX_METADATA=$(python3 << EOF
import json
import sys

with open("$MATRICES_FILE", 'r') as f:
    data = json.load(f)

metadata = {}
for m in data.get('matrices', []):
    metadata[m['id']] = {
        'dimension': m['dimension'],
        'number_ess': m.get('number_ess', 0),
        'is_cs': m.get('is_cs', False)
    }

# Output as JSON for bash to parse
import json
print(json.dumps(metadata))
EOF
)

# Initialize summary
cat > "$SUMMARY_FILE" << 'EOF'
# Callgrind Profiling Summary

Generated from Callgrind profiling of all verification matrices.

EOF

# Extract statistics for each matrix
echo "Extracting statistics from callgrind files..."

STATS_FILE=$(mktemp)
declare -A TOTAL_INSTRUCTIONS
declare -A FILE_SIZES
declare -A EXECUTION_TIMES

for CALLGRIND_FILE in "${CALLGRIND_FILES[@]}"; do
    if [ ! -f "$CALLGRIND_FILE" ]; then
        continue
    fi
    
    # Extract matrix ID from filename
    ID=$(basename "$CALLGRIND_FILE" | sed 's/callgrind\.out\.//')
    
    # Get file size
    FILE_SIZE=$(du -h "$CALLGRIND_FILE" | cut -f1)
    FILE_SIZES[$ID]=$FILE_SIZE
    
    # Extract total instructions from callgrind file
    # The "summary:" line contains the total instruction count
    SUMMARY_LINE=$(grep "^summary:" "$CALLGRIND_FILE" | head -1)
    if [ -n "$SUMMARY_LINE" ]; then
        # Summary format: "summary: <instruction_count>"
        INSTRUCTIONS=$(echo "$SUMMARY_LINE" | awk '{print $2}')
        if [[ "$INSTRUCTIONS" =~ ^[0-9]+$ ]]; then
            TOTAL_INSTRUCTIONS[$ID]=$INSTRUCTIONS
        fi
    fi
    
    # Fallback: try callgrind_annotate if summary line not found
    if [ -z "${TOTAL_INSTRUCTIONS[$ID]}" ] && command -v callgrind_annotate &> /dev/null; then
        ANNOTATE_OUTPUT=$(callgrind_annotate --auto=yes "$CALLGRIND_FILE" 2>/dev/null | head -20)
        TOTAL_LINE=$(echo "$ANNOTATE_OUTPUT" | grep -i "total:" | head -1)
        if [ -n "$TOTAL_LINE" ]; then
            INSTRUCTIONS=$(echo "$TOTAL_LINE" | awk '{print $2}' | tr -d ',')
            if [[ "$INSTRUCTIONS" =~ ^[0-9]+$ ]]; then
                TOTAL_INSTRUCTIONS[$ID]=$INSTRUCTIONS
            fi
        fi
    fi
done

# Generate per-matrix statistics
echo "" >> "$SUMMARY_FILE"
echo "## Per-Matrix Statistics" >> "$SUMMARY_FILE"
echo "" >> "$SUMMARY_FILE"
echo "| Matrix ID | Dimension | ESS Count | Type | Instructions | File Size |" >> "$SUMMARY_FILE"
echo "|-----------|-----------|----------|------|--------------|-----------|" >> "$SUMMARY_FILE"

# Sort by matrix ID
for ID in $(printf '%s\n' "${!TOTAL_INSTRUCTIONS[@]}" | sort -n); do
    METADATA_JSON=$(echo "$MATRIX_METADATA" | python3 -c "import sys, json; data = json.load(sys.stdin); print(json.dumps(data.get('$ID', {})))" 2>/dev/null)
    
    if [ -n "$METADATA_JSON" ] && [ "$METADATA_JSON" != "null" ]; then
        DIMENSION=$(echo "$METADATA_JSON" | python3 -c "import sys, json; print(json.load(sys.stdin).get('dimension', 'N/A'))")
        NUMBER_ESS=$(echo "$METADATA_JSON" | python3 -c "import sys, json; print(json.load(sys.stdin).get('number_ess', 'N/A'))")
        IS_CS=$(echo "$METADATA_JSON" | python3 -c "import sys, json; print('Circular Symmetric' if json.load(sys.stdin).get('is_cs', False) else 'Symmetric')")
    else
        DIMENSION="N/A"
        NUMBER_ESS="N/A"
        IS_CS="N/A"
    fi
    
    INSTRUCTIONS="${TOTAL_INSTRUCTIONS[$ID]:-N/A}"
    FILE_SIZE="${FILE_SIZES[$ID]:-N/A}"
    
    # Format instructions with commas
    if [ "$INSTRUCTIONS" != "N/A" ] && [[ "$INSTRUCTIONS" =~ ^[0-9]+$ ]]; then
        INSTRUCTIONS=$(printf "%'d" "$INSTRUCTIONS")
    fi
    
    echo "| $ID | $DIMENSION | $NUMBER_ESS | $IS_CS | $INSTRUCTIONS | $FILE_SIZE |" >> "$SUMMARY_FILE"
done

# Generate aggregated statistics
echo "" >> "$SUMMARY_FILE"
echo "## Aggregated Statistics" >> "$SUMMARY_FILE"
echo "" >> "$SUMMARY_FILE"

# Calculate min/max/avg instructions
VALID_INSTRUCTIONS=()
for ID in "${!TOTAL_INSTRUCTIONS[@]}"; do
    if [[ "${TOTAL_INSTRUCTIONS[$ID]}" =~ ^[0-9]+$ ]]; then
        VALID_INSTRUCTIONS+=("${TOTAL_INSTRUCTIONS[$ID]}")
    fi
done

if [ ${#VALID_INSTRUCTIONS[@]} -gt 0 ]; then
    STATS=$(printf '%s\n' "${VALID_INSTRUCTIONS[@]}" | python3 << 'PYEOF'
import sys

numbers = [int(line.strip()) for line in sys.stdin if line.strip().isdigit()]
if numbers:
    print(f"Min: {min(numbers):,}")
    print(f"Max: {max(numbers):,}")
    print(f"Avg: {sum(numbers) // len(numbers):,}")
    print(f"Total: {sum(numbers):,}")
PYEOF
)
    
    echo "$STATS" >> "$SUMMARY_FILE"
    echo "" >> "$SUMMARY_FILE"
    
    # Find fastest and slowest matrices
    MIN_ID=""
    MAX_ID=""
    MIN_VAL=999999999999
    MAX_VAL=0
    
    for ID in "${!TOTAL_INSTRUCTIONS[@]}"; do
        if [[ "${TOTAL_INSTRUCTIONS[$ID]}" =~ ^[0-9]+$ ]]; then
            VAL=${TOTAL_INSTRUCTIONS[$ID]}
            if [ $VAL -lt $MIN_VAL ]; then
                MIN_VAL=$VAL
                MIN_ID=$ID
            fi
            if [ $VAL -gt $MAX_VAL ]; then
                MAX_VAL=$VAL
                MAX_ID=$ID
            fi
        fi
    done
    
    if [ -n "$MIN_ID" ]; then
        echo "**Fastest matrix**: ID $MIN_ID with $(printf "%'d" $MIN_VAL) instructions" >> "$SUMMARY_FILE"
    fi
    if [ -n "$MAX_ID" ]; then
        echo "**Slowest matrix**: ID $MAX_ID with $(printf "%'d" $MAX_VAL) instructions" >> "$SUMMARY_FILE"
    fi
    echo "" >> "$SUMMARY_FILE"
fi

# Generate top functions across all matrices
echo "## Top 10 Most Expensive Functions (Aggregated)" >> "$SUMMARY_FILE"
echo "" >> "$SUMMARY_FILE"

# Aggregate function costs across all callgrind files
FUNCTION_STATS=$(mktemp)

for CALLGRIND_FILE in "${CALLGRIND_FILES[@]}"; do
    if [ ! -f "$CALLGRIND_FILE" ]; then
        continue
    fi
    
    # Extract function costs from callgrind file
    # Look for "fn=" lines which indicate function definitions
    grep "^fn=" "$CALLGRIND_FILE" | while IFS= read -r line; do
        # Extract function name and cost
        FUNC_NAME=$(echo "$line" | sed 's/^fn=//' | cut -d' ' -f1)
        COST=$(echo "$line" | awk '{print $NF}')
        
        if [ -n "$FUNC_NAME" ] && [[ "$COST" =~ ^[0-9]+$ ]]; then
            echo "$FUNC_NAME|$COST" >> "$FUNCTION_STATS"
        fi
    done
    
    # Also look for "cfn=" (caller function) and "cfl=" (call file) lines
    # and aggregate costs from "calls=" lines
done

# Aggregate and sort functions
if [ -f "$FUNCTION_STATS" ] && [ -s "$FUNCTION_STATS" ]; then
    AGGREGATED=$(sort "$FUNCTION_STATS" | python3 << 'PYEOF'
import sys
from collections import defaultdict

func_costs = defaultdict(int)

for line in sys.stdin:
    if '|' in line:
        parts = line.strip().split('|')
        if len(parts) == 2:
            func_name = parts[0]
            try:
                cost = int(parts[1])
                func_costs[func_name] += cost
            except ValueError:
                pass

# Sort by total cost
sorted_funcs = sorted(func_costs.items(), key=lambda x: x[1], reverse=True)

for func_name, total_cost in sorted_funcs[:10]:
    print(f"{func_name}|{total_cost}")
PYEOF
)
    
    if [ -n "$AGGREGATED" ]; then
        echo "| Function | Total Instructions |" >> "$SUMMARY_FILE"
        echo "|----------|-------------------|" >> "$SUMMARY_FILE"
        echo "$AGGREGATED" | while IFS='|' read -r FUNC_NAME TOTAL_COST; do
            FORMATTED_COST=$(printf "%'d" "$TOTAL_COST")
            echo "| \`$FUNC_NAME\` | $FORMATTED_COST |" >> "$SUMMARY_FILE"
        done
    else
        echo "*Unable to extract function statistics*" >> "$SUMMARY_FILE"
    fi
else
    echo "*Unable to extract function statistics*" >> "$SUMMARY_FILE"
fi

rm -f "$FUNCTION_STATS"

# Generate per-matrix top functions
echo "" >> "$SUMMARY_FILE"
echo "## Per-Matrix Top 5 Functions" >> "$SUMMARY_FILE"
echo "" >> "$SUMMARY_FILE"

for ID in $(printf '%s\n' "${!TOTAL_INSTRUCTIONS[@]}" | sort -n | head -10); do
    CALLGRIND_FILE="$OUTPUT_DIR/callgrind.out.$ID"
    
    if [ ! -f "$CALLGRIND_FILE" ]; then
        continue
    fi
    
    echo "### Matrix ID $ID" >> "$SUMMARY_FILE"
    echo "" >> "$SUMMARY_FILE"
    
    if command -v callgrind_annotate &> /dev/null; then
        # Use callgrind_annotate to get top functions
        TOP_FUNCS=$(callgrind_annotate --auto=yes "$CALLGRIND_FILE" 2>/dev/null | grep -A 20 "^Profile data file" | grep -E "^\s+[0-9,]+" | head -5)
        
        if [ -n "$TOP_FUNCS" ]; then
            echo "| Function | Instructions |" >> "$SUMMARY_FILE"
            echo "|----------|--------------|" >> "$SUMMARY_FILE"
            echo "$TOP_FUNCS" | while read -r line; do
                COST=$(echo "$line" | awk '{print $1}' | tr -d ',')
                FUNC=$(echo "$line" | sed 's/^[[:space:]]*[0-9,]*[[:space:]]*//')
                if [ -n "$COST" ] && [ -n "$FUNC" ]; then
                    FORMATTED_COST=$(printf "%'d" "$COST" 2>/dev/null || echo "$COST")
                    echo "| \`$FUNC\` | $FORMATTED_COST |" >> "$SUMMARY_FILE"
                fi
            done
        else
            echo "*Unable to extract function statistics*" >> "$SUMMARY_FILE"
        fi
    else
        echo "*callgrind_annotate not available*" >> "$SUMMARY_FILE"
    fi
    
    echo "" >> "$SUMMARY_FILE"
done

echo "Summary generated: $SUMMARY_FILE"
