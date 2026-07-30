#!/bin/bash
# Master script to run callgrind profiling and generate summary

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

echo "=========================================="
echo "Callgrind Profiling for All Matrices"
echo "=========================================="
echo ""

# Check prerequisites
echo "Checking prerequisites..."

# Check valgrind
if ! command -v valgrind &> /dev/null; then
    echo "Error: valgrind is not installed."
    echo "Install with: sudo apt-get install valgrind"
    exit 1
fi
echo "✓ valgrind found"

# Check executable
EXECUTABLE="$PROJECT_ROOT/cpp/build/fracessa"
if [ ! -f "$EXECUTABLE" ]; then
    echo "Error: Executable not found at $EXECUTABLE"
    echo "Please build the project first: $PROJECT_ROOT/build.sh"
    exit 1
fi
echo "✓ Executable found: $EXECUTABLE"

# Check matrices file
MATRICES_FILE="$PROJECT_ROOT/python/verification/verification_matrices.json"
if [ ! -f "$MATRICES_FILE" ]; then
    echo "Error: Matrices file not found at $MATRICES_FILE"
    exit 1
fi
echo "✓ Matrices file found: $MATRICES_FILE"

echo ""
echo "Output directory: $SCRIPT_DIR"
echo ""

# Run profiling
echo "=========================================="
echo "Step 1: Running callgrind profiling..."
echo "=========================================="
echo ""

"$SCRIPT_DIR/run_callgrind_all.sh"

if [ $? -ne 0 ]; then
    echo "Error: Profiling failed"
    exit 1
fi

echo ""
echo "=========================================="
echo "Step 2: Generating summary..."
echo "=========================================="
echo ""

"$SCRIPT_DIR/generate_summary.sh"

if [ $? -ne 0 ]; then
    echo "Error: Summary generation failed"
    exit 1
fi

echo ""
echo "=========================================="
echo "Complete!"
echo "=========================================="
echo ""
echo "Summary report: $PROJECT_ROOT/aidocs/profiling/SUMMARY.md"
echo ""
echo "To view individual profiles:"
echo "  kcachegrind $SCRIPT_DIR/callgrind.out.<id>"
echo ""
