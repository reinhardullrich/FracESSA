#!/bin/bash
# Script to regenerate codebase_cpp.txt and codebase_python.txt
# Usage: ./generate_codebase.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

cd "$PROJECT_ROOT"

# Create C++ codebase file
echo "Generating codebase_cpp.txt..."
cat > analysis/codebase_cpp.txt << 'EOF'
=== cpp/CMakeLists.txt ===
EOF
cat cpp/CMakeLists.txt >> analysis/codebase_cpp.txt
echo "" >> analysis/codebase_cpp.txt
echo "=== .github/workflows/release.yml ===" >> analysis/codebase_cpp.txt
cat .github/workflows/release.yml >> analysis/codebase_cpp.txt

# Add all header files
find cpp/include -name "*.hpp" | sort | while read file; do
    echo "" >> analysis/codebase_cpp.txt
    echo "=== $file ===" >> analysis/codebase_cpp.txt
    cat "$file" >> analysis/codebase_cpp.txt
done

# Add all source files
find cpp/src -name "*.cpp" | sort | while read file; do
    echo "" >> analysis/codebase_cpp.txt
    echo "=== $file ===" >> analysis/codebase_cpp.txt
    cat "$file" >> analysis/codebase_cpp.txt
done

echo "✅ codebase_cpp.txt generated"

# Create Python codebase file
echo "Generating codebase_python.txt..."
cat > analysis/codebase_python.txt << 'EOF'
EOF

# Add all Python files (excluding results/ directory)
find python -name "*.py" ! -path "*/results/*" | sort | while read file; do
    echo "=== $file ===" >> analysis/codebase_python.txt
    cat "$file" >> analysis/codebase_python.txt
    echo "" >> analysis/codebase_python.txt
done

# Add verification_matrices.json
echo "=== python/verification/verification_matrices.json ===" >> analysis/codebase_python.txt
cat python/verification/verification_matrices.json >> analysis/codebase_python.txt
echo "" >> analysis/codebase_python.txt

# Add baseline_result.json
echo "=== python/verification/baseline_result.json ===" >> analysis/codebase_python.txt
cat python/verification/baseline_result.json >> analysis/codebase_python.txt

echo "✅ codebase_python.txt generated"
echo ""
echo "Files generated:"
echo "  - analysis/codebase_cpp.txt"
echo "  - analysis/codebase_python.txt"

