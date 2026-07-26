#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 SOURCE_DIR BUILD_DIR LOG_FILE" >&2
  exit 2
fi

source_dir=$1
build_dir=$2
log_file=$3
experiment_dir=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$experiment_dir/../.." && pwd)
sysroot=/tmp/fracessa-sysroot/usr
export CCACHE_DISABLE=1
export PKG_CONFIG_SYSROOT_DIR=/tmp/fracessa-sysroot
export PKG_CONFIG_LIBDIR="$sysroot/lib64/pkgconfig"

mkdir -p "$build_dir" "$(dirname "$log_file")"

{
  echo "source=$source_dir"
  echo "build=$build_dir"
  echo "compiler=$(c++ --version | head -n 1)"
  echo "flags=-O3 -DNDEBUG -mcpu=native -march=native -flto + project Release flags"

  cmake -S "$source_dir/cpp" -B "$build_dir" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_FLAGS=-mcpu=native \
    -DCMAKE_CXX_FLAGS=-mcpu=native \
    -DCMAKE_PREFIX_PATH="$sysroot" \
    -DCMAKE_INCLUDE_PATH="$sysroot/include" \
    -DCMAKE_LIBRARY_PATH="$sysroot/lib64" \
    -DFRACESSA_NATIVE_ARCH=ON \
    -DFETCHCONTENT_SOURCE_DIR_SPDLOG="$repo/cpp/build/_deps/spdlog-src" \
    -DFETCHCONTENT_SOURCE_DIR_ARGPARSE="$repo/cpp/build/_deps/argparse-src" \
    -DFETCHCONTENT_SOURCE_DIR_GOOGLETEST="$repo/cpp/build/_deps/googletest-src" \
    -DFETCHCONTENT_SOURCE_DIR_PYBIND11="$repo/cpp/build/_deps/pybind11-src" \
    -DCMAKE_C_COMPILER_LAUNCHER= \
    -DCMAKE_CXX_COMPILER_LAUNCHER= \
    -DCMAKE_EXPORT_COMPILE_COMMANDS=ON

  cmake --build "$build_dir" --target fracessa -j"$(nproc)"
} 2>&1 | tee "$log_file"
