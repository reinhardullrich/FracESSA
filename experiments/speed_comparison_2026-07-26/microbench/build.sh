#!/usr/bin/env bash
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
exp=$(cd "$here/.." && pwd)
repo=$(cd "$exp/../.." && pwd)
sysroot=/tmp/fracessa-sysroot/usr

mkdir -p "$here/build"
export CCACHE_DISABLE=1

build_variant() {
  local name=$1
  local source_root=$2
  local project_build=$3

  c++ -std=c++17 -O3 -DNDEBUG -DSPDLOG_COMPILED_LIB -flto -funroll-loops -march=native -mcpu=native \
    -I"$source_root/cpp/include" \
    -I"$source_root/cpp/include/fracessa" \
    -I"$sysroot/include" \
    -I"$repo/cpp/build/_deps/spdlog-src/include" \
    "$here/benchmark_one.cpp" \
    "$project_build/libfracessa_lib.a" \
    "$project_build/_deps/spdlog-build/libspdlog.a" \
    "$sysroot/lib64/libflint.so" \
    "$sysroot/lib64/libmpfr.so" \
    "$sysroot/lib64/libgmp.so" \
    -Wl,-rpath,"$sysroot/lib64" \
    -o "$here/build/benchmark_$name"
}

build_variant current "$repo" "$exp/builds/current"
build_variant exact_pd_only "$exp/sources/exact_pd_only" "$exp/builds/exact_pd_only"
