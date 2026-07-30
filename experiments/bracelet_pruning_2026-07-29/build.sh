#!/usr/bin/env bash
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
sysroot=/tmp/fracessa-sysroot/usr
benchmark_source="$repo/experiments/speed_comparison_2026-07-26/microbench/benchmark_one.cpp"
spdlog_build="$repo/cpp/build/_deps/spdlog-build"

common_flags=(
  -std=c++17 -O3 -DNDEBUG -flto -funroll-loops -march=native -mcpu=native
  -DSPDLOG_COMPILED_LIB -Wall -Wextra -Wpedantic
  -I"$sysroot/include" -I"$repo/cpp/build/_deps/spdlog-src/include"
  -I"$repo/cpp/build/_deps/argparse-src/include"
)
link_flags=(
  "$spdlog_build/libspdlog.a"
  "$sysroot/lib64/libflint.so" "$sysroot/lib64/libmpfr.so" "$sysroot/lib64/libgmp.so"
  -pthread "-Wl,-rpath,$sysroot/lib64"
)

build_variant() {
  local name=$1
  local fracessa_source=$2
  shift 2
  local include_flags=("$@" -I"$repo/cpp/include")
  local build="$here/builds/$name"
  mkdir -p "$build/obj"

  c++ "${common_flags[@]}" "${include_flags[@]}" -c "$repo/cpp/src/matrix_parser.cpp" -o "$build/obj/matrix_parser.o"
  c++ "${common_flags[@]}" "${include_flags[@]}" -c "$repo/cpp/src/findeq.cpp" -o "$build/obj/findeq.o"
  c++ "${common_flags[@]}" "${include_flags[@]}" -c "$repo/cpp/src/checkstab.cpp" -o "$build/obj/checkstab.o"
  c++ "${common_flags[@]}" "${include_flags[@]}" -c "$fracessa_source" -o "$build/obj/fracessa.o"
  ar rcs "$build/libfracessa_lib.a" \
    "$build/obj/matrix_parser.o" "$build/obj/findeq.o" \
    "$build/obj/checkstab.o" "$build/obj/fracessa.o"

  c++ "${common_flags[@]}" "${include_flags[@]}" \
    "$repo/cpp/src/main.cpp" "$build/libfracessa_lib.a" \
    "${link_flags[@]}" -o "$build/fracessa"
  c++ "${common_flags[@]}" "${include_flags[@]}" \
    "$benchmark_source" "$build/libfracessa_lib.a" \
    "${link_flags[@]}" -o "$build/benchmark"
}

build_variant current "$repo/cpp/src/fracessa.cpp"
build_variant bracelet "$here/src/fracessa.cpp" -I"$here/include"
