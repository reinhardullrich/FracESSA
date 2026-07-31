#!/usr/bin/env bash
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
binary="$here/builds/subset_generation_benchmark"
result="$here/results/one_hot_run_cpu2.tsv"
profile_result="$here/results/eager_profile_30_cpu2.tsv"
small_result="$here/results/all_generators_5_median_100_cpu2.tsv"
range_result="$here/results/all_generators_3_to_25_median_100_cpu2.tsv"

mkdir -p "$here/builds" "$here/results"
ccache c++ -std=c++17 -Wall -Wextra -Wpedantic -O3 -DNDEBUG -flto \
  -funroll-loops -march=native -mcpu=native \
  "$here/benchmark.cpp" -o "$binary"

taskset -c 2 "$binary" | tee "$result"
taskset -c 2 "$binary" profile-eager-30 | tee "$profile_result"
taskset -c 2 "$binary" median-5 | tee "$small_result"
taskset -c 2 "$binary" median-3-25 | tee "$range_result"
