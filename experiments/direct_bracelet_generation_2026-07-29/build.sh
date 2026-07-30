#!/usr/bin/env bash
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
source="$here/compare_bracelet_generators.cpp"
flags=(-std=c++17 -Wall -Wextra -Wpedantic)

mkdir -p "$here/builds"
c++ "${flags[@]}" -O3 -DNDEBUG -flto -funroll-loops -march=native -mcpu=native \
  "$source" -o "$here/builds/compare_bracelet_generators"
c++ "${flags[@]}" -O1 -g -fno-omit-frame-pointer \
  -fsanitize=address,undefined "$source" -o "$here/builds/compare_bracelet_generators_san"
