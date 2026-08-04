#!/usr/bin/env bash
set -euo pipefail

if [[ "$#" -ne 4 ]]; then
    echo "Usage: $0 VERSION PLATFORM TRIPLET VCPKG_ROOT" >&2
    exit 2
fi

version="$1"
platform="$2"
triplet="$3"
vcpkg_root="$4"
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd -- "$script_dir/../.." && pwd)"
build_dir="$repo_root/cpp/build-release-$platform"
output="$repo_root/dist/fracessa-$version-$platform"

"$script_dir/install-vcpkg.sh" "$triplet" "$vcpkg_root"

linker_flags=""
if [[ "$platform" == linux-* ]]; then
    linker_flags="-static-libgcc -static-libstdc++"
fi

cmake -B "$build_dir" -S "$repo_root/cpp" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_TOOLCHAIN_FILE="$vcpkg_root/scripts/buildsystems/vcpkg.cmake" \
    -DVCPKG_TARGET_TRIPLET="$triplet" \
    -DVCPKG_OVERLAY_TRIPLETS="$repo_root/.github/triplets" \
    -DBUILD_TESTING=OFF \
    -DFRACESSA_NATIVE_ARCH=OFF \
    -DFRACESSA_BUILD_CLI=ON \
    -DFRACESSA_BUILD_PYTHON=OFF \
    -DCMAKE_EXE_LINKER_FLAGS_RELEASE="$linker_flags"
cmake --build "$build_dir" --config Release --target fracessa -j 4

binary="$build_dir/fracessa"
actual_version="$("$binary" --version)"
if [[ "$actual_version" != "$version" ]]; then
    echo "Release version $version does not match binary version $actual_version" >&2
    exit 1
fi

strip "$binary"
mkdir -p "$repo_root/dist"
cp "$binary" "$output"
