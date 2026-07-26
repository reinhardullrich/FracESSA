#!/usr/bin/env bash
set -uo pipefail

exp=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$exp/../.." && pwd)
status_file="$exp/results/tag_build_status.tsv"
mkdir -p "$exp/metadata/cmake_patches"

printf 'representative\tequivalent_tags\tcommit\tstatus\n' > "$status_file"
tail -n +2 "$exp/tag_manifest.tsv" | while IFS=$'\t' read -r tag aliases; do
  [[ -n "$tag" ]] || continue
  source_dir="$exp/sources/$tag"
  build_dir="$exp/builds/$tag"
  mkdir -p "$source_dir"
  git -C "$repo" archive "$tag" | tar -x -C "$source_dir"
  commit=$(git -C "$repo" rev-parse "$tag^{commit}")
  cp "$source_dir/cpp/CMakeLists.txt" "$source_dir/cpp/CMakeLists.txt.original"
  "$exp/patch_archived_cmake.py" "$source_dir/cpp/CMakeLists.txt"
  diff -u "$source_dir/cpp/CMakeLists.txt.original" "$source_dir/cpp/CMakeLists.txt" \
    > "$exp/metadata/cmake_patches/${tag}.patch" || true

  echo "=== building $tag ($aliases) ==="
  if "$exp/build_one.sh" "$source_dir" "$build_dir" "$exp/logs/build_${tag}.log"; then
    printf '%s\t%s\t%s\tok\n' "$tag" "$aliases" "$commit" >> "$status_file"
  else
    printf '%s\t%s\t%s\tfailed\n' "$tag" "$aliases" "$commit" >> "$status_file"
  fi
done
