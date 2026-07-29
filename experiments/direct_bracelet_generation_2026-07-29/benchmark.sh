#!/usr/bin/env bash
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
binary="$here/builds/compare_bracelet_generators"
target_seconds=${1:-3}
cpu=${CPU:-2}
result="$here/results/paired_${target_seconds}s.tsv"

mkdir -p "$here/results/raw"
printf 'dimension\tbracelets\tfkm_ns\tdirect_ns\tchange_percent\n' > "$result"

for dimension in $(seq 1 24); do
  variants=(fkm direct)
  if (( dimension % 2 == 0 )); then
    variants=(direct fkm)
  fi
  for variant in "${variants[@]}"; do
    taskset -c "$cpu" "$binary" benchmark "$variant" "$dimension" "$target_seconds" \
      > "$here/results/raw/n${dimension}_${variant}.txt"
  done

  bracelets=$(sed -n 's/^bracelets=//p' "$here/results/raw/n${dimension}_fkm.txt")
  fkm_ns=$(sed -n 's/^median_ns=//p' "$here/results/raw/n${dimension}_fkm.txt")
  direct_ns=$(sed -n 's/^median_ns=//p' "$here/results/raw/n${dimension}_direct.txt")
  change=$(awk -v fkm="$fkm_ns" -v direct="$direct_ns" \
    'BEGIN { printf "%.6f", 100.0 * (direct / fkm - 1.0) }')
  printf '%s\t%s\t%.17g\t%.17g\t%s\n' \
    "$dimension" "$bracelets" "$fkm_ns" "$direct_ns" "$change" >> "$result"
  printf 'n=%02d bracelets=%d fkm=%.3fus direct=%.3fus change=%+.2f%%\n' \
    "$dimension" "$bracelets" \
    "$(awk -v value="$fkm_ns" 'BEGIN { print value / 1000.0 }')" \
    "$(awk -v value="$direct_ns" 'BEGIN { print value / 1000.0 }')" \
    "$change"
done
