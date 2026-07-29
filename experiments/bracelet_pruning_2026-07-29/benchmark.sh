#!/usr/bin/env bash
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
target_seconds=${1:-3}
cpu=${CPU:-2}
out="$here/results/benchmark"
mkdir -p "$out/raw"

printf 'matrix_id\tdimension\tcurrent_ns\tbracelet_ns\tchange_percent\tess_count\tcandidate_count\n' > "$out/paired_${target_seconds}s.tsv"
index=0
while IFS=$'\t' read -r id dimension matrix; do
  spec="$dimension#$matrix"
  variants=(current bracelet)
  if (( index % 2 != 0 )); then
    variants=(bracelet current)
  fi

  for variant in "${variants[@]}"; do
    taskset -c "$cpu" "$here/builds/$variant/benchmark" \
      "$spec" "$target_seconds" > "$out/raw/${id}_${variant}.txt"
  done

  current_ns=$(sed -n 's/^median_ns=//p' "$out/raw/${id}_current.txt")
  bracelet_ns=$(sed -n 's/^median_ns=//p' "$out/raw/${id}_bracelet.txt")
  ess_count=$(sed -n 's/^ess_count=//p' "$out/raw/${id}_current.txt")
  candidate_count=$(sed -n 's/^candidate_count=//p' "$out/raw/${id}_current.txt")
  change=$(awk -v current="$current_ns" -v bracelet="$bracelet_ns" \
    'BEGIN { printf "%.6f", 100.0 * (bracelet / current - 1.0) }')
  printf '%s\t%s\t%.17g\t%.17g\t%s\t%s\t%s\n' \
    "$id" "$dimension" "$current_ns" "$bracelet_ns" "$change" \
    "$ess_count" "$candidate_count" >> "$out/paired_${target_seconds}s.tsv"
  printf 'id=%02d n=%02d current=%.3fus bracelet=%.3fus change=%+.2f%%\n' \
    "$id" "$dimension" \
    "$(awk -v value="$current_ns" 'BEGIN { print value / 1000.0 }')" \
    "$(awk -v value="$bracelet_ns" 'BEGIN { print value / 1000.0 }')" \
    "$change"
  index=$((index + 1))
done < <(
  jq -r '.matrices[] | select(.is_cs == true) | [.id,.dimension,.matrix] | @tsv' \
    "$repo/python/verification/verification_matrices.json"
)
