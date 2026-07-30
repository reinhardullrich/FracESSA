#!/usr/bin/env bash
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
out="$here/results/correctness"
mkdir -p "$out"
work=$(mktemp -d "$out/work.XXXXXX")
trap 'rm -rf "$work"' EXIT

printf 'matrix_id\tdimension\tcurrent_candidates\tbracelet_candidates\tcurrent_ess\tbracelet_ess\tcontent_match\tbyte_identical\n' > "$out/summary.tsv"
: > "$out/mismatches.txt"

while IFS=$'\t' read -r id dimension matrix; do
  spec="$dimension#$matrix"
  current="$work/${id}_current.txt"
  bracelet="$work/${id}_bracelet.txt"
  "$here/builds/current/fracessa" --candidates "$spec" > "$current" 2>/dev/null
  "$here/builds/bracelet/fracessa" --candidates "$spec" > "$bracelet" 2>/dev/null

  current_ess=$(head -n1 "$current")
  bracelet_ess=$(head -n1 "$bracelet")
  current_candidates=$(( $(wc -l < "$current") - 2 ))
  bracelet_candidates=$(( $(wc -l < "$bracelet") - 2 ))

  for variant in current bracelet; do
    tail -n +3 "$work/${id}_${variant}.txt" |
      awk -F';' '{print $2 FS $3 FS $4 FS $5 FS $6 FS $8 FS $9 FS $10 FS $11}' |
      LC_ALL=C sort > "$work/${id}_${variant}.normalized"
  done

  content_match=no
  if [[ "$current_ess" == "$bracelet_ess" ]] &&
     [[ "$current_candidates" == "$bracelet_candidates" ]] &&
     cmp -s "$work/${id}_current.normalized" "$work/${id}_bracelet.normalized"; then
    content_match=yes
  else
    {
      echo "matrix $id (n=$dimension)"
      diff -u "$work/${id}_current.normalized" "$work/${id}_bracelet.normalized" |
        head -100 || true
    } >> "$out/mismatches.txt"
  fi

  byte_identical=no
  cmp -s "$current" "$bracelet" && byte_identical=yes
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$id" "$dimension" "$current_candidates" "$bracelet_candidates" \
    "$current_ess" "$bracelet_ess" "$content_match" "$byte_identical" \
    >> "$out/summary.tsv"
  printf 'checked id=%s n=%s candidates=%s match=%s\n' \
    "$id" "$dimension" "$current_candidates" "$content_match"
done < <(
  jq -r '.matrices[] | select(.is_cs == true) | [.id,.dimension,.matrix] | @tsv' \
    "$repo/python/verification/verification_matrices.json"
)

cat "$out/summary.tsv"
if [[ -s "$out/mismatches.txt" ]]; then
  echo '--- mismatches ---'
  cat "$out/mismatches.txt"
  exit 1
fi
