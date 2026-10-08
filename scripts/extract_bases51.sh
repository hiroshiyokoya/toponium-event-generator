#!/bin/bash
# Extract the BASES/SPRING V5.1 library from the CPC Program Library deck.
#
# BASES/SPRING is distributed under the CPC licence (no redistribution), so it
# is not part of this repository. Download it yourself from
#   https://elsevier.digitalcommonsdata.com/datasets/bsdm9422gc/1
#   (S. Kawabata, Comput. Phys. Commun. 88 (1995) 309, catalogue AAFW_v2_0)
# and put the downloaded zip (or aafw_v2_0.gz / aafw_v2_0) into third_party/.
#
# usage: scripts/extract_bases51.sh <output dir>
#   writes <output dir>/bases51.f (library part only, columns 1-72)
set -euo pipefail
OUT=${1:?output dir}
REPO=$(cd "$(dirname "$0")/.." && pwd)
TP=$REPO/third_party
mkdir -p "$OUT"
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT

deck=""
for z in "$TP"/*.zip; do
  [ -f "$z" ] || continue
  if unzip -l "$z" 2>/dev/null | grep -q 'aafw_v2_0'; then
    unzip -q -o -j "$z" '*aafw_v2_0*' -d "$TMP"; break
  fi
done
[ -f "$TP/aafw_v2_0.gz" ] && cp "$TP/aafw_v2_0.gz" "$TMP/"
[ -f "$TP/aafw_v2_0" ] && cp "$TP/aafw_v2_0" "$TMP/"
[ -f "$TMP/aafw_v2_0.gz" ] && gunzip -f "$TMP/aafw_v2_0.gz"
[ -f "$TMP/aafw_v2_0" ] && deck=$TMP/aafw_v2_0
if [ -z "$deck" ]; then
  echo "BASES/SPRING V5.1 (aafw_v2_0) not found in third_party/." >&2
  echo "Download it from https://elsevier.digitalcommonsdata.com/datasets/bsdm9422gc/1" >&2
  exit 1
fi

# library part: from the 'BASES/SPRING library' banner to the test-output data
s=$(grep -n 'BASES/SPRING library' "$deck" | head -1 | cut -d: -f1)
e=$(grep -n '^CPC FREE FORMAT' "$deck" | head -1 | cut -d: -f1)
[ -n "$s" ] && [ -n "$e" ] || { echo "unexpected deck layout" >&2; exit 1; }
sed -n "$((s-1)),$((e-1))p" "$deck" | cut -c1-72 | sed 's/[[:space:]]*$//' \
  | sed 's/REAL FUNCTION DRN\*8(ISEED)/REAL*8 FUNCTION DRN(ISEED)/' \
  > "$OUT/bases51.f"
grep -q 'REAL\*8 FUNCTION DRN(ISEED)' "$OUT/bases51.f" \
  || { echo "DRN fix not applied" >&2; exit 1; }
echo "extracted $(grep -ciE '^ +subroutine' "$OUT/bases51.f") subroutines to $OUT/bases51.f"
