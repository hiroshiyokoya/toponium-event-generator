#!/bin/bash
# Build and run a 2015 stand-alone event generator inside the Docker image.
#
# usage (inside the container, repo mounted at /work):
#   scripts/run_standalone.sh <bwbw|blvblv|blvbjj|bjjbjj> [tag] [KEY=VALUE ...]
#     KEY=VALUE go to the namelist /TOPBS/ (see standalone/common/steer.f), e.g.
#     scripts/run_standalone.sh blvblv test IGRN=1 ITR=1 INR=0 NCALL=100000 \
#         ITMX1=5 ITMX2=5 NEVENT=500
#   IBSS=1 (integrate with BASES) is always set unless given explicitly.
# The run directory is build/run/<proc>_<tag>/ (git-ignored); the LHE file
# and the logs are copied to results/standalone_<proc>_<tag>/.
set -euo pipefail
PROC=${1:?bwbw|blvblv|blvbjj|bjjbjj}; shift
TAG=${1:-default}; [ $# -gt 0 ] && shift
REPO=${REPO:-/work}
TABLES=${TABLES:-$REPO/legacy/paper-2010/ppblvblv_Grn_v42/GrnFnc_blvblv}

make -s -C "$REPO/standalone" PROC="$PROC" BASES="${BASES:-51}"
EXE=$REPO/build/standalone/$PROC/EG_$PROC.exe
RUN=$REPO/build/run/${PROC}_${TAG}
rm -rf "$RUN"; mkdir -p "$RUN"; cd "$RUN"

# Green-function tables (m_t=173 GeV, Gamma_t=1.49 GeV) and SM parameters
cp "$TABLES"/grnep{1,8}{thre,high}.tbl .
cp "$REPO/standalone/common/param_card_sm.dat" .

# namelist: IBSS=1 by default, then the user's settings
{
  echo "&TOPBS"
  case " $* " in *" IBSS="*) ;; *) echo "  IBSS=1," ;; esac
  for kv in "$@"; do echo "  $kv,"; done
  echo "/"
} > topbs.nml
cat topbs.nml

"$EXE" > run.log 2>&1 || { tail -30 run.log; exit 1; }
grep -E "SIGMA =|STEER|IGRN" run.log || true

OUT=$REPO/results/standalone_${PROC}_${TAG}; mkdir -p "$OUT"
cp topbs.nml run.log "$OUT/"
for f in events*.lhe; do [ -f "$f" ] && gzip -c "$f" > "$OUT/$f.gz"; done
echo "results in $OUT"
