#!/bin/bash
# Run the 2010 MadEvent process (pp -> b mu+ nu b~ mu- nu~ / b l nu b j j)
# with the Green-function (bound-state) modification, inside the Docker image.
#
# usage (inside the container, repo mounted at /work):
#   scripts/run_paper2010.sh <blvblv|blvbjj> <igrn> <inr> <itr> <ikf> [nevents] [ncores] [tag]
#     thr.inc switches: IGRN INR ITR IKF  (see GrnFnc_*/thr11.inc), e.g. 1 1 1 1
#     igrn=0 inr=0 is the conventional (no Coulomb/non-resonant) set-up.
# The process directory is copied to $WORKDIR (container-local, default /tmp/tbs)
# because MadEvent relies on symbolic links; Events/ summaries are copied to
# /work/results/<proc>_<tag>/ (git-ignored).
set -euo pipefail
PROC=${1:?blvblv|blvbjj}; IGRN=${2:?}; INR=${3:?}; ITR=${4:?}; IKF=${5:?}
NEV=${6:-1000}; NCORE=${7:-$(nproc)}; TAG=${8:-g${IGRN}n${INR}t${ITR}k${IKF}}
REPO=${REPO:-/work}
WORKDIR=${WORKDIR:-/tmp/tbs}
SRC=$REPO/legacy/paper-2010/ppblvblv_Grn_v42; GRN=GrnFnc_blvblv
if [ "$PROC" = blvbjj ]; then SRC=$REPO/legacy/paper-2010/ppblvbjj_Grn_v42; GRN=GrnFnc_blvbjj; fi
rm -rf "$WORKDIR"; mkdir -p "$WORKDIR"; cp -r "$SRC" "$WORKDIR/proc"; cd "$WORKDIR/proc"
export FFLAGS_LEGACY="${FFLAGS_LEGACY:--std=legacy -fallow-argument-mismatch -ffixed-line-length-132 -O}"
# MadEvent makefiles use $(FFLAGS) set inside; patch them for gfortran
sed -i 's/^FFLAGS *=.*/FFLAGS= -std=legacy -fallow-argument-mismatch -ffixed-line-length-132 -O/' \
    Source/makefile SubProcesses/makefile Source/*/makefile 2>/dev/null || true
for m in SubProcesses/P*/makefile; do
  sed -i 's/^FFLAGS *=.*/FFLAGS= -std=legacy -fallow-argument-mismatch -ffixed-line-length-132 -O/' "$m"
done
# thr.inc for the requested switches
cat > "$GRN/thr_user.inc" <<EOT
	INTEGER IGRN,ITR,INR,IKF
	PARAMETER (IGRN=$IGRN, INR=$INR, ITR=$ITR, IKF=$IKF)
EOT
cat "$GRN"/thr11.inc | sed -n '3,$p' >> "$GRN/thr_user.inc"
ln -sf thr_user.inc "$GRN/thr.inc"
( cd "$GRN" && bash ./link.sh )
sed -i "s/^ *[0-9]* *= *nevents.*/ $NEV = nevents/" Cards/run_card.dat
# optional run_card overrides, e.g. RUNCARD="ptl=0 etal=1d2" (removes lepton cuts)
for kv in ${RUNCARD:-}; do
  key=${kv%%=*}; val=${kv#*=}
  grep -qE "^\s*\S+\s*=\s*${key}(\s|!|$)" Cards/run_card.dat || { echo "no $key in run_card"; exit 1; }
  sed -i -E "s/^\s*\S+(\s*=\s*${key})(\s|!|$)/ ${val}\1\2/" Cards/run_card.dat
  grep -E "=\s*${key}(\s|!|$)" Cards/run_card.dat
done
./bin/generate_events 2 "$NCORE" "$TAG" 2>&1 | tail -40
OUT=$REPO/results/${PROC}_${TAG}; mkdir -p "$OUT"
cp Events/${TAG}_unweighted_events.lhe.gz "$OUT/" 2>/dev/null || true
cp Events/${TAG}_banner.txt "$OUT/" 2>/dev/null || true
cp -r HTML/${TAG}* "$OUT/" 2>/dev/null || true
echo "results in $OUT"
