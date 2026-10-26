#!/bin/bash
# Regression check: run the stand-alone generators and the Green-function table
# programs with fixed, short settings and print cross sections and checksums.
# Two runs of this script on different source versions must print identical
# output if the changes do not affect the results (e.g. comments, removal of
# unused declarations).
#
# usage (inside the Docker image, repo mounted at /work):
#   scripts/check/regression.sh > before.txt
#   ... change the sources ...
#   scripts/check/regression.sh > after.txt;  diff before.txt after.txt
set -euo pipefail
REPO=${REPO:-/work}
export BASES=${BASES:-51}
Q="IGRN=1 ITR=1 INR=0 IKF=1 NCALL=20000 ITMX1=3 ITMX2=3 NEVENT=200"

# --- stand-alone event generators (fresh build)
make -s -C "$REPO/standalone" clean >/dev/null
for p in bwbw blvblv blvbjj bjjbjj; do
  "$REPO/scripts/run_standalone.sh" "$p" regression $Q >/dev/null 2>&1
  r=$REPO/results/standalone_${p}_regression
  s=$(grep "SIGMA =" "$r/run.log" | head -1 | awk '{print $3, $6}')
  h=$(zcat "$r"/events*.lhe.gz | md5sum | cut -c1-32)
  echo "EG $p  sigma $s  lhe-md5 $h"
done

# --- Green-function tables: every 200th (thre) / 100th (high) energy point
W=$(mktemp -d); cp -r "$REPO/legacy/paper-2010/GrnEP_pro/." "$W/"; cd "$W"
cp Vqcd_1LFO.f Vqcd.f
sed -i "s/^      ASB = ALPHAS (MUB,MZ,ASZ,NF,NORD)/      ASB = 0.1534D0/" GrnEPthre.f GrnEPhigh.f
sed -i "s/^      DO 50 IE = 0, NE1$/      DO 50 IE = 0, NE1, 200/" GrnEPthre.f
sed -i "s/^      DO 50 IE = 0, NE2$/      DO 50 IE = 0, NE2, 100/" GrnEPhigh.f
make -s LD=gfortran FFLAGS="-O3 -ffixed-line-length-132 -std=legacy -w" \
     GrnEPthre.exe GrnEPhigh.exe >/dev/null 2>&1
./GrnEPthre.exe >/dev/null 2>&1; ./GrnEPhigh.exe >/dev/null 2>&1
for t in grnep1thre grnep8thre grnep1high grnep8high; do
  echo "GrnEP $t.tbl md5 $(md5sum < $t.tbl | cut -c1-32)"
done
rm -rf "$W"
