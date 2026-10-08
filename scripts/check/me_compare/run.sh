#!/bin/bash
# Compare the singlet gg matrix elements of the 2010 (MadEvent) and 2015
# (stand-alone) versions at the same phase-space points (issue #9).
# Run inside the Docker image with the repository at /work, after
#   scripts/run_standalone.sh blvblv <tag> ...   (provides the LHE points)
# usage: scripts/check/me_compare/run.sh <events.lhe.gz> [npoints] [idprup]
# On ITR=0 events, the mean of ITR1/ITR0 estimates sigma(ITR=1)/sigma(ITR=0)
# independently of the integrator.
set -euo pipefail
LHE=${1:?events.lhe.gz}; NPT=${2:-200}; IDP=${3:-}
REPO=/work; HERE=$REPO/scripts/check/me_compare
W=/tmp/mecmp; rm -rf $W; mkdir -p $W; cd $W
python3 -I $HERE/make_points.py "$LHE" $W/points.dat "$NPT" $IDP
TBL=$REPO/legacy/paper-2010/ppblvblv_Grn_v42/GrnFnc_blvblv
FL="-O2 -std=legacy -w -fallow-argument-mismatch"

# ---- 2010: one build per ITR (thr.inc is compile-time)
for ITR in 0 1; do
  D=$W/v2010_itr$ITR; cp -r $REPO/legacy/paper-2010/ppblvblv_Grn_v42 $D
  G=$D/GrnFnc_blvblv
  { printf '\tINTEGER IGRN,ITR,INR,IKF\n\tPARAMETER (IGRN=1, INR=0, ITR=%s, IKF=0)\n' $ITR
    sed -n '3,$p' $G/thr11.inc; } > $G/thr_user.inc
  ln -sf thr_user.inc $G/thr.inc
  ( cd $G && bash ./link.sh >/dev/null )
  ( cd $D/Source && make -s FFLAGS="$FL -ffixed-line-length-132" >/dev/null 2>&1 )
  P1=$D/SubProcesses/P1_gg_mu+vmmu-vmxbbx
  ( cd $P1 && gfortran $FL -ffixed-line-length-132 -I. -c matrix.f -o matrix.o \
      && gfortran $FL -ffixed-line-length-132 -I. $HERE/t2010.f matrix.o \
         $D/Source/kin_functions.f -L$D/lib -lmodel -ldhelas3 -o $W/t2010_itr$ITR )
done
mkdir -p $W/Cards; cp $D/Cards/param_card.dat $W/param_card.dat; cp $D/Cards/param_card.dat $W/Cards/
cp $TBL/grnep*.tbl $W/

# ---- 2015: one build, switches at run time
S=$REPO/standalone; L=$REPO/build/standalone/lib
make -s -C $S PROC=blvblv >/dev/null
gfortran $FL -I$S/blvblv -I$S/common $HERE/t2015.f $S/blvblv/me_ggblvblvE.f \
  $S/common/{ggtt,grnmod,readgrnep,corr,kin_functions,couplings_sm}.f \
  $S/common/fvixxz.F $S/blvblv/kinema6.f \
  -L$L -ldhelas3 -o $W/t2015
cp $S/common/param_card_sm.dat $W/

# ---- run and compare
./t2010_itr0 > a0.txt; ./t2010_itr1 > a1.txt
./t2015 1 0 0 2>/dev/null | grep -E '^ +[0-9]+ ' > b0.txt
./t2015 1 1 0 2>/dev/null | grep -E '^ +[0-9]+ ' > b1.txt
python3 -I - <<'EOF'
import statistics as st
def r(f):
    out = []
    for l in open(f):
        t = l.split()
        if len(t) == 2 and t[0].isdigit():
            out.append(float(t[1]))
    return out
a0, a1, b0, b1 = r("a0.txt"), r("a1.txt"), r("b0.txt"), r("b1.txt")
n = min(map(len, (a0, a1, b0, b1)))
abs0 = [b0[i] / a0[i] for i in range(n) if a0[i] > 0]
ra = [a1[i] / a0[i] for i in range(n) if a0[i] > 0]
rb = [b1[i] / b0[i] for i in range(n) if b0[i] > 0]
dr = [rb[i] / ra[i] for i in range(min(len(ra), len(rb)))]
q = lambda x: f"median {st.median(x):.6f}  min {min(x):.6f}  max {max(x):.6f}"
print(f"points: {n}")
print("2015/2010 |M|^2, ITR=0      :", q(abs0))
print("ITR1/ITR0 in 2010            :", q(ra))
print("ITR1/ITR0 in 2015            :", q(rb))
print("(2015 ratio)/(2010 ratio)    :", q(dr))
m = lambda x: f"{st.mean(x):.4f} +- {st.stdev(x)/len(x)**0.5:.4f}"
print("mean ITR1/ITR0 (2010 ME)     :", m(ra))
print("mean ITR1/ITR0 (2015 ME)     :", m(rb))
EOF
