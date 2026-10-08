# toponium-event-generator

Event generator for top-quark pair production at hadron colliders including
bound-state (toponium) effects near the threshold (code name: **TopBS**).

Based on: Y. Sumino and H. Yokoya, *Bound-state effects on kinematical
distributions of top quarks at hadron colliders*, JHEP 09 (2010) 034,
[arXiv:1007.0075](https://arxiv.org/abs/1007.0075).

> Status: work in progress. Legacy code (2010 paper version, 2015 extended
> version with all final states) is being restored from backups and made to run
> in Docker. See the tracking issue.

## Docker

All code is built and run inside Docker (Ubuntu 24.04 + gfortran; `f77` is
provided as an alias of `gfortran` because the 2010 TopBS code is a
MadEvent 4.4.42 process directory).

```bash
docker build -t toponium-eg:dev docker
docker run --rm -it -v "$PWD":/work toponium-eg:dev
```

Inside the container, `$FFLAGS_LEGACY` holds the flags needed for the old
Fortran 77 sources (`-std=legacy -fallow-argument-mismatch
-ffixed-line-length-132 -O`), e.g. `make FFLAGS="$FFLAGS_LEGACY"` in `Source/`.

## Layout

| Path | Content |
|---|---|
| `legacy/paper-2010/ppblvblv_Grn_v42/` | MadEvent 4.4.42 process `pp -> b mu+ vm b~ mu- vm~` (via `t t~`) with the Green-function modification (`GrnFnc_blvblv/`) |
| `legacy/paper-2010/ppblvbjj_Grn_v42/` | same for `pp -> b l nu b~ j j` (`GrnFnc_blvbjj/`) |
| `legacy/paper-2010/GrnEP_pro/` | programs that compute the momentum-space Green function tables (`grnep{1,8}{thre,high}.tbl`) |
| `standalone/` | 2015 stand-alone generators (BASES/SPRING), all final states: `bwbw`, `bwblv`, `blvblv`, `blvbjj`, `bjjbjj`. Shared sources in `standalone/common/` |
| `scripts/run_standalone.sh` | build + run wrapper for `standalone/` (inside Docker) |
| `docs/versions-2010-2015.md` | how the two versions relate (switches, inputs, outputs, cross-checks) |
| `scripts/run_paper2010.sh` | build + run wrapper (inside Docker) |
| `docker/` | Docker image |
| `REFERENCES.md` | papers to cite (this work, MadGraph/MadEvent, ...) |

The 2010 sources were restored from the original tarballs (symbolic links
replaced by real files). Run products (`Events/`, `results/`, object files)
are git-ignored.

## Quick test (2010 version)

```bash
docker build -t toponium-eg:dev docker
docker run --rm -v "$PWD":/work toponium-eg:dev \
  scripts/run_paper2010.sh blvblv 1 1 1 1 500 4 test11
# args: <blvblv|blvbjj> IGRN INR ITR IKF [nevents] [ncores] [tag]
```

`IGRN`: 0 no Coulomb correction / 1 overall prescription / 2 total-angular-momentum
prescription; `INR`: non-factorizable diagrams; `ITR`: PS-suppression removal;
`IKF`: NLO normalisation. Output (LHE + banner) goes to `results/<proc>_<tag>/`.

## Quick test (2015 stand-alone version)

```bash
docker run --rm -v "$PWD":/work toponium-eg:dev \
  scripts/run_standalone.sh blvblv test IGRN=1 ITR=1 INR=0 IKF=1 \
  NCALL=20000 ITMX1=3 ITMX2=3 NEVENT=200
```

`KEY=VALUE` arguments are written to `topbs.nml`, which overrides the defaults
hard-coded in `standalone/<proc>/EG_<proc>.f` (see `standalone/common/steer.f`).
`make -C standalone all-procs` builds every final state into `build/standalone/`.

BASES/SPRING 5.0 is needed in `third_party/bases50/`. It is not distributed
here while its licence is being checked (#6).

## Which version to use

- **MadEvent style (2010)**: the reference implementation of arXiv:1007.0075.
  It is run like any MadEvent process (`run_card.dat`, `param_card.dat`,
  LHE + banner), so it is the natural choice for experimental studies.
- **Stand-alone (2015)**: more final states and lighter dependencies. It writes
  LHE in the same layout (`<init>`, 12-particle records with t, W, b, leptons,
  mothers and colour). Small differences remain: unit event weights, no
  banner, mass column 0 for b quarks, and different process ids.

See [docs/versions-2010-2015.md](docs/versions-2010-2015.md).

## References

MadGraph/MadEvent version used, and papers to cite: see
[REFERENCES.md](REFERENCES.md). The 2010 code is based on **MadGraph/MadEvent
v4 (4.4.42)**, not on MadGraph5_aMC@NLO.

## Verified

`pp -> b mu+ vm b~ mu- vm~` at 14 TeV (default `run_card.dat`, 500 events,
gfortran 13 in Docker): total cross section **7.98 pb** with the bound-state
(Green-function) correction (`IGRN=1 INR=1 ITR=1 IKF=1`) vs **7.14 pb** without
it (`IGRN=0 INR=0 ITR=0 IKF=1`), i.e. +12%. (Smoke test only; not yet compared
with the numbers in arXiv:1007.0075.)
