# toponium-event-generator

Event generator for top-quark pair production and decay at hadron colliders, including the
bound-state (toponium) effects near the $t\bar t$ threshold.

Near threshold, the $t\bar t$ pair interacts through the QCD (Coulomb-like) potential before it
decays. The production amplitude is multiplied by a correction built from the momentum-space Green
function $G(E,p)$ of the Schrödinger equation, separately for the colour-singlet and colour-octet
states. The correction is applied to the full $pp\to bW^+\bar bW^-$ (+ decays) matrix elements.
The fully differential cross section is then correct at LO both near threshold and at high energy,
and unweighted events are written in Les Houches format.

- Y. Sumino and H. Yokoya, *Bound-state effects on kinematical distributions of top quarks at
  hadron colliders*, [arXiv:1007.0075](https://arxiv.org/abs/1007.0075), JHEP 09 (2010) 034.
  This repository contains the code of that paper and its 2015 extension.

The repository has two versions of the generator:

- **MadEvent style (2010)**, `legacy/paper-2010/`. These are the MadGraph/MadEvent v4.4.42 process
  directories of the paper, with the Green-function modification. They run like any MadEvent
  process (`run_card.dat`, `param_card.dat`, LHE + banner), so this is the natural choice for
  experimental studies. Final states: `blvblv`, `blvbjj`.
- **Stand-alone (2015)**, `standalone/`. This is a Fortran program with BASES/SPRING and
  MadGraph II/HELAS matrix elements. It has more final states (`bwbw`, `bwblv`, `blvblv`, `blvbjj`,
  `bjjbjj`) and writes LHE in the same layout.

> **Status (2026-10).**
>
> - Both versions build and run in Docker (gfortran 13). All final states produce events.
> - With BASES/SPRING V5.1 from the CPC Program Library, the stand-alone results are
>   bit-for-bit identical to those of the 2015 setup.
> - **Open:** for the same settings, the two versions differ by about 10% in the `blvblv`
>   cross section (#9). Without the threshold correction and K factors the difference is about
>   4%. Until this is understood, treat absolute normalisations with care.
> - The numbers of arXiv:1007.0075 have not been reproduced yet.

## Repository layout

| Path | Content |
|---|---|
| `legacy/paper-2010/ppblvblv_Grn_v42/` | MadEvent 4.4.42 process `pp -> b mu+ vm b~ mu- vm~` with the Green-function modification (`GrnFnc_blvblv/`) |
| `legacy/paper-2010/ppblvbjj_Grn_v42/` | the same for `pp -> b l nu b~ j j` (`GrnFnc_blvbjj/`) |
| `legacy/paper-2010/GrnEP_pro/` | programs that compute the Green-function tables `grnep{1,8}{thre,high}.tbl` |
| `standalone/common/` | sources shared by all final states of the stand-alone version |
| `standalone/<proc>/` | final-state specific sources: `bwbw`, `bwblv`, `blvblv`, `blvbjj`, `bjjbjj` |
| `standalone/bases51/` | compatibility layer for BASES/SPRING V5.1 (BASES itself is not included) |
| `standalone/Makefile` | `make PROC=<proc>` or `make all-procs`; output in `build/standalone/` |
| `scripts/run_paper2010.sh` | builds and runs the 2010 MadEvent process |
| `scripts/run_standalone.sh` | builds and runs a stand-alone generator |
| `scripts/extract_bases51.sh` | extracts the BASES/SPRING V5.1 library from the CPC deck you downloaded |
| `scripts/lhe_masswindow.py` | fractions of LHE events inside resonance mass windows (cross-checks) |
| `docs/versions-2010-2015.md` | how the two versions relate: switches, inputs, outputs, cross-checks |
| `docs/madgraph-status.md` | current status of MG5_aMC and MadGraph7 |
| `docker/Dockerfile` | toolchain image: gfortran (also as `f77`), LHAPDF 6.5.5 + CTEQ6L1 |
| `REFERENCES.md`, `THIRD_PARTY.md` | papers to cite; licences of third-party code and data |

The 2010 sources were restored from the original tarballs, with symbolic links replaced by real
files. Run products (`Events/`, `build/`, `results/`, object files) are git-ignored.

## Quick start (Docker)

Everything runs in a container, with the repository mounted at `/work`.

```bash
docker build -t toponium-eg:dev docker                # once
docker run --rm -v "$PWD:/work" toponium-eg:dev \
  scripts/run_paper2010.sh blvblv 1 1 1 1 500 4 test11
```

On Windows Git Bash, prefix the commands with `MSYS_NO_PATHCONV=1` and use an absolute path such as
`-v "D:/Physics/toponium-event-generator:/work"`. For an interactive shell, use
`docker run --rm -it -v "$PWD:/work" toponium-eg:dev`.

## Switches

Both versions use the same switches (details in
[docs/versions-2010-2015.md](docs/versions-2010-2015.md)):

| Switch | Meaning |
|---|---|
| `IGRN` | 0: no bound-state correction, 1: overall correction prescription, 2: total-angular-momentum prescription |
| `INR` | non-factorizable (non-resonant) diagrams: 0 none, 1 fixed-width scheme (2015 also: 2 overall factor) |
| `ITR` | 0: naive diagrammatic separation, 1: remove the phase-space suppression in the resonant amplitude |
| `IKF` | 1: normalise each subprocess (gg singlet, gg octet, qq̄) to NLO with K factors for LHC 14 TeV |

The Green-function tables in the repository are for $m_t=173$ GeV, $\Gamma_t=1.4911$ GeV and
$\alpha_s(\mu_B)=0.1534$. Other parameters need new tables from `legacy/paper-2010/GrnEP_pro/`
(several hours per table).

## MadEvent style (2010)

```bash
docker run --rm -v "$PWD:/work" toponium-eg:dev \
  scripts/run_paper2010.sh blvblv 1 1 1 1 500 4 test11
# args: <blvblv|blvbjj> IGRN INR ITR IKF [nevents] [ncores] [tag]
```

The script copies the process directory into the container, writes the switches into
`GrnFnc_*/thr.inc`, runs `link.sh` and `bin/generate_events`. The LHE file and the banner go to
`results/<proc>_<tag>/`. Use `RUNCARD="ptl=0 etal=1d2"` to override `run_card.dat` entries.

The default `run_card.dat` of the paper uses 14 TeV, CTEQ6L1, a fixed scale $\mu=m_t$, lepton cuts
($p_T>10$ GeV, $|\eta|<2.5$) and `bwcutoff=15`. Because the W's are declared as decay chains,
`bwcutoff` keeps only events with $|m_{\ell\nu}-m_W|<15\,\Gamma_W$.

## Stand-alone (2015)

```bash
docker run --rm -v "$PWD:/work" toponium-eg:dev \
  scripts/run_standalone.sh blvblv test IGRN=1 ITR=1 INR=0 IKF=1 \
  NCALL=20000 ITMX1=3 ITMX2=3 NEVENT=200
```

`KEY=VALUE` arguments go to `topbs.nml`, which overrides the defaults hard-coded in
`standalone/<proc>/EG_<proc>.f` (see `standalone/common/steer.f`). Without the file, the programs
behave as the 2015 code. The LHE file and the log go to `results/standalone_<proc>_<tag>/`.

### BASES/SPRING (download it yourself)

The stand-alone version needs **BASES/SPRING V5.1** (S. Kawabata, Comput. Phys. Commun. 88 (1995)
309). It is distributed by the CPC Program Library under the CPC licence, which does not allow
redistribution, so it is not included here.

1. Download the zip from Mendeley Data:
   <https://elsevier.digitalcommonsdata.com/datasets/bsdm9422gc/1>
2. Put it into `third_party/` as is (git-ignored).

The build (`scripts/extract_bases51.sh`) takes the library part of the CPC deck. It changes one
line, `REAL FUNCTION DRN*8` becomes `REAL*8 FUNCTION DRN`, for gfortran. It then links it with
`standalone/bases51/bases51_compat.f`, which provides the routines the 2015 code expects:
`BSSETD`, `BSSETP`, `XHSAVE2`, and the CERNLIB `DATIME`, `UCOPY`, `TIMEX`, `TIMEST`.

With V5.1 the cross sections and the LHE files are bit-for-bit identical to those obtained with
the KEK V5.0 copy used in 2015. All five final states were checked. If you have that V5.0 copy,
`BASES=50` uses `third_party/bases50/` instead.

### LHE output compared with MadEvent

The record layout is the same: `<init>`, 12-particle records with t, W, b, leptons, mothers and
colour. The differences are unit event weights, no banner, mass column 0 for b quarks, and
different process ids.

## References

See [REFERENCES.md](REFERENCES.md). It lists the paper of this work, the MadGraph/MadEvent v4
papers (the 2010 code is based on MadGraph/MadEvent **4.4.42**, not on MadGraph5_aMC@NLO), HELAS,
CTEQ6 and BASES/SPRING.

## Licence

[MIT](LICENSE), © 2010–2026 Hiroshi Yokoya. You may use, modify and redistribute the code freely,
as long as the copyright notice and the licence text are kept. If you use it in a publication,
please cite arXiv:1007.0075. If you use the stand-alone version, also cite BASES/SPRING, which
the CPC licence requires.

The MIT licence does not cover the third-party code and data in the MadEvent process directories
of `legacy/paper-2010/`: MadGraph/MadEvent v4, HELAS, CERNLIB routines (GPL) and PDF tables.
These are shipped unchanged and remain under the terms of their authors; see
[THIRD_PARTY.md](THIRD_PARTY.md). BASES/SPRING is not included.

## Author

Hiroshi Yokoya
