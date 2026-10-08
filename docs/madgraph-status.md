# MadGraph status (surveyed 2026-10-08)

Why this matters: the 2010 code is a **MadGraph/MadEvent v4.4.42** process
directory. Experimental users today work with MG5_aMC (or soon MadGraph7), so
this note records where the MadGraph family stands, to decide how a
MadEvent-style interface of this generator should evolve.

## MG5_aMC@NLO — current production tool

| Item | Status |
|---|---|
| Latest stable | **v3.8.0** (2026-09-14) |
| LTS | **v3.5.x** (v3.5.16 offered); replaced 2.9.x, which reached end of life in Dec 2025 |
| Development | GitHub `mg5amcnlo/mg5amcnlo`, default branch `3.x` (Launchpad still used for downloads, Q&A and bug reports, active as of Sep 2026) |
| Requirements | Python >= 3.7, gfortran/gcc >= 4.6 |
| Licence | adapted University of Illinois/NCSA (permissive; keep notices, no endorsement) |
| Standard reference | J. Alwall et al., JHEP 07 (2014) 079, arXiv:1405.0301 (about 25,000 citations on INSPIRE) |

So MG5_aMC is still the standard tool, under active maintenance, and is the
natural target for a "MadEvent-style" interface aimed at experimentalists.

## MadGraph7 — next generation (alpha)

- First public **alpha** announced on 2026-09-16 (GitHub `MadGraphTeam/MadGraph7`, release v0.2.0).
- Needs Python >= 3.12 and a C++ compiler. Aims at vectorised/GPU matrix
  elements and phase-space sampling, ML-assisted sampling (MadNIS) and better
  spin correlations in MadSpin.
- LO event generation works. The MG5-style LO workflow is still reachable via
  `output madevent`. The README itself recommends MG5_aMC for stable
  production.
- Plugin API, process-directory format and LHE details are not documented
  yet. Revisit once it leaves alpha.

## Consequences for this repository

1. Keep the **2010 MadEvent v4** version as the reference implementation of
   arXiv:1007.0075 (frozen under `legacy/paper-2010/`).
2. Keep the **stand-alone Fortran** version (2015, all final states) runnable,
   with LHE output in the same format as MadEvent.
3. Candidate for later: port the Green-function correction to an **MG5_aMC
   3.x** process (modified `matrix.f` + tables, or an event-by-event
   reweighting). Do this only when there is demand; MadGraph7 is not yet
   stable enough to target.

Sources: launchpad.net/mg5amcnlo (+announcements), github.com/mg5amcnlo/mg5amcnlo,
github.com/MadGraphTeam/MadGraph7, INSPIRE.
