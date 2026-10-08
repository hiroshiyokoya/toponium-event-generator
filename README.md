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
