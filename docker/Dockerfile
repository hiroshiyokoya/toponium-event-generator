# toponium-event-generator toolchain: gfortran (as f77) + make + ROOT-free base.
# The 2010 TopBS code is a MadEvent 4.4.42 process directory (old Fortran 77,
# built with `f77`), so gfortran is exposed as `f77` and legacy-friendly flags
# are passed via the environment.  LHAPDF is optional (CTEQ is built in).
# The repository itself is mounted at /work at run time (see README).
FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        gfortran g++ make perl python3 wget ca-certificates \
        tar gzip unzip file less \
    && rm -rf /var/lib/apt/lists/* \
    && ln -sf "$(command -v gfortran)" /usr/local/bin/f77

# Old MadEvent sources are not standard-conforming Fortran
ENV FFLAGS_LEGACY="-std=legacy -fallow-argument-mismatch -ffixed-line-length-132 -O"

WORKDIR /work
CMD ["bash"]
