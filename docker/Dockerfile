# toponium-event-generator toolchain: gfortran (also as f77) + make + LHAPDF 6.
# - 2010 TopBS (legacy/paper-2010): a MadEvent 4.4.42 process directory built
#   with `f77`; it uses its built-in CTEQ6 PDFs.
# - 2015 stand-alone generators (standalone/): use LHAPDF through the LHAGlue
#   Fortran interface (InitPDFsetByName, evolvePDF, alphasPDF).
# The repository itself is mounted at /work at run time (see README).
FROM ubuntu:24.04

ARG LHAPDF_VERSION=6.5.5
ARG PDF_SETS="cteq6l1"

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        gfortran g++ make perl python3 wget ca-certificates \
        tar gzip unzip file less \
    && rm -rf /var/lib/apt/lists/* \
    && ln -sf "$(command -v gfortran)" /usr/local/bin/f77

# LHAPDF 6 (C++ library + LHAGlue Fortran interface; no Python bindings)
RUN wget -q -O /tmp/lhapdf.tar.gz "https://lhapdf.hepforge.org/downloads/?f=LHAPDF-${LHAPDF_VERSION}.tar.gz" \
    && tar -xzf /tmp/lhapdf.tar.gz -C /tmp \
    && cd /tmp/LHAPDF-${LHAPDF_VERSION} \
    && ./configure --prefix=/opt/lhapdf --disable-python \
    && make -j"$(nproc)" && make install \
    && rm -rf /tmp/lhapdf.tar.gz /tmp/LHAPDF-${LHAPDF_VERSION}
ENV PATH=/opt/lhapdf/bin:$PATH \
    LD_LIBRARY_PATH=/opt/lhapdf/lib \
    LHAPDF_DATA_PATH=/opt/lhapdf/share/LHAPDF

# PDF sets (CTEQ6L1 = "cteq6ll.LHpdf" of LHAPDF 5, used by the 2015 code)
RUN for s in ${PDF_SETS}; do \
        wget -q -O - "https://lhapdfsets.web.cern.ch/lhapdfsets/current/${s}.tar.gz" \
        | tar -xz -C /opt/lhapdf/share/LHAPDF; \
    done

# Old MadEvent sources are not standard-conforming Fortran
ENV FFLAGS_LEGACY="-std=legacy -fallow-argument-mismatch -ffixed-line-length-132 -O"

WORKDIR /work
CMD ["bash"]
