# Stand-alone TopBS generators (2015 version), one Makefile for all final states.
#
#   make PROC=blvblv            # EG_blvblv.exe (+ other programs of that process)
#   make all-procs              # every final state
#
# Shared sources live in common/, process-specific ones in <proc>/.
# Objects and executables go to $(BUILD)/<proc>/; nothing is written to the
# source tree. Libraries:
#   - DHELAS: legacy/paper-2010/.../Source/DHELAS (MG/ME 4.4.42; same routines
#     as the 2015 copy, plus fixes for negative masses), built once
#   - BASES/SPRING 5.0: third_party/bases50 (incl. XHSAVE2), built once
#   - LHAPDF 6 (LHAGlue Fortran interface), from the Docker image
#
# Build inside the Docker image (see README).

FC      ?= gfortran
# The 2015 code was compiled with g77 (fixed form, 72 columns, static locals)
FFLAGS  ?= -O2 -std=legacy -fno-automatic -finit-local-zero -w
LHAPDF_CONFIG ?= lhapdf-config

TOP     := $(abspath ..)
BUILD   ?= $(TOP)/build/standalone
LIBDIR  := $(BUILD)/lib
DHELAS_SRC := $(TOP)/legacy/paper-2010/ppblvblv_Grn_v42/Source/DHELAS
BASES_SRC  := $(TOP)/third_party/bases50

PROCS   := bwbw bwblv blvblv blvbjj bjjbjj
PROC    ?= blvblv

MADG    := steer.o kin_functions.o gamt.o corr.o ggtt.o fvixxz.o grnmod.o \
           readgrnep.o couplings_sm.o
PDF     := ptfl.o pdfport_lhapdf.o

# programs and their objects, per final state
ifeq ($(PROC),bwbw)
  PROGS            := EG_bwbw
  OBJS_EG_bwbw     := EG_bwbw.o intvegas.o kinema4.o sigh_bwbwE.o me_ggbwbwE.o \
                      me_qqbwbwE.o pdfport_lhapdf.o lheout_bwbw.o $(MADG)
endif
ifeq ($(PROC),bwblv)
  PROGS            := main_bwblv tot_bwblv
  UTL5             := kinema5.o sigh_bwblv.o me_ggbwblv.o me_qqbwblv.o
  OBJS_main_bwblv  := main_bwblv.o $(PDF) intvegas.o $(UTL5) $(MADG)
  OBJS_tot_bwblv   := tot_bwblv.o $(PDF) intvegas.o $(UTL5) $(MADG)
endif
ifeq ($(PROC),blvblv)
  PROGS            := EG_blvblv main_blvblv tot_blvblv
  UTL6             := kinema6.o sigh_blvblv.o me_ggblvblv.o me_qqblvblv.o
  OBJS_EG_blvblv   := EG_blvblv.o intvegas.o lheout_12.o kinema6.o \
                      sigh_blvblvE.o me_ggblvblvE.o me_qqblvblvE.o \
                      pdfport_lhapdf.o dilep.o dilep_solve.o $(MADG)
  OBJS_main_blvblv := main_blvblv.o $(PDF) intvegas.o lheout_12.o $(UTL6) $(MADG)
  OBJS_tot_blvblv  := tot_blvblv.o $(PDF) intvegas.o lheout_12.o $(UTL6) $(MADG)
endif
ifeq ($(PROC),blvbjj)
  PROGS            := EG_blvbjj
  OBJS_EG_blvbjj   := EG_blvbjj.o intvegas.o lheout_12.o kinema6.o \
                      sigh_blvbjjE.o me_ggblvbjjE.o me_qqblvbjjE.o \
                      pdfport_lhapdf.o $(MADG)
endif
ifeq ($(PROC),bjjbjj)
  PROGS            := EG_bjjbjj
  OBJS_EG_bjjbjj   := EG_bjjbjj.o intvegas.o lheout.o kinema6.o \
                      sigh_bjjbjjE.o me_ggbjjbjjE.o me_qqbjjbjjE.o \
                      pdfport_lhapdf.o $(MADG)
endif
ifeq ($(PROGS),)
  $(error unknown PROC=$(PROC); choose one of $(PROCS))
endif

ODIR    := $(BUILD)/$(PROC)
LIBS    := -L$(LIBDIR) -ldhelas3 -lbases50 \
           $(shell $(LHAPDF_CONFIG) --ldflags 2>/dev/null || echo -lLHAPDF)

# process directory first, so that its nexternal.inc etc. win
vpath %.f $(PROC) common
vpath %.F $(PROC) common
INCS    := -I$(PROC) -Icommon

.PHONY: all all-procs libs clean
all: libs $(addprefix $(ODIR)/,$(addsuffix .exe,$(PROGS)))

all-procs:
	@for p in $(PROCS); do $(MAKE) --no-print-directory PROC=$$p || exit 1; done

.SECONDEXPANSION:
$(ODIR)/%.exe: $$(addprefix $(ODIR)/,$$(OBJS_$$*)) | libs
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

$(ODIR)/%.o: %.f | $(ODIR)
	$(FC) $(FFLAGS) $(INCS) -c $< -o $@
$(ODIR)/%.o: %.F | $(ODIR)
	$(FC) $(FFLAGS) $(INCS) -c $< -o $@

$(ODIR) $(LIBDIR) $(LIBDIR)/dhelas $(LIBDIR)/bases:
	mkdir -p $@

# ---- libraries (built once, shared by all final states) ----
libs: $(LIBDIR)/libdhelas3.a $(LIBDIR)/libbases50.a

DHELAS_OBJ := $(patsubst %,$(LIBDIR)/dhelas/%.o,$(basename $(notdir \
              $(wildcard $(DHELAS_SRC)/*.f $(DHELAS_SRC)/*.F))))
$(LIBDIR)/dhelas/%.o: $(DHELAS_SRC)/%.F | $(LIBDIR)/dhelas
	$(FC) -O2 -std=legacy -w -I$(DHELAS_SRC) -c $< -o $@
$(LIBDIR)/dhelas/%.o: $(DHELAS_SRC)/%.f | $(LIBDIR)/dhelas
	$(FC) -O2 -std=legacy -w -I$(DHELAS_SRC) -c $< -o $@
$(LIBDIR)/libdhelas3.a: $(DHELAS_OBJ)
	ar rcs $@ $^

BASES_OBJ := $(patsubst %,$(LIBDIR)/bases/%.o,$(basename $(notdir \
             $(wildcard $(BASES_SRC)/*.f))))
$(LIBDIR)/bases/%.o: $(BASES_SRC)/%.f | $(LIBDIR)/bases
	$(FC) -O2 -std=legacy -w -fno-automatic -fno-second-underscore \
	      -I$(BASES_SRC) -c $< -o $@
$(LIBDIR)/libbases50.a: $(BASES_OBJ)
	ar rcs $@ $^

clean:
	rm -rf $(BUILD)
