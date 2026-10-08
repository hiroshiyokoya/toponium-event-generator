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
#   - BASES/SPRING, built once, chosen by BASES=51 (default) or BASES=50:
#       51: official CPC V5.1. Put the zip downloaded from Mendeley Data into
#           third_party/ (not redistributable; scripts/extract_bases51.sh
#           extracts the library) + bases51/bases51_compat.f (BSSETD,
#           BSSETP, XHSAVE2, CERNLIB DATIME/UCOPY/TIMEX/TIMEST)
#       50: the KEK V5.0 copy used in 2015, if present in third_party/bases50
#   - LHAPDF 6 (LHAGlue Fortran interface), from the Docker image
#
# Build inside the Docker image (see README).

ifeq ($(origin FC),default)
  FC := gfortran
endif
# The 2015 code was compiled with g77 (fixed form, 72 columns, static locals)
FFLAGS  ?= -O2 -std=legacy -fno-automatic -finit-local-zero -w
LHAPDF_CONFIG ?= lhapdf-config

TOP     := $(abspath ..)
BUILD   ?= $(TOP)/build/standalone
LIBDIR  := $(BUILD)/lib
DHELAS_SRC := $(TOP)/legacy/paper-2010/ppblvblv_Grn_v42/Source/DHELAS
BASES   ?= 51
BASES50_SRC := $(TOP)/third_party/bases50

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
LIBS    := -L$(LIBDIR) -ldhelas3 -lbases$(BASES) \
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
$(ODIR)/%.exe: $$(addprefix $(ODIR)/,$$(OBJS_$$*)) $(LIBDIR)/libbases$(BASES).a $(LIBDIR)/libdhelas3.a
	$(FC) $(FFLAGS) -o $@ $(filter %.o,$^) $(LIBS)

$(ODIR)/%.o: %.f | $(ODIR)
	$(FC) $(FFLAGS) $(INCS) -c $< -o $@
$(ODIR)/%.o: %.F | $(ODIR)
	$(FC) $(FFLAGS) $(INCS) -c $< -o $@

$(ODIR) $(LIBDIR) $(LIBDIR)/dhelas $(LIBDIR)/bases:
	mkdir -p $@

# ---- libraries (built once, shared by all final states) ----
libs: $(LIBDIR)/libdhelas3.a $(LIBDIR)/libbases$(BASES).a

DHELAS_OBJ := $(patsubst %,$(LIBDIR)/dhelas/%.o,$(basename $(notdir \
              $(wildcard $(DHELAS_SRC)/*.f $(DHELAS_SRC)/*.F))))
$(LIBDIR)/dhelas/%.o: $(DHELAS_SRC)/%.F | $(LIBDIR)/dhelas
	$(FC) -O2 -std=legacy -w -I$(DHELAS_SRC) -c $< -o $@
$(LIBDIR)/dhelas/%.o: $(DHELAS_SRC)/%.f | $(LIBDIR)/dhelas
	$(FC) -O2 -std=legacy -w -I$(DHELAS_SRC) -c $< -o $@
$(LIBDIR)/libdhelas3.a: $(DHELAS_OBJ)
	ar rcs $@ $^

BASESFLAGS := -O2 -std=legacy -w -fno-automatic -fno-second-underscore

# V5.1 (CPC): extracted from the downloaded deck + compatibility layer
$(LIBDIR)/bases51/bases51.f: | $(LIBDIR)
	$(TOP)/scripts/extract_bases51.sh $(LIBDIR)/bases51
$(LIBDIR)/bases51/%.o: $(LIBDIR)/bases51/%.f
	$(FC) $(BASESFLAGS) -c $< -o $@
$(LIBDIR)/bases51/bases51_compat.o: bases51/bases51_compat.f \
                                    $(LIBDIR)/bases51/bases51.f
	$(FC) $(BASESFLAGS) -c $< -o $@
$(LIBDIR)/libbases51.a: $(LIBDIR)/bases51/bases51.o \
                        $(LIBDIR)/bases51/bases51_compat.o
	ar rcs $@ $^

# V5.0 (KEK, as used in 2015): local copy in third_party/bases50
BASES50_OBJ := $(patsubst %,$(LIBDIR)/bases/%.o,$(basename $(notdir \
               $(wildcard $(BASES50_SRC)/*.f))))
$(LIBDIR)/bases/%.o: $(BASES50_SRC)/%.f | $(LIBDIR)/bases
	$(FC) $(BASESFLAGS) -I$(BASES50_SRC) -c $< -o $@
$(LIBDIR)/libbases50.a: $(BASES50_OBJ)
	@test -n "$(BASES50_OBJ)" || { echo "third_party/bases50 not found"; exit 1; }
	ar rcs $@ $^

clean:
	rm -rf $(BUILD)
