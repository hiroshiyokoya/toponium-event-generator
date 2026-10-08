.SUFFIXES.:
.SUFFIXES.: .x .o .f .F

LD		= g77
FFLAGS		= -O3 -ffixed-line-length-132
INC		= tblgrnep.inc incgrnep.inc filedata.inc

ALPHAS		= alphas_pegasus.o asrgkt.o betafct.o
GAMT		= gamt.o
EXT		= POTQCD.o EI.o
.f.o:
	$(LD) -c $(FFLAGS) $<

all:	GrnEPthre.exe GrnEPhigh.exe plotgrnep.exe

GrnEPthre.exe: GrnEPthre.o Vqcd.o $(INC) $(ALPHAS) $(GAMT) $(EXT)
	$(LD) $(FFLAGS) GrnEPthre.o Vqcd.o \
	$(ALPHAS) $(GAMT) $(EXT) \
	-o GrnEPthre.exe

GrnEPhigh.exe: GrnEPhigh.o Vqcd.o $(INC) $(ALPHAS) $(GAMT) $(EXT)
	$(LD) $(FFLAGS) GrnEPhigh.o Vqcd.o \
	$(ALPHAS) $(GAMT) $(EXT) \
	-o GrnEPhigh.exe

plotgrnep.exe: plotgrnep.o readgrnep.o $(INC)
	$(LD) $(FFLAGS) plotgrnep.o readgrnep.o -o plotgrnep.exe

clean:; @rm -vf *.o a.out \
	GrnEPthre.exe GrnEPhigh.exe plotgrnep.exe \
