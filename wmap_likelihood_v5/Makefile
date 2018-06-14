#SYSTYPE="IUCAA_CETUS"
#SYSTYPE="IUCAA_CRAY"
SYSTYPE="IUCAA_PERSEUS"


F90=ifort
FFLAGS+= -openmp -O2 -ip -W0 -WB -fpp2 -vec_report0

ifeq ($(SYSTYPE),"IUCAA_PERSEUS")
CFITSIO=/data1/sysadmin/cjayanti/Software/cfitsio/
LAPACK=/data1/sysadmin/cjayanti/Software/lapack-3.4.2
MKL=/opt/software/intel/mkl/lib/intel64/
LGCC=/usr/lib/gcc/x86_64-redhat-linux/4.4.6/
endif


INCS   = -I. -I$(CFITSIO)/include

WMAPLIB = libwmap9.a
OBJS = read_archive_map.o \
	read_fits.o \
	healpix_types.o \
	br_mod_dist.o \
	WMAP_9yr_options.o \
	WMAP_9yr_util.o \
	WMAP_9yr_gibbs.o \
	WMAP_9yr_tt_pixlike.o \
	WMAP_9yr_tt_beam_ptsrc_chisq.o \
	WMAP_9yr_teeebb_pixlike.o \
	WMAP_9yr_tetbeebbeb_pixlike.o \
	WMAP_9yr_likelihood.o
#
#			General Commands.
#
DIFF = diff -w
RM = rm -f

# See the CHANGES files for a description of these options
WMAPFLAGS  = -DOPTIMIZE
#WMAPFLAGS += -DUSE_LOWELL_TBEB    # turns on low-l maxlike TB/EB
#WMAPFLAGS += -DUSE_HIGHELL_TB     # turns on high-l master TB
#WMAPFLAGS += -DFASTERTT           # speed up matrix low-l TT
#WMAPFLAGS += -DTIMING             # print out timing stats for profiling

#
#			Compiler/linker configuration.  Several samples
#			are supplied.
#

FFLAGS += -O2 -fpic $(WMAPFLAGS) 
#LIBS = -L. -L$(MKL)  -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lmkl_mc3 -lmkl_def -lmkl_lapack95_ilp64 -L$(CFITSIO)/lib  -lcfitsio
LIBS =  -L$(LAPACK)  -llapack -lblas -L$(CFITSIO)/lib  -lcfitsio -L$(LGCC) -lgfortran 

PROGRAMS = test 

all: $(PROGRAMS) $(WMAPLIB)  

check: test
	./test

$(WMAPLIB): $(OBJS)
	ar r $@ $^

%: $(OBJS) %.o
	$(F90) $(FFLAGS) -o $@ $^ $(LIBS)

%.o: %.f90
	$(F90) $(FFLAGS) $(INCS) -c -o $@ $<

%.o: %.F90
	$(F90) $(FFLAGS) $(INCS) -c -o $@ $<

clean:
	$(RM) *.o  *.log *~  *.mod 

distclean: clean
	$(RM) *.a $(PROGRAMS)




