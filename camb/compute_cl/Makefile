SYSTYPE="IUCAA_PERSEUS"
#SYSTYPE="IUCAA_CRAY"

F90C=ifort
FFLAGS=-openmp -O2 -ip -W0 -WB -fpp2 -vec_report0 -O2 -fpic -DOPTIMIZE  

ifeq ($(SYSTYPE),"IUCAA_PERSEUS")
GCC=/usr/lib/gcc/x86_64-redhat-linux/4.1.2
CFITSIO=/data1/pdf/cjayanti/Software/cfitsio
LPACK=/data1/pdf/cjayanti/Software/lapack-3.4.2
MPI=/opt/software/intel/impi/4.1.0.024/intel64
FFLAGS+=-I$(MPI)/include
LDFLAGS=-L$(MPI)/lib -lmpi -lmpiif 
CAMB=../../camb
WMAP9=/data1/pdf/cjayanti/data/wmap_likelihood_v5/
CLIK=/data1/pdf/cjayanti/data/Planck2013/plc-1.0
endif

ifeq ($(SYSTYPE),"IUCAA_CRAY")
GCC=/usr/lib/gcc/x86_64-redhat-linux/4.1.2
CFITSIO=/data1/software/cfitsio
LPACK=/data1/software/lapack-3.5.0
MPI=/opt/platform_mpi
FFLAG+=-I$(MPI)/include -L$(MPI)/lib/linux_amd64 -lmpi
CAMB=../camb
WMAP9=../wmap_likelihood_v5
endif

WMAPLIB = libwmap9.a
CAMBLIB = libcamb.a

FFLAGS  +=-I$(CAMB) -I$(WMAP9) -I$(MPI)/include -I$(CLIK)/include -I. 
LDFLAGS += -L.  -lcamb  -lwmap9  -L$(CLIK)/lib -lclik -lclik_f90
LDFLAGS +=-L$(CFITSIO)/lib -lcfitsio  -L$(LPACK) -llapack -lblas -L$(GCC) -lgfortran -L$(MPI)/lib/linux_amd64 -lmpi 

CAMBOBJS=$(CAMB)/constants.o $(CAMB)/utils.o  $(CAMB)/subroutines.o $(CAMB)/inifile.o $(CAMB)/power_tilt.o $(CAMB)/recfast.o $(CAMB)/reionization.o $(CAMB)/modules.o $(CAMB)/bessels.o $(CAMB)/equations.o $(CAMB)/halofit.o $(CAMB)/lensing.o $(CAMB)/SeparableBispectrum.o $(CAMB)/cmbmain.o $(CAMB)/camb.o

WMAPOBJS=$(WMAP9)/read_archive_map.o $(WMAP9)/read_fits.o $(WMAP9)/healpix_types.o $(WMAP9)/br_mod_dist.o $(WMAP9)/WMAP_9yr_options.o $(WMAP9)/WMAP_9yr_util.o $(WMAP9)/WMAP_9yr_gibbs.o $(WMAP9)/WMAP_9yr_tt_pixlike.o $(WMAP9)/WMAP_9yr_tt_beam_ptsrc_chisq.o $(WMAP9)/WMAP_9yr_teeebb_pixlike.o $(WMAP9)/WMAP_9yr_tetbeebbeb_pixlike.o $(WMAP9)/WMAP_9yr_likelihood.o

PSOOBJS =   defs.o compute_cl.o main.o

all : camb wmap  $(CAMBLIB)   $(WMAPLIB) $(PSOOBJS) cosmopso

camb:
	cd $(CAMB)     && $(MAKE) FFLAGS="$(FFLAGS)"
	
wmap:
	cd $(WMAP9)    && $(MAKE) FFLAGS="$(FFLAGS)" 

$(WMAPLIB): $(WMAPOBJS)
	ar r $@ $^

$(CAMBLIB): $(CAMBOBJS)
	ar r $@ $^
%.o : %.f90
	$(F90C)  -c  $(FFLAGS) $< -o $@


cosmopso: $(PSOOBJS) 
	$(F90C)  $(FFLAGS) $(PSOOBJS) -o cosmopso $(LDFLAGS)   
clean_camb:
	cd $(CAMB)     && $(MAKE) clean  
clean_wmap:
	cd $(WMAP9)    && $(MAKE) clean
clean:
	rm -rf *~ *.o *.out *.mod  *.err cosmopso gbest.dat pvals.dat pbest.dat conver.dat 
allclean:
	 rm -rf *~ *.o *.out  *.err  *.mod $(WMAP9)/*.o $(WMAP9)/*.mod $(WMAP9)/*.a $(CAMB)/*.o $(CAMB)/*.mod	


	
