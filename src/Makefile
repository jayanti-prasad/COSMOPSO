SYSTYPE="IUCAA_PERSEUS"
#SYSTYPE="IUCAA_CRAY"

F90C=ifort
FFLAGS=-openmp -O2 -ip -W0 -WB -fpp2 -vec_report0 -O2 -fpic -DOPTIMIZE  

ifeq ($(SYSTYPE),"IUCAA_PERSEUS")
GCC=/usr/lib/gcc/x86_64-redhat-linux/4.1.2
CFITSIO=/data1/sysadmin/cjayanti/Software/cfitsio
LPACK=/data1/sysadmin/cjayanti/Software/lapack-3.4.2
MPI=/opt/software/intel/impi/4.1.0.024/intel64
FFLAGS+=-I$(MPI)/include
LDFLAGS=-L$(MPI)/lib -lmpi -lmpiif 
CAMB=../camb
WMAP9=../wmap_likelihood_v5/
CLIK=../plc-2.0/
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

RCAMB=$(CAMB)

WMAPLIB = libwmap9.a
CAMBLIB = libcamb.a

FFLAGS  +=-I$(RCAMB) -I$(WMAP9) -I$(MPI)/include -I$(CLIK)/include -I. 
LDFLAGS += -L.  -lcamb  -lwmap9  -L$(CLIK)/lib -lclik -lclik_f90
LDFLAGS +=-L$(CFITSIO)/lib -lcfitsio  -L$(LPACK) -llapack -lblas -L$(GCC) -lgfortran -L$(MPI)/lib/linux_amd64 -lmpi 

CAMBOBJS=$(RCAMB)/constants.o $(RCAMB)/utils.o  $(RCAMB)/subroutines.o $(RCAMB)/inifile.o $(RCAMB)/power_tilt.o $(RCAMB)/recfast.o $(RCAMB)/reionization.o $(RCAMB)/modules.o $(RCAMB)/bessels.o $(RCAMB)/equations.o $(RCAMB)/halofit.o $(RCAMB)/lensing.o $(RCAMB)/SeparableBispectrum.o $(RCAMB)/cmbmain.o $(RCAMB)/camb.o

WMAPOBJS=$(WMAP9)/read_archive_map.o $(WMAP9)/read_fits.o $(WMAP9)/healpix_types.o $(WMAP9)/br_mod_dist.o $(WMAP9)/WMAP_9yr_options.o $(WMAP9)/WMAP_9yr_util.o $(WMAP9)/WMAP_9yr_gibbs.o $(WMAP9)/WMAP_9yr_tt_pixlike.o $(WMAP9)/WMAP_9yr_tt_beam_ptsrc_chisq.o $(WMAP9)/WMAP_9yr_teeebb_pixlike.o $(WMAP9)/WMAP_9yr_tetbeebbeb_pixlike.o $(WMAP9)/WMAP_9yr_likelihood.o

PSOOBJS =   defs.o compute_cl.o  compute_wmap9.o compute_planck.o compute_likelihood.o  uniform.o setpara.o  pso.o pso_func_mpi.o main.o

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


	
