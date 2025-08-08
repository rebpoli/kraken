# bevo3.mak - Machine and compiler make include file for bevo3

#### reset SIZE so that the make would pick up error of SIZE preprocessor

SIZEIT= rm -f dum?; $(SIZE) < ech 2> dum1; cat dum1; grep 0 dum1 > dum2;

######################################################## COMPILERS and FLAGS
#### define compilers, linker etc.

CC       = mpicc
CPP      = mpicxx
FORT     = mpif90
FORT90   = mpif90
LINK     = mpif90 -L$(INTEL_HOME)/lib/intel64 -limf 

#HYPRE_DIR  = /home/bganis/hypre-2.11.0/src/hypre
HYPRE_DIR  = /home/bganis/hypre-2.11.0-intel15/src/hypre

#SUPERLU_CFLAGS = -DAdd_ -DUSE_VENDOR_BLAS -I/home/bganis/SuperLU_DIST_3.3/SRC
#SUPERLU_LFLAGS = /home/bganis/SuperLU_DIST_3.3/lib/libsuperlu_dist_3.3.a \
                 -L/home/bganis/parmetis-4.0.3/build/Linux-x86_64/libparmetis -lparmetis \
                 -L/home/bganis/parmetis-4.0.3/build/Linux-x86_64/libmetis -lmetis

VTKPATH = /usr/local/lib/vtk-5.4
VTKINCL = /usr/local/include/vtk-5.4
VTKLIBS = -lvtkIO -lvtkFiltering -lvtkzlib -lvtkCommon -lvtksys

#LFLAGS   = -Wl,-rpath,$(ICES_MKL_LIB) $(LIBS_HYPRE) $(DAGHMB_LFLAGS) \
#           $(SUPERLU_LFLAGS) $(TECLIB)
#FFLAGS   = -c -shared-intel -mcmodel=medium $(VIS_FFLAGS) $(VIS_PFLAGS)
#FFLAGS90 = 
#CFLAGS   = -c -shared-intel -mcmodel=medium -DSAMG_UNIX_LINUX -DSAMG_LCASE_USCORE \
#	    $(VIS_CFLAGS) $(SUPERLU_CFLAGS)
#CPPFLAGS = -c -shared-intel -mcmodel=medium -DSAMG_UNIX_LINUX -DSAMG_LCASE_USCORE

LFLAGS   = -lstdc++ $(MACE_LFLAGS) $(DAGHMB_LFLAGS) $(LIBS_HYPRE) \
           $(SUPERLU_LFLAGS) $(TECLIB) $(TRILINOS_LIBS) $(SAMG_LIBS)
FFLAGS   = -c $(MACE_FFLAGS) $(DAGHMB_FFLAGS) \
            -shared-intel -mcmodel=medium $(FDEFS)
CFLAGS   = -c $(SAMG_CFLAGS) $(MACE_CFLAGS) $(DAGHMB_CFLAGS) \
           $(DISCOVER_CFLAGS) $(SUPERLU_CFLAGS)
CPPFLAGS = -c $(SAMG_CFLAGS) $(MACE_CPPFLAGS) $(DAGHMB_CPPFLAGS) \
           $(DISCOVER_CPPFLAGS) $(TRILINOS_CPPFLAGS)

ifdef DEBUG_FLAGS
 FFLAGS   += -g -CB -traceback -fp-stack-check -O0 -fpe0
 FFLAGS90 += -g -CB -traceback -fp-stack-check -O0 -fpe0
 CFLAGS   += -g
 CPPFLAGS += -g
else
 FFLAGS   += -O3 -unroll-aggressive -traceback
 FFLAGS90 += -O3 -unroll-aggressive -traceback
 CFLAGS   += -O3
 CPPFLAGS += -O3
# FFLAGS   += -O2 -unroll-aggressive -traceback
# FFLAGS90 += -O2 -unroll-aggressive -traceback
# CFLAGS   += -O2
# CPPFLAGS += -O2
endif  

########################################################### BLAS/LAPACK
# use libblas.a etc. for static libs (faster) 
# or lblas, llapack for dynamically linked libs (smaller executable)
#

BLAS_LIB = -L$(ICES_MKL_LIB) -lmkl_intel_lp64 -lmkl_sequential -lmkl_core

OBJST = $(OBJS) cputime.o

############################################################
# define the default dependencies

.f.o:
	$(FORT) $(FFLAGS) $*.f

%.o:%.f90
	$(FORT90) -c $(FFLAGS90) -o $@ $<

.c.o:
	$(CC) $(CFLAGS) $*.c

.C.o:
	$(CPP) $(CPPFLAGS) $*.C

.cpp.o:
	$(CPP) $(CPPFLAGS)  $*.cpp


$(EXENAM): size logstamp $(OBJST) $(LIBS_TO_STAMP) 
	$(LINK) $(OBJST) -o $(EXENAM) $(LFLAGS) $(LIBS) $(BLAS_LIB)
	rm -f logstamp
	date > logstamp
	touch logstamp
	chmod g+w * 2>/dev/null || true

size:
	$(LINK) -o size ../size/size.f

logstamp: 
	rm -f logstamp 
	date > logstamp
 
clean:
	rm -f $(WORK)/*.f
	rm -f $(WORK)/*.c
	rm -f $(WORK)/*.C
	rm -f $(WORK)/*.cpp
	rm -f $(WORK)/*.hpp
	rm -f $(WORK)/*.h
	rm -f $(WORK)/*.o
	rm -f $(WORK)/*.i
	rm -fr $(WORK)/ii_files/
	rm -f $(WORK)/*.lst
	rm -f $(WORK)/ech
	rm -f $(WORK)/dum*
	rm -f $(WORK)/*.f90
	rm -f $(WORK)/*.mod
	rm -f size
	rm -f $(EXENAM)
