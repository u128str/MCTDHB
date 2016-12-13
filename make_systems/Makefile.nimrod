cmp=sun
cmp=intel
#cmp=gcc
#cmp=Pathscale

FFT=FFTW
FFT=MKLFFT

ifeq ($(cmp),gcc)
FF=mpif90
FC =$(FF) -O4 -fopenmp 
BLAS=-lmkl_gf_lp64 -lmkl_lapack -lmkl_core -lmkl_gnu_thread 
ifeq ($(FFT),FFTW)
LIBRARIES=-L./FFTW/lib -lfftw3_threads -lfftw3 -lm
FFTFLAG= -DFFTW=1 -DMKLFFT=0
DIRS1 :=
endif
# link to MKL_FFT libraries
ifeq ($(FFT),MKLFFT)
FFTFLAG= -DFFTW=0 -DMKLFFT=1
DIRS1 := ./MKL_FFT
LIBRARIES =
endif
MY_LIB= -L./user_guesslib -lguess
F90_LD= gfortran -fPIC -shared  
endif

ifeq ($(cmp),Pathscale)
FF=mpif90
FC =mpif90 -O1  -fopenmp -fixedform
FC =mpif90 -O1  -fopenmp 
BLAS =  -L/opt/acml/current/pathscale64/lib -lacml -lacml_mv
LIBRARIES =
MY_LIB= -L./user_guesslib -lguess
F90_LD= mpif90 -fPIC -shared 
endif

#================= INTEL ===================================================================#
ifeq ($(cmp),intel)
FF=ifort
# link to FFTW libraries
ifeq ($(FFT),FFTW)
#LIBRARIES=-L./FFTW/lib -lfftw3_threads -lfftw3 -lm 
LIBRARIES=-L./FFTW/lib -lfftw3 -lm 
FFTFLAG= -DFFTW=.TRUE. -DMKLFFT=.FALSE.
DIRS1 := 
endif
# link to MKL_FFT libraries
ifeq ($(FFT),MKLFFT) 
FFTFLAG= -DFFTW=.FALSE. -DMKLFFT=.TRUE. 
DIRS1 := 
LIBRARIES = 
endif
BLAS= -lmkl_lapack -lmkl_core -lmkl -lguide

BLAS= -L/opt/bwgrid/compiler/intel/ct_4.0/mkl/10.2.5.035/lib/em64t/ \
/opt/bwgrid/compiler/intel/ct_4.0/mkl/10.2.5.035/lib/em64t//libmkl_solver_lp64_sequential.a \
 -Wl,--start-group -lmkl_core  -lmkl_sequential -lmkl_intel_lp64  -Wl,--end-group -openmp -lpthread
BLAS= -L/opt/bwgrid/compiler/intel/ct_4.0/mkl/10.2.5.035/lib/em64t/ \
/opt/bwgrid/compiler/intel/ct_4.0/mkl/10.2.5.035/lib/em64t/libmkl_solver_lp64.a \
 -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group -openmp -lpthread

BLAS= -L$(MKLROOT) -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lmkl_intel_lp64  
BLAS= -L$(MKLROOT)  -lmkl_core -lmkl_intel_lp64 -lpthread -lmkl_intel_thread
FC =ifort  -O3 -parallel -recursive -openmp -par-schedule-runtime -fpp3
#FC =mpif90 -xSSE4.1  -O3 -parallel -openmp -par-schedule-runtime -fpp3
#FC =mpif90 -W4 -warn all  -O1 -g -parallel -openmp 
MPI= /usr/local/mpich2/
FF= ifort
LIBRARIES= -L$(LIB_MKL) -lmkl_lapack -lmkl -lguide -lpthread -lg2c  -L$(MPI)/lib -lfmpich -lmpich
LIBRARIES= -L$(MKLROOT)/lib  -lmkl_core -lmkl_intel_lp64 -lpthread -lmkl_intel_thread   -L$(MPI)/lib -lfmpich -lmpich
FC=$(FF) -O1 -openmp -traceback  -assume 2underscore   -I$(MPI)/include
FC=$(FF) -O1 -openmp  -parallel  -I$(MPI)/include
#FC =mpif90 -O1 -g -parallel -openmp 
MY_LIB= -L./user_guesslib -lguess
F90_LD= ifort -fPIC -shared  
endif
#=================SUN ======================================================================#
ifeq ($(cmp),sun)
`export LD_LIBRARY_PATH=/opt/SUNWhpc/HPC8.0/lib:$LD_LIBRARY_PATH`
MPI= 
FF=/opt/SUNWhpc/HPC8.0/bin/mpif90
FC=$(FF)  -O3 -xopenmp=parallel -I$(MPI)/include -xlic_lib=sunperf -L/opt/SUNWhpc/HPC8.0/lib 
LIBRARIES= 
MY_LIB= -L ./user_guesslib -lguess
F90_LD= f90 -pic -G
endif
#============================================================================================#

PATH_OBJ=
MORB="files are taken from"
PRT="`pwd`"

#VPATH = ./source: ./user_guesslib


DIRS := ./source
FILES_f := $(foreach dir, $(DIRS), $(wildcard $(dir)/*.f))
FILES_F := $(foreach dir, $(DIRS), $(wildcard $(dir)/*.F))
FILES_f90 := $(foreach dir, $(DIRS), $(wildcard $(dir)/*.f90))


OBJ_f   := $(patsubst %.f, %.o, $(FILES_f))
OBJ_F   := $(patsubst %.F, %.o, $(FILES_F))
OBJ_f90 := $(patsubst %.f90, %.o, $(FILES_f90))

MKL_DFTI_MOD=
ifeq ($(FFT),MKLFFT)
MKL_DFTI_MOD=./mkl_dfti.mod
endif

MODS= $(MKL_DFT_MOD) ./source/MODULES_ALL_allocate.o ./source/MODULES_SUBR.o ./source/MODULES_CI_SUBR.o

TIME = "`date +"%j_%s"`" 

.SUFFIXES:
.SUFFIXES:.F .f .f90 .o


%.o : %.f90
	$(FC) $(FFTFLAG)  -c $< -o $@
#	$(FC) -c $< -o ./objects/$(notdir $@)

%.o : %.F
	$(FC) $(FFTFLAG) -c $< -o $@
#	$(FF) -c $< -o ./objects/$(notdir $@)

%.o : %.f
	$(FC) $(FFTFLAG) -c $< -o $@
#	$(FF) -c $< -o ./objects/$(notdir $@)


boson_MCTDHB_$(cmp): $(MODS) mkl_dfti.mod $(OBJ_f90) $(OBJ_F) $(OBJ_f) libguess.so
	$(FC) -o $@ $(OBJ_f) $(OBJ_F)  $(OBJ_f90) $(BLAS) $(MY_LIB) $(LIBRARIES) $(FFTFLAG)
	mv boson_MCTDHB_$(cmp)  ./bin/.
	cp -f ./user_guesslib/libguess.so ./bin/.
#	mv boson_MCTDHB_$(cmp)  ./bin/boson_MCTDHB_$(cmp)_$(TIME)
	echo "`ls -ltr ./bin/*`"

mkl_dfti.mod: $(MKLROOT)/include/mkl_dfti.f90 
	cp -f  $(MKLROOT)/include/mkl_dfti.f90 ./.
	$(FC) $(FFTFLAG)  -c ./mkl_dfti.f90

clean:
	rm -r ./source/*.o ./*.mod ./user_guesslib/libguess.so ./user_guesslib/*.o
	rm -f  ./mkl_dfti* ./*.o

cl:
	rm -r ./source/*.o ./*.mod ./user_guesslib/libguess.so ./user_guesslib/*.o
	rm -f  ./mkl_dfti* ./*.o
#================= Dynamic Library =================
DYN_DIRS := ./user_guesslib #All files from this directory are compiled and linced to the dynamical library
DYN_FILES_F := $(foreach dir, $(DYN_DIRS), $(wildcard $(dir)/*.F))
DYN_OBJ_F   := $(patsubst %.F, %.o, $(DYN_FILES_F))

libguess.so : $(DYN_OBJ_F)
	$(F90_LD) $(DYN_FILES_F) -o ./user_guesslib/libguess.so
#	mv $(notdir $(DYN_OBJ_F)) $(DYN_DIRS)/.


#libguess.so : ./user_guesslib/Guess_PSI.F ./user_guesslib/Guess_CI.F ./user_guesslib/Get_InterParticle.F
#	f90 -pic -G Guess_PSI.F Guess_CI.F  Get_InterParticle.F VTRAP_EXT_TD.F PotfitAll.F   -o ./user_guesslib/libguess.so



