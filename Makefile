# USAGE: make mk_file=MinGW_Win_fftw.mk
platform=unix
MCTDHB_V="V3.3.01" # ATTN!! Also correct files: source/ParaMain.F  and source/MODULES_ALL_allocate.F90
#=================== Selector for already tuned configurations  PCI, K100, BwGrid, etc
#mk=./make_systems/BwGriD_ifort.mk
#mk=./make_systems/BwGriD_gnu_fftw.mk

#mk=./make_systems/PCI_gnu_ifort.mk #cmp=gnu FFTW is a wrapper of  MKL
#mk=./make_systems/PCI_ifort.mk
#mk=./make_systems/PCI_ifort_MPICH_static.mk #STATIC version for UNIX should work without any external libs
#mk=./make_systems/PCI_ifort_MPICH.mk
#mk=./make_systems/PCI_gnu.mk
#mk=./make_systems/PCI_gnu_static.mk
#mk=./make_systems/PCI_gnu_MPICH.mk


#mk=./make_systems/MinGW_Win_fftw.mk
#mk=./make_systems/Cygwin_gnu.mk

#mk=./make_systems/Hermit_ifort.mk
#mk=./make_systems/Hermit_gnu_fftw.mk

#mk=./make_systems/Hydra_ifort_static.mk # On Hydra creation of the semi-static  boinaries with ifortmpi mkl etc..

mk=./make_systems/Ubuntu_gnu.mk
#mk=./make_systems/SUSE_gnu_fftw.mk #standard for SUSE

#mk=./make_systems/MacPCI_gnu_OpenMPI_static.mk #STATIC version for UNIX should work without any external libs
#=================== Selector for Compiler 

ifneq ($(mk_file),)
ifneq ($(wildcard ./make_systems/$(mk_file)),)
mk=./make_systems/$(mk_file)
$(info  Make222 file used: ($(wildcard $(mk))))
endif
endif

include $(mk)

$(info  Make file used: ($(wildcard $(mk))))


ifeq ($(cmp),)
cmp=gnu
endif 
$(info  cmp=$(cmp))
#==================  Selector for BLAS LAPACK ====================================
#================= If you use this option, you have to provide a directory to the libraries 
#BLAS_LAPACK=intel  # change  DIR_MKL=$(MKLROOT)/lib/ia32/ 
ifeq ($(BLAS_LAPACK),'')
BLAS_LAPACK=gnu
endif 
$(info  BLAS_LAPACK=$(BLAS_LAPACK))
#BLAS_LAPACK=gnu   # change DIR_BLAS=/usr/lib  # type: locate blas.a
#==========================================================================================

#==================  Selector for FFT  ====================================
#FFT=FFTW   #type:  locate  fftw3.f
#FFT=MKLFFT
#FFT=FFTW_MKL
#FFT=CUDACPP
#FFT=CUDAPGI
ifeq ($(FFT),'')
FFT=FFTW
endif 
$(info  FFT=$(FFT))

FF=mpif90
ifneq ($(mpi_f90),'')
FF=$(mpi_f90)         #ifort
endif 
$(info  FF=$(FF))
#FF=mpif90.openmpi #gfortran on PCI
FC =$(FF)  -qopenmp
FC =$(FF)  -openmp
FC =$(FF)  
#FC =$(FF) -O1  -fopenmp  -openmp #gun
$(info  FC=$(FC))
#======================================== My user_guesslib Dynamical library ========================
MY_LIB= -Bdynamic -L./ -lguess
#MY_LIB=  -Bdynamic libguess.so
F90_LD= $(FC) -fPIC -shared  -Bdynamic  

#MY_LIB= -Bdynamic libguess.so
#F90_LD= ftn -fPIC -shared 
#F90_LD= ftn  

#==========================================================================================

#================= INTEL ===================================================================#
ifeq ($(BLAS_LAPACK),intel)
ifeq ($(DIR_MKL),'')
DIR_MKL=$(MKLROOT)/lib/ia32/ 
endif
BLAS=  -Wl,--start-group -L$(DIR_MKL) -lmkl_intel -lmkl_sequential -lmkl_core -Wl,--end-group -lpthread -I$(MKL_INC_DIR)
BLAS=
endif
#================= GCC GNU ===================================================================#
ifeq ($(BLAS_LAPACK),gnu)
ifeq ($(DIR_BLAS),)
DIR_BLAS=/usr/lib  
DIR_BLAS=/Users/tc-user/Alexej/lib 
endif
BLAS =  -L $(DIR_BLAS)  -lblas -llapack 
ifneq ($(DIR_ACML),) 
DIR_BLAS=$(DIR_ACML)
BLAS =  -L $(DIR_BLAS) -lacml
$(info  ACML is used from=$(DIR_BLAS))
endif
endif
$(info  BLAS libs are from=$(DIR_BLAS))
#==========================================================================================
#======================================================= FFT inc and compilation keys ==========================
# link to FFT libraries via MKL
ifeq ($(FFT),MKLFFT) 
MKL_INC_DIR=$(MKLROOT)/include/
INC_ALL =   -I$(MKL_INC_DIR)  
LIBS =    -mkl
MKL_DFTI_MOD=mkl_dfti.mod 
endif

# link to FFT libraries
ifeq ($(FFT),FFTW)
ifeq ($(FFTW_INC_DIR),)
FFTW_INC_DIR=/usr/include/
endif
$(info  FFTW_INC_DIR=$(FFTW_INC_DIR))
ifeq ($(FFTW_LIB_DIR),)
FFTW_LIB_DIR=/usr/lib
endif
$(info  FFTW_INC_DIR=$(FFTW_LIB_DIR))
LIBS=  -L $(FFTW_LIB_DIR) -lfftw3  -lfftw3_threads 
ifeq ($(platform),windows)
LIBS=  -L $(FFTW_LIB_DIR) -lfftw3-3 -L $(MPI_LIB_DIR)
endif
INC_ALL =  -I $(FFTW_INC_DIR) 
endif

# link to MKL_FFT libraries
ifeq ($(FFT),FFTW_MKL) 
MKL_INC_DIR=$(MKLROOT)/include/
INC_DIR_FFT=$(MKL_INC_DIR)/fftw
INC_ALL =   -I $(INC_DIR_FFT) 
LIBS =    -lfftw3_threads -lfftw3 -mkl
MKL_DFTI_MOD=mkl_dfti.mod
endif

# link to CUDA_FFT libraries with CPP
ifeq ($(FFT),CUDACPP) 
CUDALIB := /common/cuda
LIBS = -L $(CUDALIB)/lib64 -lcufft -lcudart  -lcublas
endif
#====================================== Hard-core Co-Processor variables defining FFT
ifeq ($(FFT),FFTW_MKL) #FFTW but with MKL wrappers...
FFTFLAG= -DSFX1D=FFT_FFTW_1D -DSFX2D=FFT_FFTW_2D  -DSFX3D=FFT_FFTW_3D -DFFTIMEST=FFTWIMEST \
         -DFFTW=.TRUE.  -DCUDACPP=.FALSE.  -DMKLFFT=.FALSE.  -DCUDAPGI=.FALSE.
DIRSFFT := ./source/FFTFFTW
endif
ifeq ($(FFT),FFTW)     #Pure FFTW 
FFTFLAG= -DSFX1D=FFT_FFTW_1D -DSFX2D=FFT_FFTW_2D  -DSFX3D=FFT_FFTW_3D -DFFTIMEST=FFTWIMEST \
         -DFFTW=.TRUE.  -DCUDACPP=.FALSE.  -DMKLFFT=.FALSE.  -DCUDAPGI=.FALSE.
DIRSFFT := ./source/FFTFFTW
endif
ifeq ($(FFT),CUDACPP) 
FFTFLAG= -DSFX1D=FFT_CUDACPP_1D -DSFX2D=FFT_CUDACPP_2D  -DSFX3D=FFT_CUDACPP_3D  -DFFTIMEST=CUDACPPIMEST\
       -DCUDACPP=.TRUE.  -DFFTW=.FALSE.  -DMKLFFT=.FALSE.  -DCUDAPGI=.FALSE.
DIRSFFT := ./source/FFTCUDACPP
endif
ifeq ($(FFT),MKLFFT) 
FFTFLAG= -DSFX1D=FFT_MKL_1D -DSFX2D=FFT_MKL_2D  -DSFX3D=FFT_MKL_3D -DFFTIMEST=MKLIMEST \
         -DMKLFFT=.TRUE.  -DFFTW=.FALSE.  -DCUDACPP=.FALSE.   -DCUDAPGI=.FALSE.
DIRSFFT := ./source/FFTMKL
LIBRARIES = 
endif
ifeq ($(FFT),CUDAPGI) 
FFTFLAG= -DSFX1D=FFT_CUDAPGI_1D -DSFX2D=FFT_CUDAPGI_2D  -DSFX3D=FFT_CUDAPGI_3D   -DFFTIMEST=CUDAPGIIMEST \
       -DCUDAPGI=.TRUE. -DCUDACPP=.FALSE.  -DFFTW=.FALSE.  -DMKLFFT=.FALSE. 
DIRSFFT := ./source/FFTCUDAPGI
endif
#MCTDHB_V="V3.2.27"
#VERFLAG= -DMCTDHB_V=$(MCTDHB_V)
VERFLAG= 

$(info  FFTFLAG=$(FFTFLAG))
#============================ Check For errors ========================================================== 
#check::
#	  @$(err) $(ERR)
#============================ Here one decides on single exe or with dynamical library ..====================
#with libguess.so
comp_mode=Split
#normal mode -- single exe file
comp_mode=
$(info  comp_mode0=$(comp_mode))
#============================ Here list of files to compile and rules for executables====================
include ./make_systems/mctdhb.mk

$(info  CURDIR=$(CURDIR))
