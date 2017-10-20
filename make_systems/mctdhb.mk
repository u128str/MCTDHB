#comp_mode='Split1'
#comp_mode=''

PATH_OBJ=
MORB="files are taken from"
PRT="`pwd`"

MODF=  ./source/MODULES_ALL_allocate.F90 ./source/MODULES_SUBR.f90 ./source/MODULES_CI_SUBR.f90
#DIRS := ./source $(DIRSFFT) ./user_guesslib ./external
DIRS := ./source $(DIRSFFT) ./external

#ifeq ($(comp_mode),Split)
#LIBS+= $(MY_LIB)
#$(info  LYBS=$(LIBS))
#endif
ifeq ($(comp_mode),)
MY_LIB=
DIRS+= ./user_guesslib 
$(info  DIRS=$(DIRS))
endif

#$(info  DIRS=$(DIRS))
#$(info  LYBS=$(LIBS))


FILES_f    := $(foreach dir, $(DIRS), $(wildcard $(dir)/*.f))
FILES_FALL := $(foreach dir, $(DIRS), $(wildcard $(dir)/*.F))
FILES_f90ALL  := $(foreach dir, $(DIRS), $(wildcard $(dir)/*.f90))
FILES_F90  := $(foreach dir, $(DIRS), $(wildcard $(dir)/*.F90))


ifeq ($(platform),unix)
FILES_F=$(filter-out ./source/Func_windows.F ./source/Assist_Orb_windows.F ,$(FILES_FALL))
FILES_f90=$(filter-out ./source/PrdCIJKL1body_M_OMP_windows.f90 ./source/PrdCIJKL2body_M_OMP_windows.f90 ,$(FILES_f90ALL))
endif

ifeq ($(platform),windows)
FILES_F=$(filter-out ./source/Func_unix.F ./source/Assist_Orb_unix.F ,$(FILES_FALL))
FILES_f90=$(filter-out ./source/PrdCIJKL1body_M_OMP_unix.f90 ./source/PrdCIJKL2body_M_OMP_unix.f90 ,$(FILES_f90ALL))
endif
#$(info FFT=$(FILES_f90_ALL))
#$(info FFT=$(FILES_f90))

FILES_cu := $(foreach dir, $(DIRS), $(wildcard $(dir)/*.cu))
FILES_CUF := $(foreach dir, $(DIRS), $(wildcard $(dir)/*.CUF))


OBJ_f   := $(patsubst %.f, %.o, $(FILES_f))
OBJ_F   := $(patsubst %.F, %.o, $(FILES_F))
OBJ_f90 := $(patsubst %.f90, %.o, $(FILES_f90))
OBJ_F90 := $(patsubst %.F90, %.o, $(FILES_F90))
OBJ_cu := $(patsubst %.cu, %.o, $(FILES_cu))
OBJ_CUF := $(patsubst %.CUF, %.o, $(FILES_CUF))

FC+= $(VERFLAG)

#MODS=  ./source/MODULES_ALL_allocate.o ./source/MODULES_SUBR.o ./source/MODULES_CI_SUBR.o
MODS= $(MKL_DFTI_MOD) \
                      ci_all.mod ci_prod.mod dvr_all.mod parallel_ci.mod parallel_orb.mod rr_ww.mod shared_dms.mod  w_interpartile.mod\
                    pass_arg.mod \
                     ci_subr.mod \
                     interpreter.mod precision.mod

TIME = "`date +"%j_%s"`" 

.SUFFIXES:
.SUFFIXES: .CUF .cu .F .f .f90 .F90


%.o : %.F90 
	$(FC) $(FFTFLAG) $(INC_ALL)  -c $< -o $@
#	$(FC) -c $< -o ./objects/$(notdir $@)


%.o : %.f90 
	$(FC) $(FFTFLAG) $(INC_ALL)  -c $< -o $@
#	$(FC) -c $< -o ./objects/$(notdir $@)


%.o : %.F  
	$(FC) $(FFTFLAG) $(INC_ALL) -c $< -o $@
#	$(FF) -c $< -o ./objects/$(notdir $@)

%.o : %.f
	$(FC) $(FFTFLAG) $(INC_ALL) -c $< -o $@
#	$(FF) -c $< -o ./objects/$(notdir $@)

%.o : %.cu
	$(CUDALIB)/bin/nvcc -O1 -arch=compute_20 -c $(FFTFLAG) $(INC_ALL) $< -o $@

%.o : %.CUF 
	mpif90 -O2  -fast  -Mcuda=3.2 -c $(INC_ALL)  $< -o $@ 

add_dep=
RES=  $(MODS) ./bin/boson_MCTDHB_$(cmp)_$(FFT) ./bin/properties_LR_$(cmp)_$(FFT)  info 
ifeq ($(comp_mode),Split)
add_dep= libguess.a
endif

$(info  add_dep=$(add_dep))

all: $(add_dep) $(RES)

$(info  LIBS=$(LIBS))
$(info  MY_LIB=$(MY_LIB))

./bin/boson_MCTDHB_$(cmp)_$(FFT): $(OBJ_F) $(OBJ_F90) $(OBJ_f90) $(OBJ_f) $(OBJ_cu)  $(OBJ_CUF)  $(add_dep)
	$(FC)  -o $@ $(OBJ_F90) $(OBJ_f90) $(OBJ_F) $(OBJ_f) $(OBJ_cu) $(OBJ_CUF) $(LIBRARIES) $(BLAS) $(FFTFLAG) $(CUDA_ALL) $(INC_ALL) $(LIBS) $(MY_LIB)
#	cp ./libguess.so  ./bin/.
#	$(info  BINARY is constructed okay )

interpreter.mod precision.mod: dvr_all.mod
	$(FC) -c external/interpreter.f90  -o ./external/interpreter.o

mkl_dfti.mod: 
	$(FC) -c $(MKL_INC_DIR)/mkl_dfti.f90  -o ./source/mkl_dfti.o

pass_arg.mod: 
	$(FC) -c ./source/MODULES_SUBR.f90 -o ./source/MODULES_SUBR.o

ci_subr.mod: 
	$(FC) -c ./source/MODULES_CI_SUBR.f90 -o  ./source/MODULES_CI_SUBR.o
 
ci_all.mod ci_prod.mod dvr_all.mod parallel_ci.mod parallel_orb.mod rr_ww.mod shared_dms.mod w_interpartile.mod:  
	$(FC) -c ./source/MODULES_ALL_allocate.F90 -o  ./source/MODULES_ALL_allocate.o
 
clean:
	rm -r  *.a ./external/*.o ./external/*.mod ./source/*.o ./*.mod ./user_guesslib/libguess.so ./user_guesslib/*.o ./*.o ./bin/* ./source/FFTMKL/*.o ./source/FFTFFTW/*.o ./source/FFTCUDACPP/*.o  ./source/FFTCUDAPGI/*.o ./source/FFTCUDAPGI/*.mod  properties_LR/*.o properties_LR/*.mod ./libguess.so  

cl:
	rm -r  *.a ./external/*.o ./external/*.mod ./source/*.o ./*.mod ./user_guesslib/libguess.so ./user_guesslib/*.o ./*.o ./bin/* ./source/FFTMKL/*.o ./source/FFTFFTW/*.o ./source/FFTCUDACPP/*.o  ./source/FFTCUDAPGI/*.o ./source/FFTCUDAPGI/*.mod  properties_LR/*.o properties_LR/*.mod ./libguess.so 
clP:
	rm -r  *.a ./properties_LR/*.o ./properties_LR/*.mod

#================= Dynamic Library =================
DYN_DIRS := ./user_guesslib #All files from this directory are compiled and linked to the dynamical library 
DYN_FILES_F := $(foreach dir, $(DYN_DIRS), $(wildcard $(dir)/*.F))
DYN_OBJ_F   := $(patsubst %.F, %.o, $(DYN_FILES_F))

libguess.a: $(MODS) $(DYN_OBJ_F)  ./external/interpreter.o ./external/read_string.o 
	$(info  $(DYN_OBJ_F))
	$(F90_LD) $(DYN_FILES_F) -o  ./user_guesslib/libguess.a
	cp ./user_guesslib/libguess.a  ./bin/.
	cp ./user_guesslib/libguess.a  .

#================= Properties  =================
include ./make_systems/properties.mk

.PHONY: info
info:
	$(info MCTDHB package has been compiled with cmp=$(cmp))
	$(info BLAS and LAPACK=$(BLAS_LAPACK))
	$(info FFT=$(FFT))
	$(info with FFTFLAG=$(FFTFLAG))
	$(info  ======================================================================================)
	$(info  To use scripts do not forget to put in your .bashrc the following line)
	$(info  export mctdhb_dir=$(CURDIR))
	$(info  ======================================================================================)
	$(info  ========================= BBB: Be superB with the mctdhB =============================)
	$(info all executables are in )
#echo "`ls -ltr ./bin/*`"
