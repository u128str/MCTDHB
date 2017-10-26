P_MODS= ./source/MODULES_ALL_allocate.o  

ifeq ($(platform),windows)
FILES_ADD= ./source/Func_windows.o ./source/Assist_Orb_windows.o 
endif
ifeq ($(platform),unix)
FILES_ADD=./source/Func_unix.o ./source/Assist_Orb_unix.o
endif


P_DEPS= $(FILES_ADD)  \
       ./external/interpreter.o\
       ./external/read_string.o\
       ./external/sind_cosd_etc.o\
      ./source/BINARY_I_O.o \
     ./source/FuncLib.o \
     ./source/Read_Init.o \
     ./source/PRINT_CIc.o \
     ./source/PRINT_RHO.o \
     ./source/Init_CI.o  \
     ./source/Get_Ind_vs_ii.o \
     ./source/Init_h.o \
     ./source/xvlib.o \
     ./source/mmlib.o \
     ./source/zvode.o \
     ./source/xarg.o \
     ./source/enoi.o \
     ./source/schmidtortho.o \
     ./source/op1lib.o \
     ./source/dvrweights.o \
     ./source/init1.o \
     ./source/genphi1.o \
     ./source/get_mom_pos_fft.o \
     ./source/compvtilde.o \
     ./source/HPSI_v2.o \
     ./source/T_FFT_1D.o \
     ./source/GetCIJKL*.o \
     ./source/PrdCIJKL*.o \
     ./source/Parallel_MNGR_CI_Part_v2.o \
     ./source/Parallel_MNGR_Orbital_Part.o \
     ./source/Get_H_W.o \
     ./source/Get_r_R.o \
     ./source/Share_r_R.o \
     ./source/Share_H_W.o \
     ./source/OPERS.o \
     ./source/Get_WSL.o \
     ./source/Nadr_Rhoall.o \
     ./source/Get_WSL_omp.o \
     ./source/Get_Operators.o \
     ./source/Get_OPSI_WSL_disbalanced_OMP.o \
     ./source/Get_Op_PSI.o \
     ./source/computewsl.o \
     ./source/DIR.o \
     ./source/Get_FFTPSI.o\
     ./source/Shift_Zero_K_FFT.o \
     ./source/Guess_Diag.o  ./source/Guess_Read.o \
     ./source/Get_OPSI_WSL_balanced_OMP.o \
     ./source/Assistent_All.o  \
     ./source/Assist_CI_v1.o   \
     ./source/Get_d_PSI.o   \
     ./source/Get_Op_PSI_Gnrl.o

#================= Dynamic Library =================
P_DIRS := ./LR_ARNOLDI $(DIRSFFT)
#P_DIRS := $(DIRSFFT)
ifeq ($(comp_mode),Split)
#LIBS+=MY_LIB
else
P_DIRS+= ./user_guesslib 
endif
P_FILES_F := $(foreach dir, $(P_DIRS), $(wildcard $(dir)/*.F))
P_OBJ_F   := $(patsubst %.F, %.o, $(P_FILES_F))
P_FILES_F90 := $(foreach dir, $(P_DIRS), $(wildcard $(dir)/*.f90))
P_OBJ_F90   := $(patsubst %.f90, %.o, $(P_FILES_F90))

$(info  LR ARNOLDI Directories =$(P_DIRS))

MOD_PR= lr_arnoldi_mod.mod


./bin/lr_arnoldi_$(cmp)_$(FFT): $(MOD_PR) $(P_OBJ_F90) $(P_OBJ_F) 
	$(FC) -o $@ $(P_OBJ_F) $(P_MODS)  $(P_OBJ_F90) $(P_DEPS)  $(LIBRARIES) $(BLAS) $(LIBS) $(MY_LIB) $(PARPACK) $(MEMDEBUG)
#	mv properties_LR_$(cmp)_$(FFT)  ./bin/.

lr_arnoldi_mod.mod: 
	$(FC) -c ./LR_ARNOLDI/MODULE_lr_arnoldi.f90 -o  ./LR_ARNOLDI/MODULE_lr_arnoldi.o



