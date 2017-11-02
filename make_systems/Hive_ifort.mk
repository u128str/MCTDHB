#ifort on Hive : module add intel/2017.02   
cmp=ifort

ARNOLDI=PARPACK
DIR_PARPACK=/data/home/alon/shaldar/arpack-ng-master/build/lib
PARPACK=-L $(DIR_PARPACK) -lparpack -lblas -llapack -lparpack

BLAS_LAPACK=ifort

FFT=MKLFFT

mpi_f90=mpiifort  -mkl  -O2 -qopenmp
mpi_f90=mpiifort  -mkl  -qopenmp

DIR_MKL=$(MKLROOT)/lib/intel64/
