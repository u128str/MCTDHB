cmp=ifort
BLAS_LAPACK=intel

FFT=MKLFFT

ARNOLDI=PARPACK
DIR_PARPACK=/home/u128str/ARPACK-NG/lib/
PARPACK=-L $(DIR_PARPACK) -lparpack #-lblas -llapack

mpi_f90=mpif90  -O2 -mkl
mpi_f90=mpif90  -O2  -fopenmp
mpi_f90=mpif90  -g  -fopenmp

#DIR_MKL=$(MKLROOT)/lib/ia32/
DIR_MKL=$(MKLROOT)/lib/intel64/
DIR_MKL=/opt/intel/compilers_and_libraries_2018.0.128/linux/mkl/lib/intel64
