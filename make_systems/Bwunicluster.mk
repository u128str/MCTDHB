cmp=ifort

BLAS_LAPACK=ifort

FFT=MKLFFT

mpi_f90=mpiifort -mkl=parallel  -O2

DIR_MKL=$(MKLROOT)/lib/intel64/
