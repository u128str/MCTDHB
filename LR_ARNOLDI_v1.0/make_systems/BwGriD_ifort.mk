cmp=ifort

BLAS_LAPACK=ifort

FFT=MKLFFT

mpi_f90=mpif90 -O4   -mkl=parallel 

DIR_MKL=$(MKLROOT)/lib/intel64/
