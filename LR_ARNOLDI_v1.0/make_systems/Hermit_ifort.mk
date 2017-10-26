cmp=ifort

BLAS_LAPACK=ifort

FFT=MKLFFT

mpi_f90=ftn -g -mkl=parallel  -O0 

DIR_MKL=$(MKLROOT)/lib/intel64/
