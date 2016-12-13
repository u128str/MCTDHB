cmp=ifort

BLAS_LAPACK=ifort

FFT=MKLFFT

mpi_f90=ftn -mkl=parallel  -O2

DIR_MKL=$(MKLROOT)/lib/intel64/
