cmp=ifort

BLAS_LAPACK=ifort

FFT=MKLFFT

mpi_f90=ftn -O2 -openmp -mkl=parallel -mtune=core-avx2

DIR_MKL=$(MKLROOT)/lib/intel64/
