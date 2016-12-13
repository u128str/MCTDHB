cmp=ifort

BLAS_LAPACK=ifort

FFT=MKLFFT

#mpi_f90=/home/alexej/mpich-3.1.2_install/bin/mpif90  -O2 -mkl -static
mpi_f90=/home/alexej/mpich-3.1.2_install/bin/mpif90   -O2 -mkl -static
#mpi_f90=/home/alexej/mpich-3.1.2_install/bin/mpif90   -g -mkl -static

DIR_MKL=$(MKLROOT)/lib/ia32/

