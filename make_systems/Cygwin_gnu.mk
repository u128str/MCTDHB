cmp=gnu

BLAS_LAPACK=gnu

FFT=FFTW

#mpi_f90=mpif90  -fopenmp -O3                   -g  -fbounds-check -fopenmp -Wall -pedantic -ansi
mpi_f90=mpifort  -fopenmp -O3  #-static 
#DIR_BLAS= "c:/cygwin64/usr/lib/"
DIR_BLAS= /usr/lib/

#FFTW_INC_DIR= "C:/cygwin64/usr/include/"
FFTW_INC_DIR= /usr/include/
#FFTW_LIB_DIR= "C:/cygwin64/usr/lib/"
FFTW_LIB_DIR= /usr/lib/

platform=windows
$(info  in Cygwun FFT_INC=$(FFTW_INC__DIR))
