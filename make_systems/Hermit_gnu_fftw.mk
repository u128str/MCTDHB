#MDL=module load compiler/gnu/4.5 mpi/openmpi/1.4.3-gnu-4.5 numlib/lapack/3.3.1 numlib/fftw/3.3.2-openmpi-1.4.3-gnu-4.5


cmp=gnu

BLAS_LAPACK=gnu
DIR_ACML=$(ACML_DIR)/gfortran64/lib/

FFT=FFTW

mpi_f90=ftn -O2  -fopenmp

FFTW_INC_DIR=$(FFTW_INC)
FFTW_LIB_DIR=$(FFTW_DIR)
