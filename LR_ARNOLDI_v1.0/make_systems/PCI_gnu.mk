cmp=gnu

BLAS_LAPACK=gnu

FFT=FFTW

mpi_f90=mpif90.openmpi  -O2 -fopenmp -static
mpi_f90=mpif90.openmpi  -O2 -fopenmp 
mpi_f90=mpif90.openmpi  -g -fopenmp  -fcheck=bounds

FFTW_INC_DIR=/usr/include/
FFTW_LIB_DIR=/usr/lib
