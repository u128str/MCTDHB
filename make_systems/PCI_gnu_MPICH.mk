cmp=gnu

BLAS_LAPACK=gnu

FFT=FFTW

#mpi_f90=/home/alexej/mpich-3.1.2_install_GNU/bin/mpif90  -fopenmp -O1
mpi_f90=/home/alexej/mpich2-1.0.8p1/bin/mpif90  -fopenmp -O1
mpi_f90=/home/alexej/mpich2-1.0.8p1/bin/mpif90  -fopenmp -O1
mpi_f90=mpif90.openmpi -g  -O0 -fopenmp  -fcheck=bounds

FFTW_INC_DIR=/usr/include/
FFTW_LIB_DIR=/usr/lib
