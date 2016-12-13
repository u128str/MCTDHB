#MDL=module load compiler/gnu/4.5 mpi/openmpi/1.4.3-gnu-4.5 numlib/lapack/3.3.1 numlib/fftw/3.3.2-openmpi-1.4.3-gnu-4.5

cmp=gnu

BLAS_LAPACK=gnu

FFT=FFTW

mpi_f90=mpif90 -O4  -fopenmp 

FFTW_INC_DIR=/opt/bwgrid/numlib/fftw/3.3.2-openmpi-1.4.3-gnu-4.5/include
FFTW_LIB_DIR=/opt/bwgrid/numlib/fftw/3.3.2-openmpi-1.4.3-gnu-4.5/lib
