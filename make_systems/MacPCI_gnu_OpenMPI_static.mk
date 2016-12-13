cmp=gnu

BLAS_LAPACK=gnu

FFT=FFTW

mpi_f90=/Users/tc-user/Alexej/openmpi-1.8.4_static/bin/mpif90  -O2 -fopenmp  -Bstatic 
mpi_f90=/Users/tc-user/Alexej/openmpi-1.8.4_new/bin/mpif90 -O2 -fopenmp  -static-libgfortran -static-libgcc -lquadmath 

mpi_f90=/Users/tc-user/Alexej/openmpi-1.8.4_full_static/bin/mpif90 -O2 -fopenmp  -static-libgfortran -lquadmath -static-libgcc -L /Users/tc-user/Alexej/openmpi-1.8.4_full_static/lib_static -lmpi -lmpi_mpifh -lopen-rte -lopen-pal

FFTW_INC_DIR=/usr/local/include/
FFTW_LIB_DIR=/usr/local/lib


