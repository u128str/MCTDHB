cmp=gnu

BLAS_LAPACK=gnu
DIR_BLAS=/home/alexej/lib
DIR_BLAS=/usr/lib

FFT=FFTW

#mpi_f90=gfortran -I/usr/lib/openmpi/include -pthread -I/usr/lib/openmpi/lib -L/usr/lib/openmpi/lib -lmpi_f90 -lmpi_f77 -lmpi -O2 -fopenmp 
mpi_f90=gfortran -I/usr/lib/openmpi/include -pthread -I/usr/lib/openmpi/lib -L/usr/lib/openmpi/lib -lmpi_f90 -lmpi_f77 -lmpi -g -fopenmp 



#mpi_f90=gfortran -I/usr/lib/openmpi/include -pthread -I/usr/lib/openmpi/lib -L/usr/lib/openmpi/lib -lmpi -lmpi_f77 -lopen-rte -lopen-pal -O2 -fopenmp 
#mpi_f90=gfortran -I/usr/lib/openmpi/include -pthread  -L/home/alexej/lib -lmpi -lmpi_f77 -lopen-rte -lopen-pal -O2 -fopenmp -static

FFTW_INC_DIR=/usr/include/
FFTW_LIB_DIR=/usr/lib
#FFTW_LIB_DIR=/home/alexej/lib
