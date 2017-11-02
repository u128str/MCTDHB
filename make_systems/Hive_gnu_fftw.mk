# module load FFTW/3.3.6-gompi-2017a OpenBLAS/0.2.19-GCC-6.3.0-2.27-LAPACK-3.7.0
cmp=gnu

BLAS_LAPACK=gnu

DIR_BLAS=/data/apps/Easybuild/apps/OpenBLAS/0.2.19-GCC-6.3.0-2.27-LAPACK-3.7.0/lib 
FFT=FFTW

mpi_f90=srun mpif90 -O4  -fopenmp 
#mpi_f90=srun mpif90 -O2  -fopenmp  -L /data/apps/Easybuild/apps/OpenBLAS/0.2.15-GCC-4.9.3-2.25-LAPACK-3.6.0/lib/ -lopenblas 

FFTW_INC_DIR=/data/apps/Easybuild/apps/FFTW/3.3.6-gompi-2017a/include
FFTW_LIB_DIR=/data/apps/Easybuild/apps/FFTW/3.3.6-gompi-2017a/lib
