#MDL=module load compiler/gnu/4.5 mpi/openmpi/1.4.3-gnu-4.5 numlib/lapack/3.3.1 numlib/fftw/3.3.2-openmpi-1.4.3-gnu-4.5
#-static -lgomp -lpthread -static-libgfortran -static-libgcc -static-libstdc++  

cmp=gnu

BLAS_LAPACK=gnu
#DIR_BLAS=c:/BLAS_LAPACK64/
DIR_BLAS=c:/BLAS_LAPACK

FFT=FFTW

mpi_f90=c:/MinGW/bin/gfortran.exe -O2 -fcheck=all -fopenmp -DWIN32  -static -I"c:\IBM\Platform-MPI\include\32" 
mpi_f90=c:/MinGW/bin/gfortran.exe -O2 -fopenmp -DWIN32  -I"c:\IBM\Platform-MPI\include\32" 
#mpi_f90=c:/MinGW/bin/gfortran.exe -g -Wl,--stack,167772161 -fcheck=all -fopenmp -DWIN32  -I"c:\IBM\Platform-MPI\include\32" 
#mpi_f90=c:/TDM-GCC-64/bin/gfortran.exe -fopenmp  -I./ -I"c:\IBM\Platform-MPI\include\64" 
	
#MPI_LIB_DIR="c:\IBM\Platform-MPI\bin" -lpcmpi64 
MPI_LIB_DIR="c:\IBM\Platform-MPI\bin" -lpcmpi32

FFTW_INC_DIR=c:/fftw64
FFTW_LIB_DIR=c:/fftw64

FFTW_INC_DIR=c:/fftw
FFTW_LIB_DIR=c:/fftw



platform=windows
