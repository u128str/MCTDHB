cmp=ifort

BLAS_LAPACK=ifort

FFT=MKLFFT

mpi_f90=mpiifort -O2   -mkl  -qopenmp  
#mpi_f90=mpiifort -O2  -pthread  -mkl  -static -fPIE -pie -shared
#mpi_f90=mpiifort -O2 -lc -lpthread  -mkl  -static -static-libgcc  -fPIE -pie -fopenmp
#mpi_f90=mpiifort -O2  -static-intel -static_mpi  -static-libgcc -static-libstdc++ -static
#mpi_f90=mpiifort -O2 -static -static-intel -static_mpi -static-libgcc -static-libstdc++ 

#mpi_f90=mpiifort -O2   -debug none -debug noinline-debug-info  -qopenmp      -static-intel -static_mpi -static-libgcc -static-libstdc++ 
# -L/cvmfs/hybrilit.jinr.ru/sw/intel/2016.1.150/compilers_and_libraries_2016.1.150/linux/mpi/intel64/lib/libmpifort.a
#mpi_f90=mpiifort -O2   -mkl  

#FCFLAGS = -cpp -static-intel -static-libgcc
#LDFLAGS = -pie -static-intel -static-libgcc -static-libstdc++ -static
LDFLAGS = --static
LDFLAGS =

DIR_MKL=$(MKLROOT)/lib/intel64/
