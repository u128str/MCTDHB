cmp=gnu

BLAS_LAPACK=ifort

#FFT=FFTW_MKL

FFT=MKLFFT

mpi_f90=mpif90 -m64 -I${MKLROOT}/include -O2 

#DIR_MKL=$(MKLROOT)/lib/ia32/

LIBS = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a -Wl,--end-group -liomp5 -lpthread -lm -ldl

DIR_MKL=$(MKLROOT)/lib/intel64/
