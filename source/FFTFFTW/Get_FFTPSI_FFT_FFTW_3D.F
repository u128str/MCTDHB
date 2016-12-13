! FFT business with FFT
!  modus=1 FFT
!  modus=2 back FFT
      SUBROUTINE Get_FFTPSI_FFT_FFTW_3D(psiin,fftpsi,modus)
      USE SHARED_DIMS
      USE DVR_ALL
      USE W_INTERPARTICLE
      USE omp_lib
      IMPLICIT NONE
      INCLUDE 'mpif.h'
cccccccc FFTW stuff
      include 'fftw3.f'
      integer*8 plan,plany
      integer iret,nthr
      COMPLEX*16 :: dummyXYZ(NDX*NDY*NDZ)
      INTEGER ::  ierr
      INTEGER :: FFTW_USE_WISDOM,isuccess,modus
      COMPLEX*16,   DIMENSION(NDX*NDY*NDZ) :: fftpsi,psiin
      REAL*8 :: tm
      tm=MPI_WTIME(ierr)
c===============================================================================
! 3D FFT 
            dummyXYZ=psiin
       nthr = OMP_GET_MAX_THREADS()
       call dfftw_init_threads(iret)
       call dfftw_plan_with_nthreads(nthr)
             mds: Select case (modus)
        Case (1)
      call dfftw_plan_dft_3d(plan,NDX,NDY,NDZ,dummyXYZ,dummyXYZ,
     .                FFTW_FORWARD,FFTW_ESTIMATE,fftw_unaligned)
        Case (2)
      call dfftw_plan_dft_3d(plan,NDX,NDY,NDZ,dummyXYZ,dummyXYZ,
     .                FFTW_BACKWARD,FFTW_ESTIMATE,fftw_unaligned)
              end  select mds
      call dfftw_execute_dft(plan,dummyXYZ,dummyXYZ)
         fftpsi=dummyXYZ
      call dfftw_destroy_plan(plan)   
c================================================================================
c================================================================================
      timeCPU(3)=timeCPU(3)+MPI_WTIME(ierr)-tm
c================================================================================
      END  SUBROUTINE Get_FFTPSI_FFT_FFTW_3D
