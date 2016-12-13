! FFT business with FFT
!  modus=1 FFT
!  modus=2 back FFT
      SUBROUTINE Get_FFTPSI_FFT_FFTW_1D(psiin,fftpsi,modus)
      USE SHARED_DIMS
      USE DVR_ALL
      USE W_INTERPARTICLE
!      USE omp_lib
      IMPLICIT NONE
      INCLUDE 'mpif.h'
cccccccc FFTW stuff
      include 'fftw3.f'
      integer*8 plan,plany !attention to int*8 !!!!!!!
      integer iret,nthr
      COMPLEX*16 :: dummyX(NDX*NDY*NDZ)
      INTEGER ::  ierr,modus
      INTEGER :: FFTW_USE_WISDOM,isuccess
      COMPLEX*16,   DIMENSION(NDX*NDY*NDZ) :: fftpsi,psiin
      REAL*8 :: tm
      tm=MPI_WTIME(ierr)
c===============================================================================
! 1D FFT in X-dimension 
            dummyX=psiin
             mds: Select case (modus)
        Case (1)
      call dfftw_plan_dft_1d(plan,NDX*NDY*NDZ,dummyX,dummyX,
     .                       FFTW_FORWARD,
     .                       FFTW_ESTIMATE)
        Case (2)
      call dfftw_plan_dft_1d(plan,NDX,dummyX,dummyX,
     .                       FFTW_BACKWARD,
     .                       FFTW_ESTIMATE)
              end  select mds
      call dfftw_execute_dft(plan,dummyx,dummyx)
         fftpsi=dummyX
      call dfftw_destroy_plan(plan)   
c================================================================================
      timeCPU(3)=timeCPU(3)+MPI_WTIME(ierr)-tm
c================================================================================
      END  SUBROUTINE Get_FFTPSI_FFT_FFTW_1D
