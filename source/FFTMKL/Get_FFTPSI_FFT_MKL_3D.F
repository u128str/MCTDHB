! FFT business with MKL
!3D problem
! NdimMax = 3D(NDX*NDY*NDZ)
      SUBROUTINE Get_FFTPSI_FFT_MKL_3D(psiin,fftpsi,modus)
      USE SHARED_DIMS
!      USE DVR_ALL
!      USE W_INTERPARTICLE
      Use MKL_DFTI
!      USE omp_lib
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER ::  ierr
      COMPLEX*16,   DIMENSION(NDX*NDY*NDZ) :: fftpsi,psiin
c==================== Local Arrays for 3D Transformations 
      INTEGER :: I,J,K,I1,incr,j1,X_i,ind,modus
      REAL*8 :: tm
cccccccc MKL STUFF
      INTEGER :: INFO
      INTEGER :: L(3)
      type(DFTI_DESCRIPTOR), POINTER :: Desc_Handle_Dim1
      Integer :: Status
!      COMPLEX*16 :: ortky(NDY),ortkX(NDX),ortkz(NDZ)
      tm=MPI_WTIME(ierr)
c===============================================================================
! True  3D FFT
c#########################################################################
             fftpsi=psiin
c#########################################################################
      L(1)=NDX
      L(2)=NDY
      L(3)=NDZ
      Status = DftiCreateDescriptor(Desc_Handle_Dim1,DFTI_DOUBLE,
     &  DFTI_COMPLEX, 3, L )
      Status = DftiCommitDescriptor(Desc_Handle_Dim1) 

             mds: Select case (modus)
        Case (1)
        Status = DftiComputeForward(Desc_Handle_Dim1,fftpsi)
        Case (2)
        Status = DftiComputeBackward(Desc_Handle_Dim1,fftpsi)
              end  select mds

      Status = DftiFreeDescriptor(Desc_Handle_Dim1)
      timeCPU(3)=timeCPU(3)+MPI_WTIME(ierr)-tm
c================================================================================
      END  SUBROUTINE Get_FFTPSI_FFT_MKL_3D
