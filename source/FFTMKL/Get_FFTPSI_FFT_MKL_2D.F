! FFT business with MKL
! X_i index pointing to the T_x,T_y or T_z will be calculated - only for
! 2D and 3D problems
! NdimMax = 1D(NDX)
! NdimMax = 2D(NDX*NDY)
! NdimMax = 3D(NDX*NDY*NDZ)
      SUBROUTINE Get_FFTPSI_FFT_MKL_2D(psiin,fftpsi,modus) 
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

c  oksa add
      INTEGER ::  Ntmp
	  REAL*8 ::   tmoksa
cccccccc MKL STUFF
      INTEGER :: STRIDE(2),INFO
      INTEGER :: L(2)
      type(DFTI_DESCRIPTOR), POINTER :: Desc_Handle_Dim1
      Integer :: Status
      tm=MPI_WTIME(ierr)
  
   
!      Ntmp= omp_get_max_threads()
!        write(6,*) "I OKSA_MKL : max threads =!", Ntmp
c===============================================================================
! True  2D FFT
c#########################################################################
             fftpsi=psiin
c#########################################################################
      L(1)=NDX
      L(2)=NDY
      Status = DftiCreateDescriptor(Desc_Handle_Dim1, DFTI_DOUBLE,
     &  DFTI_COMPLEX, 2, L )
      Status = DftiCommitDescriptor( Desc_Handle_Dim1) 
            mds: Select case (modus)
        Case (1)
        Status = DftiComputeForward(Desc_Handle_Dim1,fftpsi)
        Case (2)
        Status = DftiComputeBackward(Desc_Handle_Dim1,fftpsi)
             end  select mds
      Status = DftiFreeDescriptor(Desc_Handle_Dim1)
      timeCPU(3)=timeCPU(3)+MPI_WTIME(ierr)-tm
c================================================================================
      END  SUBROUTINE Get_FFTPSI_FFT_MKL_2D
