! FFT business with MKL 1D
!  modus=1 FFT
!  modus=2 back FFT
      SUBROUTINE Get_FFTPSI_FFT_MKL_1D(psiin,fftpsi,modus) 
      USE SHARED_DIMS
      USE DVR_ALL
      USE W_INTERPARTICLE
      Use MKL_DFTI
!      USE omp_lib
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER ::  ierr
      COMPLEX*16,   DIMENSION(NDX*NDY*NDZ) :: psitmp,fftpsi,psiin
c==================== Local Arrays for 3D Transformations 
      INTEGER :: m,I,J,K,I1,incr,j1,X_i,ind,modus
      REAL*8 :: tm,norm
cccccccc MKL STUFF
      INTEGER :: STRIDE(2),INFO
      INTEGER :: L(2)
      type(DFTI_DESCRIPTOR), POINTER :: Desc_Handle_Dim1
!      type(DFTI_DESCRIPTOR), POINTER :: Desc_Handle_Dim2
      Integer :: Status
      tm=MPI_WTIME(ierr)
c===============================================================================
! 1D FFT in X-dimension 
      psitmp=psiin  
      Status = DftiCreateDescriptor(Desc_Handle_Dim1, DFTI_DOUBLE,
     &  DFTI_COMPLEX, 1, NDX )
      Status = DftiSetValue( Desc_Handle_Dim1,
     &  DFTI_NUMBER_OF_TRANSFORMS, NDY*NDZ )
      Status = DftiSetValue( Desc_Handle_Dim1, DFTI_INPUT_DISTANCE, NDX)
      Status = DftiSetValue( Desc_Handle_Dim1, DFTI_OUTPUT_DISTANCE,NDX)
      Status = DftiCommitDescriptor(Desc_Handle_Dim1)
             mds: Select case (modus)
        Case (1)
        Status = DftiComputeForward(Desc_Handle_Dim1,psitmp)
!             write(6,*)"OK FFT MKL modified0",modus
        Case (2)
        Status = DftiComputeBackward(Desc_Handle_Dim1,psitmp)
              end  select mds
      Status = DftiFreeDescriptor(Desc_Handle_Dim1)

         fftpsi=psitmp
!         if (mod(NDX,2).eq.0) then
!           do m=(NDX/2+1),NDX
!              fftpsi(m-(NDX/2))=psitmp(m)
!              fftpsi(m)=psitmp(m-(NDX/2))
!           enddo
!         else
!           do m=(NDX/2+1),NDX
!              fftpsi(m-(NDX/2))=psitmp(m)
!              fftpsi(m)=psitmp(m-(NDX/2))
!           enddo
!         endif
!       call normvxz(fftpsi,norm,NDX)
!       call xvixdzo(norm,fftpsi,NDX)
!             write(6,*)"STR  PSI", psiin(64) 
!             write(6,*)"OK FFT MKL modified1",modus
!             write(6,*)"STR  FFT", FFTPSI(64) 
c===============================================================================
      timeCPU(3)=timeCPU(3)+MPI_WTIME(ierr)-tm
c================================================================================
      END  SUBROUTINE Get_FFTPSI_FFT_MKL_1D
