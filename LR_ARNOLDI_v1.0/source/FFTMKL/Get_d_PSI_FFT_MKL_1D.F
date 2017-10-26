! Added by STR Feb 6 2015
! FFT business with MKL
! X_i index pointing to the T_x,T_y or T_z will be calculated - only for
! 2D and 3D problems
! NdimMax = 1D(NDX)
! NdimMax = 2D(NDX*NDY)
! NdimMax = 3D(NDX*NDY*NDZ)
      SUBROUTINE Get_d_PSI_FFT_MKL_1D(psi,X_i) 
      USE SHARED_DIMS
      USE DVR_ALL
      USE W_INTERPARTICLE
      Use MKL_DFTI
      USE omp_lib
      IMPLICIT NONE
      INCLUDE 'mpif.h'
cccccccc FFTW stuff
      INTEGER ::  ierr
      COMPLEX*16,   DIMENSION(NDX*NDY*NDZ) :: psi
c==================== Local Arrays for 3D Transformations 
      INTEGER :: I,J,K,I1,incr,j1,X_i,ind
      REAL*8 :: tm
cccccccc MKL STUFF
      INTEGER :: STRIDE(2),INFO
      INTEGER :: L(2)
      type(DFTI_DESCRIPTOR), POINTER :: Desc_Handle_Dim1
      type(DFTI_DESCRIPTOR), POINTER :: Desc_Handle_Dim2
      Integer :: Status
      COMPLEX*16 :: ortkx(NDX),ortky(NDY),ortkz(NDZ)
      ortkx=ort_kx
      ortky=ort_ky
      ortkz=ort_kz
      tm=MPI_WTIME(ierr)
c===============================================================================
      Xi1: SELECT CASE (X_i)
      CASE (1) ! 1D FFT in X-dimension 
ccccccc  do the mkl fft voodoo
      Status = DftiCreateDescriptor(Desc_Handle_Dim1, DFTI_DOUBLE,
     &  DFTI_COMPLEX, 1, NDX )
      Status = DftiSetValue( Desc_Handle_Dim1,
     &  DFTI_NUMBER_OF_TRANSFORMS, NDY*NDZ )
      Status = DftiSetValue( Desc_Handle_Dim1, DFTI_INPUT_DISTANCE, NDX)
      Status = DftiSetValue( Desc_Handle_Dim1, DFTI_OUTPUT_DISTANCE,NDX)
      Status = DftiCommitDescriptor(Desc_Handle_Dim1)
      Status = DftiComputeForward(Desc_Handle_Dim1,psi)
c================= Kinetic FFT: multiply by k**2
c!$OMP PARALLEL DO 
c!$OMP& FIRSTPRIVATE(NDX,ortkx,Time_mass)
c!$OMP& PRIVATE(i,ind,j)
c      Do i=1,NDX
c         ind=i
c         Do j=1,NDY*NDZ
c            ind=i+NDX*(J-1)
c            psi(ind)=psi(ind)*ortkx(i)**2/(2*Time_mass)/NDX
c         EndDo
c      EndDo
c!$OMP END PARALLEL DO
!$OMP PARALLEL DO 
!$OMP& FIRSTPRIVATE(NDX,NDY,NDZ,ortkx,Time_mass)
!$OMP& PRIVATE(i,ind,j,K)
            Do K=1,NDX
            Do I=1,NDZ
            Do J=1,NDY
            ind=k+(J-1)*NDX+(I-1)*NDX*NDY
            psi(ind)=(0.0d0,1.0d0)*psi(ind)*ortkx(K)/NDX
            EndDo
            EndDo
            EndDo
!$OMP END PARALLEL DO
c            ind=1
c            Do I=1,NDZ
c            Do J=1,NDY
c            Do K=1,NDX
c            psi(ind)=psi(ind)*ortkx(K)**2/(2*Time_mass)/NDX
c            ind=ind+1
c            EndDo
c            EndDo
c            EndDo
c MKLFFT backtransform
      Status = DftiComputeBackward(Desc_Handle_Dim1,psi)
      Status = DftiFreeDescriptor(Desc_Handle_Dim1)
      timeCPU(3)=timeCPU(3)+MPI_WTIME(ierr)-tm
c===============================================================================
      CASE (2) !1D FFT in Y-dimension
c          tm=MPI_WTIME(ierr)

      Status = DftiCreateDescriptor(Desc_Handle_Dim2, DFTI_DOUBLE,
     &  DFTI_COMPLEX,1,NDY)
c      Stride(1) = 0 
      Stride(2) = NDX
      Status = DftiSetValue( Desc_Handle_Dim2,
     &  DFTI_NUMBER_OF_TRANSFORMS, NDX)
      Status = DftiSetValue(Desc_Handle_Dim2,DFTI_INPUT_DISTANCE, 1)
      Status = DftiSetValue(Desc_Handle_Dim2,DFTI_OUTPUT_DISTANCE,1)
           DO I=1,NDX*NDY*NDZ,NDX*NDY
      Stride(1) = I-1
      Status=  DftiSetValue(Desc_Handle_Dim2,DFTI_INPUT_STRIDES,Stride) 
      Status=  DftiSetValue(Desc_Handle_Dim2,DFTI_OUTPUT_STRIDES,Stride)
      Status = DftiCommitDescriptor(Desc_Handle_Dim2)
      Status = DftiComputeForward(Desc_Handle_Dim2,psi)
           ENDDO
c================= Kinetic FFT: multiply FT by k**2
c!$OMP PARALLEL DO
c      Do i=1,NDX
c         ind=i
c         Do j=1,NDY
c           ind=i+NDX*(J-1)
c           psi(ind)=psi(ind)*ort_ky(j)**2/(2*Time_mass)/NDY
c         EndDo
c      EndDo
!$OMP PARALLEL DO 
!$OMP& FIRSTPRIVATE(NDX,NDY,NDZ,ortky,Time_mass)
!$OMP& PRIVATE(i,ind,j,K)
            Do J=1,NDY
            Do I=1,NDZ
            Do K=1,NDX
            ind=k+(J-1)*NDX+(I-1)*NDX*NDY
            psi(ind)=(0.0d0,1.0d0)*psi(ind)*ortky(J)/NDY
            EndDo
            EndDo
            EndDo
!$OMP END PARALLEL DO
ccccc Do the MKL_FFT reverse mumbo-jumbo
           DO I=1,NDX*NDY*NDZ,NDX*NDY
      Stride(1) = I-1
      Status=  DftiSetValue(Desc_Handle_Dim2,DFTI_INPUT_STRIDES,Stride) 
      Status=  DftiSetValue(Desc_Handle_Dim2,DFTI_OUTPUT_STRIDES,Stride)
      Status = DftiCommitDescriptor(Desc_Handle_Dim2)
      Status = DftiComputeBackward(Desc_Handle_Dim2,psi)
           ENDDO
      Status = DftiFreeDescriptor(Desc_Handle_Dim2)
      timeCPU(3)=timeCPU(3)+MPI_WTIME(ierr)-tm

c===============================================================================
      CASE (3)
c          tm=MPI_WTIME(ierr)
      Status = DftiCreateDescriptor(Desc_Handle_Dim1, DFTI_DOUBLE,
     &  DFTI_COMPLEX, 1, NDZ )
      Stride(1) = 0 
      Stride(2) = NDX*NDY
      Status = DftiSetValue( Desc_Handle_Dim1,
     &  DFTI_NUMBER_OF_TRANSFORMS, NDX*NDY )
      Status = DftiSetValue(Desc_Handle_Dim1,DFTI_INPUT_DISTANCE, 1)
      Status = DftiSetValue(Desc_Handle_Dim1,DFTI_OUTPUT_DISTANCE,1)
      Status=  DftiSetValue(Desc_Handle_Dim1,DFTI_INPUT_STRIDES,Stride) 
      Status=  DftiSetValue(Desc_Handle_Dim1,DFTI_OUTPUT_STRIDES,Stride)
      Status = DftiCommitDescriptor(Desc_Handle_Dim1)
      Status = DftiComputeForward(Desc_Handle_Dim1,psi)
!$OMP PARALLEL DO 
!$OMP& FIRSTPRIVATE(NDX,NDY,NDZ,ortkz,Time_mass)
!$OMP& PRIVATE(i,ind,j,K)
            Do I=1,NDZ
            Do J=1,NDY
            Do K=1,NDX
            ind=k+(J-1)*NDX+(I-1)*NDX*NDY
            psi(ind)=(0.0d0,1.0d0)*psi(ind)*ortkz(I)/NDZ
            EndDo
            EndDo
            EndDo
!$OMP END PARALLEL DO
ccccc Do the MKL_FFT reverse mumbo-jumbo
      Status = DftiComputeBackward(Desc_Handle_Dim1,psi)
      Status = DftiFreeDescriptor(Desc_Handle_Dim1)
      timeCPU(3)=timeCPU(3)+MPI_WTIME(ierr)-tm
c          stop
      CASE (4) ! True  2D FFT
          stop
      END SELECT Xi1
c================================================================================
      END  SUBROUTINE Get_d_PSI_FFT_MKL_1D
