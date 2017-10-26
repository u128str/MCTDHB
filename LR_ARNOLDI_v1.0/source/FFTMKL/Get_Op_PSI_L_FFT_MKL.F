! FFT business with MKL
! This subroutine applies the angular momentun operators
! so far only L_z for 2D problems
      SUBROUTINE Get_Op_PSI_L_FFT_MKL(psi,X_i) 
      USE SHARED_DIMS
      USE DVR_ALL
      USE W_INTERPARTICLE
      Use MKL_DFTI
      USE omp_lib
      IMPLICIT NONE
      INCLUDE 'mpif.h'
cccccccc FFTW stuff
      INTEGER ::  ierr
      COMPLEX*16,   DIMENSION(NDX*NDY*NDZ) :: psi,psi2
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
      COMPLEX*16 :: ortx(NDX),orty(NDY),ortz(NDZ)
      ortkx=ort_kx
      ortky=ort_ky
      ortkz=ort_kz
      ortx=ort_x
      orty=ort_y
      psi2=psi
      tm=MPI_WTIME(ierr)
c===============================================================================
      Xi1: SELECT CASE (X_i)
      CASE (1) 
c===============================================================================
      CASE (2) 
c===============================================================================

      CASE (3)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
ccccccc  do the mkl fft voodoo
      Status = DftiCreateDescriptor(Desc_Handle_Dim1, DFTI_DOUBLE,
     &  DFTI_COMPLEX, 1, NDX )
      Status = DftiSetValue( Desc_Handle_Dim1,
     &  DFTI_NUMBER_OF_TRANSFORMS, NDY*NDZ )
      Status = DftiSetValue( Desc_Handle_Dim1, DFTI_INPUT_DISTANCE, NDX)
      Status = DftiSetValue( Desc_Handle_Dim1, DFTI_OUTPUT_DISTANCE,NDX)
      Status = DftiCommitDescriptor(Desc_Handle_Dim1)
      Status = DftiComputeForward(Desc_Handle_Dim1,psi)
c================= FFT derivative: multiply by -i kx
c================= FFT derivative: multiply by -i kx
!$OMP PARALLEL DO 
!$OMP& FIRSTPRIVATE(NDX,NDY,NDZ,ortkx,Time_mass)
!$OMP& PRIVATE(i,ind,j,K)
      Do K=1,NDX
         Do I=1,NDZ
            Do J=1,NDY
              ind=k+(J-1)*NDX+(I-1)*NDX*NDY
              psi(ind)=psi(ind)*ortkx(K)*dcmplx(0.d0,-1.d0)/NDX
            EndDo
          EndDo
      EndDo
!$OMP END PARALLEL DO


ccccc Do the MKL_FFT reverse mumbo-jumbo
      Status = DftiComputeBackward(Desc_Handle_Dim1,psi)
      Status = DftiFreeDescriptor(Desc_Handle_Dim1)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCC multiply by y CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!$OMP PARALLEL DO 
!$OMP& FIRSTPRIVATE(NDX,NDY,NDZ,orty)
!$OMP& PRIVATE(i,ind,j,K)
      Do K=1,NDX
         Do I=1,NDZ
            Do J=1,NDY
              ind=k+(J-1)*NDX+(I-1)*NDX*NDY
              psi(ind)=psi(ind)*orty(J)*dcmplx(0.d0,1.d0)
            EndDo
          EndDo
      EndDo
!$OMP END PARALLEL DO

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

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
      Status = DftiComputeForward(Desc_Handle_Dim2,psi2)
           ENDDO
c================= Kinetic FFT: multiply FT by k**2
!$OMP PARALLEL DO 
!$OMP& FIRSTPRIVATE(NDX,NDY,NDZ,ortky,Time_mass)
!$OMP& PRIVATE(i,ind,j,K)
      Do J=1,NDY
          Do I=1,NDZ
            Do K=1,NDX
              ind=k+(J-1)*NDX+(I-1)*NDX*NDY
              psi2(ind)=psi2(ind)*dcmplx(0.d0,-1.d0)*ortky(J)/NDY
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
      Status = DftiComputeBackward(Desc_Handle_Dim2,psi2)
           ENDDO
      Status = DftiFreeDescriptor(Desc_Handle_Dim2)
      timeCPU(3)=timeCPU(3)+MPI_WTIME(ierr)-tm

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCC multiply by x CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!$OMP PARALLEL DO 
!$OMP& FIRSTPRIVATE(NDX,NDY,NDZ,ortx)
!$OMP& PRIVATE(i,ind,j,K)
      Do J=1,NDY
          Do I=1,NDZ
            Do K=1,NDX
              ind=k+(J-1)*NDX+(I-1)*NDX*NDY
              psi2(ind)=dcmplx(0.d0,-1.d0)*ortx(k)*psi2(ind)
            EndDo
          EndDo
       EndDo
!$OMP END PARALLEL DO

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCC now add the contributions x(d/dy) and y(d/dx) ccccccccccccccc

       psi=psi+psi2
c===============================================================================
ccccccc  do the mkl fft voodoo
      timeCPU(3)=timeCPU(3)+MPI_WTIME(ierr)-tm
c          stop
      CASE (4) ! True  2D FFT
          stop
      END SELECT Xi1
c================================================================================
      END  SUBROUTINE Get_Op_PSI_L_FFT_MKL
