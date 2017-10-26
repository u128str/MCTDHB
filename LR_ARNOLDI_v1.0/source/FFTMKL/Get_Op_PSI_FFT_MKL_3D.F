! FFT business with MKL
! X_i index pointing to the T_x,T_y or T_z will be calculated - only for
! 2D and 3D problems
! NdimMax = 1D(NDX)
! NdimMax = 2D(NDX*NDY)
! NdimMax = 3D(NDX*NDY*NDZ)
      SUBROUTINE Get_Op_PSI_FFT_MKL_3D(psi,X_i) 
      USE SHARED_DIMS
      USE DVR_ALL
      USE W_INTERPARTICLE
      Use MKL_DFTI
      USE omp_lib
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER ::  ierr
      COMPLEX*16,   DIMENSION(NDX*NDY*NDZ) :: psi
c==================== Local Arrays for 3D Transformations 
      INTEGER :: I,J,K,I1,incr,j1,X_i,ind
      REAL*8 :: tm
cccccccc MKL STUFF
      INTEGER :: INFO
      INTEGER :: L(3)
      type(DFTI_DESCRIPTOR), POINTER :: Desc_Handle_Dim1
      Integer :: Status
      COMPLEX*16 :: ortky(NDY),ortkX(NDX),ortkz(NDZ)
      ortkx=ort_kx
      ortky=ort_ky
      ortkz=ort_kz
      tm=MPI_WTIME(ierr)
c===============================================================================
      Xi1: SELECT CASE (X_i)
      CASE (1) ! 1D FFT in X-dimension 
          stop
      CASE (2) !1D FFT in Y-dimension
          stop
      CASE (3)
          stop
      CASE (4)
          stop
      CASE (5) ! True  3D FFT
c#########################################################################
c#########################################################################
      L(1)=NDX
      L(2)=NDY
      L(3)=NDZ
      Status = DftiCreateDescriptor(Desc_Handle_Dim1,DFTI_DOUBLE,
     &  DFTI_COMPLEX, 3, L )
      Status = DftiCommitDescriptor(Desc_Handle_Dim1) 
      Status = DftiComputeForward(Desc_Handle_Dim1,psi)
      timeCPU(3)=timeCPU(3)+MPI_WTIME(ierr)-tm
      tm=MPI_WTIME(ierr)
      ind=1
c!$OMP PARALLEL DO PRIVATE(i,j,k,ind) 
c!$OMP& Firstprivate(Time_mass,ortkx,ortky,ortkz,NDY,NDX,ndz)
      Do i=1,NDZ
         Do j=1,NDY
            do k=1,NDX            
             ind=k+(j-1)*NDX+(i-1)*NDX*NDY
            psi(ind)=(ortkx(k)**2+ ortky(j)**2+ortkz(i)**2)
     &                  /(2*Time_mass)/NDX/NDY/NDZ*
     &            psi(ind)
            enddo
         EndDo
      EndDo
c!$OMP END PARALLEL DO
      timeCPU(4)=timeCPU(4)+MPI_WTIME(ierr)-tm
      tm=MPI_WTIME(ierr)
      Status = DftiComputeBackward(Desc_Handle_Dim1,psi)
      Status = DftiFreeDescriptor(Desc_Handle_Dim1)
      timeCPU(3)=timeCPU(3)+MPI_WTIME(ierr)-tm
      tm=MPI_WTIME(ierr)
      END SELECT Xi1
c================================================================================
      END  SUBROUTINE Get_Op_PSI_FFT_MKL_3D
