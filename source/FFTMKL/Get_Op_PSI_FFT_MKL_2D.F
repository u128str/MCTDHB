! FFT business with MKL
! X_i index pointing to the T_x,T_y or T_z will be calculated - only for
! 2D and 3D problems
! NdimMax = 1D(NDX)
! NdimMax = 2D(NDX*NDY)
! NdimMax = 3D(NDX*NDY*NDZ)
      SUBROUTINE Get_Op_PSI_FFT_MKL_2D(psi,X_i) 
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

c  oksa add
      INTEGER ::  Ntmp
	  REAL*8 ::   tmoksa
cccccccc MKL STUFF
      INTEGER :: STRIDE(2),INFO
      INTEGER :: L(2)
      type(DFTI_DESCRIPTOR), POINTER :: Desc_Handle_Dim1
      Integer :: Status
      tm=MPI_WTIME(ierr)
  
      tmoksa= omp_get_wtime()
c  oksa time and omp_get_max_threads()
     
        Ntmp= omp_get_max_threads()
!       write(6,*) "I OKSA_MKL : max threads =!", Ntmp
c===============================================================================
      Xi1: SELECT CASE (X_i)
      CASE (1) ! 1D FFT in X-dimension 
          stop
      CASE (2) !1D FFT in Y-dimension
          stop
      CASE (3)
          stop
      CASE (4) ! True  2D FFT


c#########################################################################
c#########################################################################
      L(1)=NDX
      L(2)=NDY
      Status = DftiCreateDescriptor(Desc_Handle_Dim1, DFTI_DOUBLE,
     &  DFTI_COMPLEX, 2, L )
      Status = DftiCommitDescriptor( Desc_Handle_Dim1) 
      Status = DftiComputeForward(Desc_Handle_Dim1,psi)
      timeCPU(3)=timeCPU(3)+MPI_WTIME(ierr)-tm
      tm=MPI_WTIME(ierr)

!$OMP PARALLEL DO
!$OMP& PRIVATE(i,ind,j)
!$OMP& shared(psi,Time_mass,ort_kx,ort_ky,NDY,NDX)
      Do i=1,NDX
         Do j=1,NDY*NDZ
           ind=i+NDX*(J-1)
            psi(ind)=(ort_kx(i)**2+ ort_ky(j)**2)/(2*Time_mass)/NDX/NDY*
     &            psi(ind)
         EndDo
      EndDo
!$OMP END PARALLEL DO

      timeCPU(4)=timeCPU(4)+MPI_WTIME(ierr)-tm
      tm=MPI_WTIME(ierr)
      Status = DftiComputeBackward(Desc_Handle_Dim1,psi)
      Status = DftiFreeDescriptor(Desc_Handle_Dim1)
      timeCPU(3)=timeCPU(3)+MPI_WTIME(ierr)-tm
      tm=MPI_WTIME(ierr)
c             oksa add time 
      tmoksa= omp_get_wtime() -tmoksa
!      write(6,*) "I OKSA : time MKL 2d =",   tmoksa
      END SELECT Xi1
c================================================================================
      END  SUBROUTINE Get_Op_PSI_FFT_MKL_2D
