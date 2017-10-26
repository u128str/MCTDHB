! FFT business with MKL
! X_i index pointing to the T_x,T_y or T_z will be calculated - only for
! 2D and 3D problems
! NdimMax = 1D(NDX)
! NdimMax = 2D(NDX*NDY)
! NdimMax = 3D(NDX*NDY*NDZ)
      SUBROUTINE Get_Op_PSI_FFT_CUDACPP_2D(psi,X_i) 
      USE SHARED_DIMS
      USE DVR_ALL
      USE W_INTERPARTICLE

!      Use MKL_DFTI

      USE omp_lib
      IMPLICIT NONE
      INCLUDE 'mpif.h'

        INTEGER ::  ierr
      COMPLEX*16,   DIMENSION(NDX*NDY*NDZ) :: psi
      
      COMPLEX*16 :: ortkx(NDX),ortky(NDY)          ! add oksa

!==================== Local Arrays for 3D Transformations 
      INTEGER :: I,J,K,I1,incr,j1,X_i,ind
      REAL*8 :: tm, tmoksa
!cccccccc MKL STUFF
   
      tm=MPI_WTIME(ierr)


      timeCPU(3)=timeCPU(3)+MPI_WTIME(ierr)-tm
      tm=MPI_WTIME(ierr)

! ------------------------- CUDA CPP------------------------
          ortkx= ort_kx*ort_kx/NDX/2/Time_mass/NDY 
          ortky= ort_ky*ort_ky/NDX/2/Time_mass/NDY
!           write(6,*) "I am entering OKSA" 

                 tmoksa=MPI_WTIME(ierr)
!===============================================================================
      Xi1: SELECT CASE (X_i)
      CASE (1) ! 1D FFT in X-dimension 
          stop
      CASE (2) !1D FFT in Y-dimension
          stop
      CASE (3)   
         stop

!           CALL kernel_FFT_2D_XY(psi, ortkx, ortky, NDX, NDY)

      CASE (4) ! True  2D FFT

          CALL kernel_FFT_2D(psi, ortkx, ortky, NDX, NDY)
!           write(6,*) "I OKSA is done!" 
 
      END SELECT Xi1
!================================================================================
!-------------------------end CUDA CPP ---------------------
            
      timeCPU(4)=timeCPU(4)+MPI_WTIME(ierr)-tm
      tm=MPI_WTIME(ierr)
      timeCPU(3)=timeCPU(3)+MPI_WTIME(ierr)-tm
      tm=MPI_WTIME(ierr)
                      tmoksa =MPI_WTIME(ierr) -tmoksa
	   write(6,*)"Time _oksa 2d [sec]:", tmoksa

! ================================================================================
      END  SUBROUTINE Get_Op_PSI_FFT_CUDACPP_2D
