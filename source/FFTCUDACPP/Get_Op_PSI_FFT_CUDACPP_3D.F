! FFT business with MKL
! X_i index pointing to the T_x,T_y or T_z will be calculated - only for
! 2D and 3D problems
! NdimMax = 1D(NDX)
! NdimMax = 2D(NDX*NDY)
! NdimMax = 3D(NDX*NDY*NDZ)
      SUBROUTINE Get_Op_PSI_FFT_CUDACPP_3D(psi,X_i) 
      USE SHARED_DIMS
      USE DVR_ALL
      USE W_INTERPARTICLE
!      Use MKL_DFTI
      USE omp_lib
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER ::  ierr
      COMPLEX*16,   DIMENSION(NDX*NDY*NDZ) :: psi
!==================== Local Arrays for 3D Transformations 
      INTEGER :: I,J,K,I1,incr,j1,X_i,ind
      REAL*8 :: tm, tmoksa
! ccccccc MKL STUFF

      COMPLEX*16 :: ortkx(NDX),ortky(NDY),ortkz(NDZ)

      tm=MPI_WTIME(ierr)
      ortkx=ort_kx*ort_kx /(NDX*NDY*NDZ) /2/Time_mass
      ortky=ort_ky*ort_ky /(NDX*NDY*NDZ) /2/Time_mass
      ortkz=ort_kz*ort_kz /(NDX*NDY*NDZ) /2/Time_mass
         tmoksa = MPI_WTIME(ierr)
! ===============================================================================
      CALL kernel_FFT_3D(psi, ortkx, ortky, ortkz, NDX, NDY,NDZ) 
      timeCPU(3)=timeCPU(3)+MPI_WTIME(ierr)-tm
      tm=MPI_WTIME(ierr)
          tmoksa =MPI_WTIME(ierr) -tmoksa
           write(6,*)"Time _oksa 3d [sec]:", tmoksa
! ================================================================================
      END  SUBROUTINE Get_Op_PSI_FFT_CUDACPP_3D
