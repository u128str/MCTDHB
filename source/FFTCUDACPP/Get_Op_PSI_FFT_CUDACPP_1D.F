! FFT business with CUDA
! X_i index pointing to the T_x,T_y or T_z will be calculated - only for
! 2D and 3D problems
! NdimMax = 1D(NDX)
! NdimMax = 2D(NDX*NDY)
! NdimMax = 3D(NDX*NDY*NDZ)
      SUBROUTINE Get_Op_PSI_FFT_CUDACPP_1D(psi,X_i) 
      USE SHARED_DIMS
      USE DVR_ALL
      USE W_INTERPARTICLE
!      Use MKL_DFTI
 
!      USE precision_m
!      USE cufft_m
!      USE GPU_ALL

    
      USE omp_lib
      IMPLICIT NONE
      INCLUDE 'mpif.h'


      INTEGER ::  ierr
      COMPLEX*16,   DIMENSION(NDX*NDY*NDZ) :: psi
!==================== Local Arrays for 3D Tnsformations 
      INTEGER :: I,J,K,I1,incr,j1,X_i,ind,NDfdim
      REAL*8 :: tm
!================== MKL STUFF
      INTEGER :: STRIDE(2),INFO
      INTEGER :: L(2)
!      type(DFTI_DESCRIPTOR), POINTER :: Desc_Handle_Dim1
!      type(DFTI_DESCRIPTOR), POINTER :: Desc_Handle_Dim2
      Integer :: Status
   
      COMPLEX*16 :: ortkx(NDX),ortky(NDY),ortkz(NDZ)

  
         integer :: plan, planType

      ortkx=ort_kx
      ortky=ort_ky
      ortkz=ort_kz
   
      tm=MPI_WTIME(ierr)

!-------------------- variabls on divece(GPU) ---------------------------
                 
!STR I have ADDED NDfdim to define the size, dimension of the PSI array
         NDfdim=ndx*ndy*ndz

!======================================================================== 
                      
        
           

      Xi1: SELECT CASE (X_i)
      CASE (1) ! 1D FFT in X-dimension 
!================================================================================ 
               
              ortkx= ort_kx*ort_kx/NDX/2/Time_mass  
                          
              CALL kernel_FFT_1D(psi, ortkx,  NDX, NDY,NDZ, 1)

             timeCPU(3)=timeCPU(3)+MPI_WTIME(ierr)-tm 

!c===============================================================================
      CASE (2) !1D FFT in Y-dimension
           
        ortky= ort_ky*ort_ky/NDY/2/Time_mass

        CALL kernel_FFT_1D(psi, ortky, NDX, NDY,NDZ, 2)
     
        timeCPU(3)=timeCPU(3)+MPI_WTIME(ierr)-tm

!c===============================================================================
      CASE (3)

        ortkz= ort_kz*ort_kz/NDZ/2/Time_mass

        CALL kernel_FFT_1D(psi, ortkz, NDX, NDY,NDZ, 3)

      timeCPU(3)=timeCPU(3)+MPI_WTIME(ierr)-tm

!c          stop
      CASE (4) ! True  2D FFT
          stop
      END SELECT Xi1
!c================================================================================
      END  SUBROUTINE Get_Op_PSI_FFT_CUDACPP_1D
