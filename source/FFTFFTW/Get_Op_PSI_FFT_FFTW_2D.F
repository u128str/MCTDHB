! FFT business with MKL
! X_i index pointing to the T_x,T_y or T_z will be calculated - only for
! 2D and 3D problems
! NdimMax = 1D(NDX)
! NdimMax = 2D(NDX*NDY)
! NdimMax = 3D(NDX*NDY*NDZ)
      SUBROUTINE Get_Op_PSI_FFT_FFTW_2D(psi,X_i) 
      USE SHARED_DIMS
      USE DVR_ALL
      USE W_INTERPARTICLE
      USE omp_lib
      IMPLICIT NONE
      INCLUDE 'mpif.h'
cccccccc FFTW stuff
      include 'fftw3.f'
      integer*8 plan
      integer iret,nthr
      COMPLEX*16 :: dummy(NDY),dummyX(NDX)
      INTEGER ::  ierr
      COMPLEX*16,   DIMENSION(NDX*NDY*NDZ) :: psi
c==================== Local Arrays for 3D Transformations 
      INTEGER :: I,J,K,I1,incr,j1,X_i,ind
      REAL*8 :: tm
cccccccc MKL STUFF
      tm=MPI_WTIME(ierr)
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
cccccccccc       do the FFTW forward transform
      nthr = OMP_GET_MAX_THREADS()
      call dfftw_init_threads(iret)
      call dfftw_plan_with_nthreads(nthr)
      call dfftw_plan_dft_2d(plan,NDX,NDY,PSI,PSI,
     .                FFTW_FORWARD,FFTW_ESTIMATE,fftw_unaligned)
      call dfftw_execute_dft(plan,psi,psi)
      call dfftw_destroy_plan(plan)
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
cccccccccc       do the FFTW forward transform
      nthr = OMP_GET_MAX_THREADS()
      call dfftw_init_threads(iret)
      call dfftw_plan_with_nthreads(nthr)
      call dfftw_plan_dft_2d(plan,NDX,NDY,PSI,PSI,
     .                FFTW_BACKWARD,FFTW_ESTIMATE,fftw_unaligned)
      call dfftw_execute_dft(plan,psi,psi)
      call dfftw_destroy_plan(plan)
      timeCPU(3)=timeCPU(3)+MPI_WTIME(ierr)-tm
      tm=MPI_WTIME(ierr)
      END SELECT Xi1
c================================================================================
      END  SUBROUTINE Get_Op_PSI_FFT_FFTW_2D
