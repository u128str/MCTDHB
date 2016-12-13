! FFT business with MKL
! X_i index pointing to the T_x,T_y or T_z will be calculated - only for
! 2D and 3D problems
! NdimMax = 1D(NDX)
! NdimMax = 2D(NDX*NDY)
! NdimMax = 3D(NDX*NDY*NDZ)
      SUBROUTINE Get_Op_PSI_FFT_FFTW_3D(psi,X_i) 
      USE SHARED_DIMS
      USE DVR_ALL
      USE W_INTERPARTICLE
      USE omp_lib
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      include 'fftw3.f'
      integer*8 plan
      integer iret,nthr
      COMPLEX*16 :: dummy(NDY),dummyX(NDX)
cccccccc FFTW stuff
      INTEGER ::  ierr
      COMPLEX*16,   DIMENSION(NDX*NDY*NDZ) :: psi
c==================== Local Arrays for 3D Transformations 
      INTEGER :: I,J,K,I1,incr,j1,X_i,ind
      REAL*8 :: tm
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
      nthr = OMP_GET_MAX_THREADS()
      call dfftw_init_threads(iret)
      call dfftw_plan_with_nthreads(nthr)
      call dfftw_plan_dft_3d(plan,NDX,NDY,NDZ,PSI,PSI,
     .                FFTW_FORWARD,FFTW_ESTIMATE,fftw_unaligned)
      call dfftw_execute_dft(plan,psi,psi)
      call dfftw_destroy_plan(plan)
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
cccccccccc       do the FFTW forward transform
      nthr = OMP_GET_MAX_THREADS()
      call dfftw_init_threads(iret)
      call dfftw_plan_with_nthreads(nthr)
      call dfftw_plan_dft_3d(plan,NDX,NDY,NDZ,PSI,PSI,
     .                FFTW_BACKWARD,FFTW_ESTIMATE,fftw_unaligned)
      call dfftw_execute_dft(plan,psi,psi)
      call dfftw_destroy_plan(plan)
      timeCPU(3)=timeCPU(3)+MPI_WTIME(ierr)-tm
      tm=MPI_WTIME(ierr)
      END SELECT Xi1
c================================================================================
      END  SUBROUTINE Get_Op_PSI_FFT_FFTW_3D
