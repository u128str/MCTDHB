! FFT business with MKL
! X_i index pointing to the T_x,T_y or T_z will be calculated - only for
! 2D and 3D problems
! NdimMax = 1D(NDX)
! NdimMax = 2D(NDX*NDY)
! NdimMax = 3D(NDX*NDY*NDZ)
      SUBROUTINE Get_Op_PSI_FFT_FFTW_1D(psi,X_i) 
      USE SHARED_DIMS
      USE DVR_ALL
      USE W_INTERPARTICLE
      USE omp_lib
      IMPLICIT NONE
      INCLUDE 'mpif.h'
cccccccc FFTW stuff
      include 'fftw3.f'
      integer*8 plan,plany
      integer iret,nthr
      COMPLEX*16 :: dummyX(NDX),dummyY(NDY),dummyZ(NDZ)
      INTEGER ::  ierr
      INTEGER :: FFTW_USE_WISDOM,isuccess
      COMPLEX*16,   DIMENSION(NDX*NDY*NDZ) :: psi

c==================== Local Arrays for 3D Transformations 
      INTEGER :: I,J,K,I1,incr,j1,X_i,ind
      REAL*8 :: tm
      COMPLEX*16 :: ortky(NDY),ortkX(NDX),ortkz(NDZ)
      ortkx=ort_kx
      ortky=ort_ky
      ortkz=ort_kz
c===============================================================================
!       write(6,*) "In FFTW0",X_i
      Xi1: SELECT CASE (X_i)
      CASE (1) ! 1D FFT in X-dimension 
      tm=MPI_WTIME(ierr)
      call dfftw_plan_dft_1d(plan,NDX,dummyX,dummyX,
     .                       FFTW_FORWARD,
     .                       FFTW_ESTIMATE)
c!$OMP PARALLEL DO
c!$OMP& PRIVATE(dummyx,ind)
c!$OMP& FIRSTPRIVATE(plan,NDY,NDX,NDZ)
      do k=1,NDZ
      do i=1,NDY
         ind=NDX*NDY*(k-1)+i
         do j=1,NDX
            ind=NDX*NDY*(k-1)+(i-1)*NDX+j
            dummyX(j)=psi(ind)
         end do
         call dfftw_execute_dft(plan,dummyX,dummyX)
         do j=1,NDX
            ind=NDX*NDY*(k-1)+(i-1)*NDX+j
            PSI(ind)=dummyX(j)
         end do
      end do
      end do
c!$OMP END PARALLEL DO
      call dfftw_destroy_plan(plan)
      timeCPU(3)=timeCPU(3)+MPI_WTIME(ierr)-tm
      tm=MPI_WTIME(ierr)
c================= Kinetic FFT: multiply by k**2
!$OMP PARALLEL DO 
!$OMP& FIRSTPRIVATE(NDX,NDY,NDZ,ortkx,Time_mass)
!$OMP& PRIVATE(i,ind,j,K)
            Do K=1,NDX
            Do I=1,NDZ
            Do J=1,NDY
            ind=k+(J-1)*NDX+(I-1)*NDX*NDY
            psi(ind)=psi(ind)*ortkx(K)**2/(2*Time_mass)/NDX
            EndDo
            EndDo
            EndDo
!$OMP END PARALLEL DO
      timeCPU(4)=timeCPU(4)+MPI_WTIME(ierr)-tm
      tm=MPI_WTIME(ierr)
      call dfftw_plan_dft_1d(plan,NDX,dummyX,dummyX,
     .                       FFTW_BACKWARD,
     .                       FFTW_ESTIMATE)
      do k=1,NDZ
      do i=1,NDY
         ind=NDX*NDY*(k-1)+i
         do j=1,NDX
            ind=NDX*NDY*(k-1)+(i-1)*NDX+j
            dummyX(j)=psi(ind)
         end do
         call dfftw_execute_dft(plan,dummyX,dummyX)
         do j=1,NDX
            ind=NDX*NDY*(k-1)+(i-1)*NDX+j
            PSI(ind)=dummyX(j)
         end do
      end do
      end do
      call dfftw_destroy_plan(plan)
      timeCPU(3)=timeCPU(3)+MPI_WTIME(ierr)-tm
c=========================Y dimensionYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY
      CASE (2) !1D FFT in Y-dimension
          tm=MPI_WTIME(ierr)
      call dfftw_plan_dft_1d(plan,NDY,dummyY,dummyY,
     .                       FFTW_FORWARD,
     .                       FFTW_ESTIMATE)
      do k=1,NDZ
      do i=1,NDX
         do j=1,NDY
            ind=i+NDX*(J-1)+NDX*NDY*(k-1)
            dummyY(j)=psi(ind)
         end do  
         call dfftw_execute_dft(plan,dummyY,dummyY)
         do j=1,NDY
            ind=i+NDX*(J-1)+NDX*NDY*(k-1)
            PSI(ind)=dummyY(j)
         end do           
      end do  
      end do  
      call dfftw_destroy_plan(plan)   

      timeCPU(3)=timeCPU(3)+MPI_WTIME(ierr)-tm
      tm=MPI_WTIME(ierr)
c================= Kinetic FFT: multiply FT by k**2
!$OMP PARALLEL DO 
!$OMP& FIRSTPRIVATE(NDX,NDY,NDZ,ortky,Time_mass)
!$OMP& PRIVATE(i,ind,j,K)
            Do J=1,NDY
            Do I=1,NDZ
            Do K=1,NDX
            ind=k+(J-1)*NDX+(I-1)*NDX*NDY
            psi(ind)=psi(ind)*ortky(J)**2/(2*Time_mass)/NDY
            EndDo
            EndDo
            EndDo
!$OMP END PARALLEL DO
      call dfftw_plan_dft_1d(plan,NDY,dummyY,dummyY,
     .                       FFTW_BACKWARD,
     .                       FFTW_ESTIMATE)
      do k=1,NDZ
      do i=1,NDX
         do j=1,NDY
            ind=i+NDX*(J-1)+NDX*NDY*(k-1)
            dummyY(j)=psi(ind)
         end do  
         call dfftw_execute_dft(plan,dummyY,dummyY)
         do j=1,NDY
            ind=i+NDX*(J-1)+NDX*NDY*(k-1)
            PSI(ind)=dummyY(j)
         end do           
      end do  
      end do  
      call dfftw_destroy_plan(plan)
      timeCPU(3)=timeCPU(3)+MPI_WTIME(ierr)-tm

!========================= Z-direction
      CASE (3) !1D FFT in Z-dimension
          tm=MPI_WTIME(ierr)
      call dfftw_plan_dft_1d(plan,NDZ,dummyZ,dummyZ,
     .                       FFTW_FORWARD,
     .                       FFTW_ESTIMATE)
      do i=1,NDX
      do k=1,NDY
         do j=1,NDZ
            ind=i+NDX*(k-1)+NDX*NDY*(j-1)
            dummyZ(j)=psi(ind)
         end do  
         call dfftw_execute_dft(plan,dummyZ,dummyZ)
         do j=1,NDZ
            ind=i+NDX*(k-1)+NDX*NDY*(j-1)
            PSI(ind)=dummyZ(j)
         end do           
      end do  
      end do  
      call dfftw_destroy_plan(plan)   
      timeCPU(3)=timeCPU(3)+MPI_WTIME(ierr)-tm
      tm=MPI_WTIME(ierr)
c================= Kinetic FFT: multiply FT by k**2
!$OMP PARALLEL DO 
!$OMP& FIRSTPRIVATE(NDX,NDY,NDZ,ortkz,Time_mass)
!$OMP& PRIVATE(i,ind,j,K)
            Do I=1,NDZ
            Do J=1,NDY
            Do K=1,NDX
            ind=k+(J-1)*NDX+(I-1)*NDX*NDY
            psi(ind)=psi(ind)*ortkz(I)**2/(2*Time_mass)/NDZ
            EndDo
            EndDo
            EndDo
!$OMP END PARALLEL DO


      call dfftw_plan_dft_1d(plan,NDZ,dummyZ,dummyZ,
     .                       FFTW_BACKWARD,
     .                       FFTW_ESTIMATE)
      do k=1,NDY
      do i=1,NDX
         do j=1,NDZ
            ind=i+NDX*(k-1)+NDX*NDY*(j-1)
            dummyZ(j)=psi(ind)
         end do  
         call dfftw_execute_dft(plan,dummyZ,dummyZ)
         do j=1,NDZ
            ind=i+NDX*(k-1)+NDX*NDY*(j-1)
            PSI(ind)=dummyZ(j)
         end do           
      end do  
      end do  
      call dfftw_destroy_plan(plan)
      timeCPU(3)=timeCPU(3)+MPI_WTIME(ierr)-tm
      END SELECT Xi1
c================================================================================
      END  SUBROUTINE Get_Op_PSI_FFT_FFTW_1D
