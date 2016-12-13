      subroutine computewsl_FFTWIMEST(wsl,psi_s,psi_l,vtilde)
      use SHARED_DIMS 
      use DVR_ALL
!      Use MKL_DFTI
      USE omp_lib

      implicit none
      double complex     :: Wsl(NDX*NDY*NDZ)
      double complex  :: psi_s(NDX*NDY*NDZ)
      double complex  :: psi_l(NDX*NDY*NDZ)
      double complex  :: vtilde(NDX*NDY*NDZ)

      integer                 :: i,j,k,m,n
      double complex          :: fsltilde(NDX*NDY*NDZ)
!---------------FFT-------------------------------
       include 'fftw3.f'
       integer*8 plan,plany
       integer iret,nthr
       INTEGER ::  ierr
       INTEGER :: FFTW_USE_WISDOM,isuccess
c=================================
!       write(6,*)"FFTW IMEST", DIM_MCTDHB
!compute f_sl = dx phi_s^* phi_l
       do m=1,NDX*NDY*NDZ
         fsltilde(m)=DCONJG(psi_s(m))*psi_l(m)
       enddo

!get f_sl-tilde dx from f_sl dx
!   compute forward fft of fsltilde
      if (DIM_MCTDHB.eq.1) then
      call dfftw_plan_dft_1d(plan,NDX,fsltilde,fsltilde,
     .                       FFTW_FORWARD,
     .                       FFTW_ESTIMATE)
      call dfftw_execute_dft(plan,fsltilde,fsltilde)
      call dfftw_destroy_plan(plan)

      else if (DIM_MCTDHB.eq.2) then
      nthr = OMP_GET_MAX_THREADS()
!       write(6,*)"OMP_GET_MAX_THREADS", nthr
      call dfftw_init_threads(iret)
      call dfftw_plan_with_nthreads(nthr)
      call dfftw_plan_dft_2d(plan,NDX,NDY,fsltilde,fsltilde,
     .                FFTW_FORWARD,FFTW_ESTIMATE,fftw_unaligned)
      call dfftw_execute_dft(plan,fsltilde,fsltilde)
      call dfftw_destroy_plan(plan)



      else if (DIM_MCTDHB.eq.3) then
      call dfftw_plan_dft_3d(plan,NDX,NDY,NDZ,fsltilde,fsltilde,
     .                FFTW_FORWARD,FFTW_ESTIMATE,fftw_unaligned)
      call dfftw_execute_dft(plan,fsltilde,fsltilde)
      call dfftw_destroy_plan(plan)
      endif

c multiply it with the FFT of the potential
      do m=1,NDX*NDY*NDZ
         Wsl(m)=fsltilde(m)*vtilde(m)/NDX/NDY/NDZ
      enddo

c       transform back 
      if (DIM_MCTDHB.eq.1) then
      call dfftw_plan_dft_1d(plan,NDX,wsl,wsl,
     .                       FFTW_BACKWARD,
     .                       FFTW_ESTIMATE)
      call dfftw_execute_dft(plan,wsl,wsl)
      call dfftw_destroy_plan(plan)

      else if (DIM_MCTDHB.eq.2) then
      nthr = OMP_GET_MAX_THREADS()
      call dfftw_init_threads(iret)
      call dfftw_plan_with_nthreads(nthr)
      call dfftw_plan_dft_2d(plan,NDX,NDY,wsl,wsl,
     .                FFTW_BACKWARD,FFTW_ESTIMATE,fftw_unaligned)
      call dfftw_execute_dft(plan,wsl,wsl)
      call dfftw_destroy_plan(plan)

      else if (DIM_MCTDHB.eq.3) then
      call dfftw_plan_dft_3d(plan,NDX,NDY,NDZ,wsl,wsl,
     .                FFTW_BACKWARD,FFTW_ESTIMATE,fftw_unaligned)
      call dfftw_execute_dft(plan,wsl,wsl)
      call dfftw_destroy_plan(plan)
      endif


       return
       end subroutine computewsl_FFTWIMEST

