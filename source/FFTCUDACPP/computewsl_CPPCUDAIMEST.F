      subroutine computewsl_CPPCUDAIMEST(Wsl,psi_s,psi_l,vtilde)
      use SHARED_DIMS 
      use DVR_ALL
!      Use MKL_DFTI
      USE omp_lib

       !---{{{
      implicit none
      double complex     :: Wsl(NDX*NDY*NDZ)
      double complex  :: psi_s(NDX*NDY*NDZ)
      double complex  :: psi_l(NDX*NDY*NDZ)
      double complex  :: vtilde(NDX*NDY*NDZ)

      integer                 :: i,j,k,m,n
      double complex          :: fsltilde(NDX*NDY*NDZ)

!---------------FFT-------------------------------
cccccccc MKL STUFF
!      type(DFTI_DESCRIPTOR), POINTER :: Desc_Handle_Dim1
      Integer :: Status
      INTEGER :: STRIDE(2),INFO
      INTEGER :: L(2),G(3),GSTRIDE(3)
c=================================
!#if FFTW
!       write(6,*) "'THE IMEST FOR CUDA!!!"
!#endif
!compute f_sl = dx phi_s^* phi_l
       do m=1,NDX*NDY*NDZ
         fsltilde(m)=DCONJG(psi_s(m))*psi_l(m)
       enddo

!get f_sl-tilde dx from f_sl dx
!--- {{{ 1d 
!c   compute forward fft of fsltilde
      if (DIM_MCTDHB.eq.1) then

         call  kernel_imestfft(Wsl,fsltilde, NDX,NDY,NDZ, 1, vtilde )
 

      else if (DIM_MCTDHB.eq.2) then
 
         call   kernel_imestfft(Wsl,fsltilde, NDX,NDY,NDZ, 2,vtilde)


       else if (DIM_MCTDHB.eq.3) then

         call    kernel_imestfft(Wsl,fsltilde, NDX,NDY,NDZ, 3, vtilde )

       endif

!c multiply it with the FFT of the potential
!      do m=1,NDX*NDY*NDZ
! 
!        Wsl(m)=fsltilde(m)*vtilde(m)/NDX/NDY/NDZ
!      enddo

       return
!---}}}
       end subroutine computewsl_CPPCUDAIMEST

