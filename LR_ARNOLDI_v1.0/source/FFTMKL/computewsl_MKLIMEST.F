      subroutine computewsl_MKLIMEST(wsl,psi_s,psi_l,vtilde)
      use SHARED_DIMS 
      use DVR_ALL
      Use MKL_DFTI
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
      type(DFTI_DESCRIPTOR), POINTER :: Desc_Handle_Dim1
      Integer :: Status
      INTEGER :: STRIDE(2),INFO
      INTEGER :: L(2),G(3),GSTRIDE(3)
c=================================
#if FFTW
       write(6,*) 'THE IMEST ALGORITHM IS NOT IMPLEMENTED FOR FFTW!!!'
       STOP
#endif
!compute f_sl = dx phi_s^* phi_l
       do m=1,NDX*NDY*NDZ
         fsltilde(m)=DCONJG(psi_s(m))*psi_l(m)
       enddo

!get f_sl-tilde dx from f_sl dx
!--- {{{ 1d 
c   compute forward fft of fsltilde
      if (DIM_MCTDHB.eq.1) then
        Status = DftiCreateDescriptor(Desc_Handle_Dim1, DFTI_DOUBLE,
     &      DFTI_COMPLEX, 1, NDX )
        Status = DftiCommitDescriptor(Desc_Handle_Dim1)
        Status = DftiComputeForward(Desc_Handle_Dim1,fsltilde)
      else if (DIM_MCTDHB.eq.2) then
        L(1)=NDX
        L(2)=NDY
        Status = DftiCreateDescriptor(Desc_Handle_Dim1, DFTI_DOUBLE,
     &    DFTI_COMPLEX, 2, L )
        Status = DftiCommitDescriptor( Desc_Handle_Dim1) 
        Status = DftiComputeForward(Desc_Handle_Dim1,fsltilde)
      else if (DIM_MCTDHB.eq.3) then
        G(1)=NDX
        G(2)=NDY
        G(3)=NDZ
        Status = DftiCreateDescriptor(Desc_Handle_Dim1, DFTI_DOUBLE,
     &    DFTI_COMPLEX, 3, G )
        Status = DftiCommitDescriptor( Desc_Handle_Dim1) 
        Status = DftiComputeForward(Desc_Handle_Dim1,fsltilde)
      endif

c multiply it with the FFT of the potential
      do m=1,NDX*NDY*NDZ
         Wsl(m)=fsltilde(m)*vtilde(m)/NDX/NDY/NDZ
      enddo

c       transform back 
      if (DIM_MCTDHB.eq.1) then
         Status = DftiComputeBackward(Desc_Handle_Dim1,wsl)
         Status = DftiFreeDescriptor(Desc_Handle_Dim1)
      else if (DIM_MCTDHB.eq.2) then
         Status = DftiComputeBackward(Desc_Handle_Dim1,wsl)
         Status = DftiFreeDescriptor(Desc_Handle_Dim1)
      else if (DIM_MCTDHB.eq.3) then
         Status = DftiComputeBackward(Desc_Handle_Dim1,wsl)
         Status = DftiFreeDescriptor(Desc_Handle_Dim1)
      endif


       return
!---}}}
       end subroutine computewsl_MKLIMEST

