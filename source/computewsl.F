      subroutine computewsl(wsl,psi_s,psi_l,vtilde)
      use SHARED_DIMS 
      implicit none
      double complex     :: Wsl(NDX*NDY*NDZ)
      double complex  :: psi_s(NDX*NDY*NDZ)
      double complex  :: psi_l(NDX*NDY*NDZ)
      double complex  :: vtilde(NDX*NDY*NDZ)

#ifdef FFTIMEST             
       call  computewsl_ FFTIMEST(wsl,psi_s,psi_l,vtilde)
#endif
       return
       end subroutine computewsl

