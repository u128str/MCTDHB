       subroutine computevtilde(vtilde,wtilde)
       use SHARED_DIMS 
       implicit none
       double complex                  :: vtilde(ndx*ndy*ndz)
       double complex                  :: wtilde(ndx*ndy*ndz)
c=========================================================
#ifdef FFTIMEST             
       call  computevtilde_ FFTIMEST(vtilde,wtilde)
#endif
       end subroutine computevtilde
