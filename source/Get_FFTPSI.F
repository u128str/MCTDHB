! Makes FFT on incoming set of the orbitals and returns the result
! modus=1 FFFT forward
! modus=2 FFFT backward
! Usually afterward to visualize the results one has to call Shift_Zero_K_FFT 
! It resorst the outcome written in  a STANDARD (MKL,FFTW) way (0,dp,..,N/2 dp, (-N/2+1)*dp,...,dp)
! to user -"adapted" way (-dp,0,dp) (To get from FFT of Gaussian Gaussian ;-))
      subroutine Get_FFTPSI(PSIIN,FFTPSI,modus)
      use SHARED_DIMS 
      implicit none
      COMPLEX*16,   DIMENSION(NDX*NDY*NDZ) :: tmpIN,tmpOUT 
      COMPLEX*16,   DIMENSION(NDX*NDY*NDZ,Morb) :: PSIIN, FFTPSI
      INTEGER :: modus,I
!============================================== Cycle over Morb ===================================
               Do I=1,Morb
                tmpIN=PSIIN(:,I)
         Dim_1_2_3D: SELECT CASE (DIM_MCTDHB)
         CASE (1)
#ifdef SFX1D              
                   CALL Get_FFTPSI_ SFX1D (tmpIN,tmpOUT,modus)  !FFT 1D in X direction
#endif
         CASE (2)
#ifdef SFX2D              
                   CALL Get_FFTPSI_ SFX2D (tmpIN,tmpOUT,modus)  !FFT 2D MKL or FFTW
#endif
         CASE (3)
#ifdef SFX3D              
                   CALL Get_FFTPSI_ SFX3D (tmpIN,tmpOUT,modus)  !TRUE 3D FFT 3D MKL FFTW or CUDA CPP 
#endif
         END SELECT Dim_1_2_3D
                FFTPSI(:,I)=tmpOUT
               EndDO
       return
       end subroutine Get_FFTPSI

