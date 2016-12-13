! After   FFT  one needs to Shift zero-frequency component to center of spectrum for visualization
! by default FFT returns the result corresponding to this order of momenta (0,dp,..,N/2 dp, (-N/2+1)*dp,...,dp)
! so  we have to reort the outcome accordingly and also normalize it !!!!
! see also http://www.mathworks.de/de/help/matlab/ref/fftshift.html
      SUBROUTINE Get_Shift_Zero_K_FFT(psiin) 
      USE SHARED_DIMS
!      USE omp_lib
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER ::  ierr
      COMPLEX*16,   DIMENSION(NDX*NDY*NDZ,Morb) :: psiin
!      COMPLEX*16,   DIMENSION(NDX*NDY*NDZ) :: tmp
      COMPLEX*16,   DIMENSION(NDX,NDY,NDZ) :: retemp
      INTEGER :: m,I
      REAL*8 :: norm
c===================================================
c===================================================
                       Do I=1,Morb
         retemp=reshape( psiin(:,I), (/ NDX,NDY,NDZ /) )
         Dim_1_2_3D: SELECT CASE (DIM_MCTDHB)
         CASE (1) ! 1D  ! 1D  ! 1D  ! 1D  ! 1D  ! 1D  ! 1D  ! 1D  ! 1D 
!         if (mod(NDX,2).eq.0) then
!           do m=(NDX/2+1),NDX
!              tmp(m-(NDX/2))=psiin(m,I)
!              tmp(m)=psiin(m-(NDX/2),I)
!           enddo
!         else
!From Axel -- wrong!!!!         do m=(NDX/2+1),NDX
!              tmp(m-(NDX/2))=psiin(m,I)
!              tmp(m)=psiin(m-(NDX/2),I)
!           enddo
!         endif
!ORG OK      if (mod(NDX,2).eq.0) tmp = cshift(psiin(:,I), SHIFT=NDX/2)
!ORG OK      if (mod(NDX,2).ne.0) tmp = cshift(psiin(:,I), SHIFT=NDX/2+1)
      if (mod(NDX,2).eq.0) retemp = cshift(retemp, SHIFT=NDX/2,DIM=1)
      if (mod(NDX,2).ne.0) retemp = cshift(retemp, SHIFT=NDX/2+1,DIM=1)
         CASE (2) ! 2D  2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D
      if (mod(NDX,2).eq.0) retemp = cshift(retemp, SHIFT=NDX/2,DIM=1)
      if (mod(NDX,2).ne.0) retemp = cshift(retemp, SHIFT=NDX/2+1,DIM=1)
      if (mod(NDY,2).eq.0) retemp = cshift(retemp, SHIFT=NDY/2,DIM=2)
      if (mod(NDY,2).ne.0) retemp = cshift(retemp, SHIFT=NDY/2+1,DIM=2)
         CASE (3) ! 3D ! 3D ! 3D ! 3D ! 3D ! 3D ! 3D ! 3D ! 3D ! 3D ! 3D
      if (mod(NDX,2).eq.0) retemp = cshift(retemp, SHIFT=NDX/2,DIM=1)
      if (mod(NDX,2).ne.0) retemp = cshift(retemp, SHIFT=NDX/2+1,DIM=1)
      if (mod(NDY,2).eq.0) retemp = cshift(retemp, SHIFT=NDY/2,DIM=2)
      if (mod(NDY,2).ne.0) retemp = cshift(retemp, SHIFT=NDY/2+1,DIM=2)
      if (mod(NDZ,2).eq.0) retemp = cshift(retemp, SHIFT=NDZ/2,DIM=3)
      if (mod(NDZ,2).ne.0) retemp = cshift(retemp, SHIFT=NDZ/2+1,DIM=3)

         END SELECT Dim_1_2_3D
!=================Normalization============================
                psiin(:,I)=reshape(retemp , (/ NDX*NDY*NDZ /) )
         call normvxz(psiin(:,I),norm,NDX*NDY*NDZ) !Renormalization of the FFT results
         call xvixdzo(norm,psiin(:,I),NDX*NDY*NDZ)
                        End Do ! End cycle over all orbitals
c===============================================================================
c================================================================================
      END  SUBROUTINE Get_Shift_Zero_K_FFT
