! Added by STR Feb 6 (2015)
! Subroutine to provide dx dy or dz*Psi -> psi Op - is operator Matrix in the
! respective DVR presentation DIM(Op)=(NDX,NDX) in the FFT case it is treated differently
! X_i index pointing to the T_x,T_y or T_z will be calculated - only for
! 2D and 3D problems
! NdimMax = 1D(NDX)
! NdimMax = 2D(NDX*NDY)
! NdimMax = 3D(NDX*NDY*NDZ)
          SUBROUTINE Get_d_PSI(psi,X_i) 
          USE SHARED_DIMS
          USE DVR_ALL
          IMPLICIT NONE
          INTEGER                     :: X_i
          COMPLEX*16,   DIMENSION(NDX*NDY*NDZ) :: psi
c==================== Local Arrays for 3D Transformations 
          INTEGER :: I
c           Cop=ZONER*Op
c===============================================================================
            Xi1: SELECT CASE (X_i)
         CASE (1)
             IF(Time_DVRMETHODX.ne.4) THEN
         call    Get_Op_PSI_Gnrl(psi,dif1matX,NDX,1)
             ELSE
#ifdef SFX1D              
             CALL Get_d_PSI_ SFX1D (psi,1)  !FFT 1D in X direction
#endif
             ENDIF
         CASE (2)
             IF(Time_DVRMETHODY.ne.4) THEN
           call  Get_Op_PSI_Gnrl(psi,dif1matY,NDY,2)
             ELSE
#ifdef SFX1D              
             CALL Get_d_PSI_ SFX1D (psi,2)  !FFT 1D in Y direction
#endif
             ENDIF
         CASE (3)
             IF(Time_DVRMETHODZ.ne.4) THEN
            call Get_Op_PSI_Gnrl(psi,dif1matZ,NDZ,3)
             ELSE
#ifdef SFX1D              
             CALL Get_d_PSI_ SFX1D (psi,3)  !FFT 1D in Z direction
#endif
             ENDIF
         END SELECT Xi1
c================================================================================
          END  SUBROUTINE Get_d_PSI
