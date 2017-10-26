MODULE PASS_ARG
INTERFACE attI
    SUBROUTINE Guess_PSI(PSI)
    USE SHARED_DIMS
    COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb) :: Psi
    END SUBROUTINE Guess_PSI

    SUBROUTINE Guess_CI(VIN)
        ! CALL Guess_CI(PSI)
       COMPLEX*16 ::  VIN(:)
    END SUBROUTINE Guess_CI

    SUBROUTINE Get_h_W(PSI,time)
    USE SHARED_DIMS
!       COMPLEX*16 ::  PSI(:,:)
    COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb) :: Psi
       real*8 :: time
    END SUBROUTINE Get_h_W

    SUBROUTINE Get_r_R(VIN,VOUT,time)
        ! CALL Get_r_R(VIN,VOUT,time)
       COMPLEX*16 ::  VIN(:),VOUT(:)
       real*8 :: time
    END SUBROUTINE Get_r_R

!    SUBROUTINE HPSI(VIN,VOUT)
!        ! CALL HPSI(VOUT1,VOUT)
!       COMPLEX*16 ::  VIN(:),VOUT(:)
!    END SUBROUTINE HPSI

    SUBROUTINE SIL_PROPG(time,VIN,VOUT,Error_SIL,MAXIT)
       COMPLEX*16 ::  VIN(:),VOUT(:)
       real*8 :: time,Error_SIL
       integer :: MAXIT
    END SUBROUTINE  SIL_PROPG

    SUBROUTINE     Get_OPSI_WSL_balanced_OMP (PSI,OPSI, FromN,TillN,WSLTime)
        ! CALL Get_OPSI_WSL_based_OMP (PSI,OPSI,FromN,TillN,WSLTime) !AIS 15JUL2014
    USE SHARED_DIMS
      REAL*8  :: WSLTime
    COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb) :: Psi,OPSI
        integer :: FromN,TillN
    END SUBROUTINE Get_OPSI_WSL_balanced_OMP

!    SUBROUTINE     Get_OPSI_WSL_disbalanced_OMP (PSI,OPSI,FromN,TillN,WSLTime)
!    REAL*8  :: WSLTime
!    USE SHARED_DIMS
!    COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb) :: Psi,OPSI
!    integer :: FromN,TillN
!    END SUBROUTINE Get_OPSI_WSL_disbalanced_OMP

END INTERFACE attI
END MODULE PASS_ARG
