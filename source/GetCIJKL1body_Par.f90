      subroutine  GetCIJKL1body_Par(MYID,VIN) 
      USE PASS_ARG  
      USE CI_SUBR 
      USE SHARED_DIMS 
      USE CI_All 
      USE CI_Prod 
      implicit NONE 
!c========================================================= 
      INTEGER ::  MYID 
!c      COMPLEX*16, DIMENSION(Nconf) :: VIN 
      COMPLEX*16 :: VIN(:) 
!c========================================================= 
!     write(6,*)MYID,"GetCIJKL1body_Par 0" 
        MorbChoice: SELECT CASE (Morb) 
                    CASE (1) 
                    CALL GetCIJKL1body_1(MYID,VIN) 
                    CASE (2) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_2_OMP(MYID,VIN) 
                    CASE (3) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_3_OMP(MYID,VIN) 
                    CASE (4) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_4_OMP(MYID,VIN) 
                    CASE (5) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_5_OMP(MYID,VIN) 
                    CASE (6) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_6_OMP(MYID,VIN) 
                    CASE (7) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_7_OMP(MYID,VIN) 
                    CASE (8) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_8_OMP(MYID,VIN) 
                    CASE (9) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_9_OMP(MYID,VIN) 
                    CASE (10) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_10_OMP(MYID,VIN) 
                    CASE (11) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_11_OMP(MYID,VIN) 
                    CASE (12) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_12_OMP(MYID,VIN) 
                    CASE (13) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_13_OMP(MYID,VIN) 
                    CASE (14) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_14_OMP(MYID,VIN) 
                    CASE (15) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_15_OMP(MYID,VIN) 
                    CASE (16) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_16_OMP(MYID,VIN) 
                    CASE (17) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_17_OMP(MYID,VIN) 
                    CASE (18) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_18_OMP(MYID,VIN) 
                    CASE (19) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_19_OMP(MYID,VIN) 
                    CASE (20) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_20_OMP(MYID,VIN) 
                    CASE (21) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_21_OMP(MYID,VIN) 
                    CASE (22) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_22_OMP(MYID,VIN) 
                    CASE (23) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_23_OMP(MYID,VIN) 
                    CASE (24) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_24_OMP(MYID,VIN) 
                    CASE (25) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_25_OMP(MYID,VIN) 
                    CASE (26) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_26_OMP(MYID,VIN) 
                    CASE (27) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_27_OMP(MYID,VIN) 
                    CASE (28) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_28_OMP(MYID,VIN) 
                    CASE (29) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_29_OMP(MYID,VIN) 
                    CASE (30) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_30_OMP(MYID,VIN) 
                          CASE (31:100)
            write(6,*)" Still NOT generated!!!!!!!!!" 
                          CASE DEFAULT 
            write(6,*)" Something wrong in Number of orbitals!!!!!!!!!" 
                          END SELECT MorbChoice 
      !     write(6,*) MYID,"GetCIJKL1body_Par 1",SUM(VIN) 
            end subroutine GetCIJKL1body_Par 
