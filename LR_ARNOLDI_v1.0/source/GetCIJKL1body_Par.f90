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
                    CASE (31) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_31_OMP(MYID,VIN) 
                    CASE (32) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_32_OMP(MYID,VIN) 
                    CASE (33) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_33_OMP(MYID,VIN) 
                    CASE (34) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_34_OMP(MYID,VIN) 
                    CASE (35) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_35_OMP(MYID,VIN) 
                    CASE (36) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_36_OMP(MYID,VIN) 
                    CASE (37) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_37_OMP(MYID,VIN) 
                    CASE (38) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_38_OMP(MYID,VIN) 
                    CASE (39) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_39_OMP(MYID,VIN) 
                    CASE (40) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_40_OMP(MYID,VIN) 
                    CASE (41) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_41_OMP(MYID,VIN) 
                    CASE (42) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_42_OMP(MYID,VIN) 
                    CASE (43) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_43_OMP(MYID,VIN) 
                    CASE (44) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_44_OMP(MYID,VIN) 
                    CASE (45) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_45_OMP(MYID,VIN) 
                    CASE (46) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_46_OMP(MYID,VIN) 
                    CASE (47) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_47_OMP(MYID,VIN) 
                    CASE (48) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_48_OMP(MYID,VIN) 
                    CASE (49) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_49_OMP(MYID,VIN) 
                    CASE (50) 
      IF(CI_Production_1b.eqv..TRUE.) CALL PrdCIJKL1body_M_OMP(MYID,VIN)
      IF(CI_Production_1b.eqv..FALSE.) &
                       CALL GetCIJKL1body_50_OMP(MYID,VIN) 
                          CASE (51:100)
            write(6,*)" Still NOT generated!!!!!!!!!" 
                          CASE DEFAULT 
            write(6,*)" Something wrong in Number of orbitals!!!!!!!!!" 
                          END SELECT MorbChoice 
      !     write(6,*) MYID,"GetCIJKL1body_Par 1",SUM(VIN) 
            end subroutine GetCIJKL1body_Par 
