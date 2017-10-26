MODULE CI_SUBR 
    INTERFACE 
!       SUBROUTINE  GetCIJKL1body_Par(MYID,VIN) 
!        ! CALL  GetCIJKL1body_Par(MYID,VIN) 
!          COMPLEX*16 ::  VIN(:) 
!          INTEGER    :: MYID 
!       END SUBROUTINE  GetCIJKL1body_Par 
!       SUBROUTINE  GetCIJKL2body_Par(MYID,VIN) 
!           ! CALL  GetCIJKL2body_Par(MYID,VIN) 
!          COMPLEX*16 ::  VIN(:) 
!          INTEGER    :: MYID 
!       END SUBROUTINE  GetCIJKL2body_Par 
    SUBROUTINE  PrdCIJKL1body_M_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  PrdCIJKL1body_M_OMP
    SUBROUTINE  PrdCIJKL2body_M_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  PrdCIJKL2body_M_OMP
    SUBROUTINE  GetCIJKL1body_1(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_1
    SUBROUTINE  GetCIJKL2body_1(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_1
    SUBROUTINE  GetCIJKL1body_2_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_2_OMP
    SUBROUTINE  GetCIJKL2body_2_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_2_OMP
    SUBROUTINE  GetCIJKL1body_3_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_3_OMP
    SUBROUTINE  GetCIJKL2body_3_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_3_OMP
    SUBROUTINE  GetCIJKL1body_4_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_4_OMP
    SUBROUTINE  GetCIJKL2body_4_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_4_OMP
    SUBROUTINE  GetCIJKL1body_5_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_5_OMP
    SUBROUTINE  GetCIJKL2body_5_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_5_OMP
    SUBROUTINE  GetCIJKL1body_6_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_6_OMP
    SUBROUTINE  GetCIJKL2body_6_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_6_OMP
    SUBROUTINE  GetCIJKL1body_7_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_7_OMP
    SUBROUTINE  GetCIJKL2body_7_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_7_OMP
    SUBROUTINE  GetCIJKL1body_8_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_8_OMP
    SUBROUTINE  GetCIJKL2body_8_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_8_OMP
    SUBROUTINE  GetCIJKL1body_9_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_9_OMP
    SUBROUTINE  GetCIJKL2body_9_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_9_OMP
    SUBROUTINE  GetCIJKL1body_10_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_10_OMP
    SUBROUTINE  GetCIJKL2body_10_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_10_OMP
    SUBROUTINE  GetCIJKL1body_11_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_11_OMP
    SUBROUTINE  GetCIJKL2body_11_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_11_OMP
    SUBROUTINE  GetCIJKL1body_12_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_12_OMP
    SUBROUTINE  GetCIJKL2body_12_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_12_OMP
    SUBROUTINE  GetCIJKL1body_13_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_13_OMP
    SUBROUTINE  GetCIJKL2body_13_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_13_OMP
    SUBROUTINE  GetCIJKL1body_14_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_14_OMP
    SUBROUTINE  GetCIJKL2body_14_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_14_OMP
    SUBROUTINE  GetCIJKL1body_15_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_15_OMP
    SUBROUTINE  GetCIJKL2body_15_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_15_OMP
    SUBROUTINE  GetCIJKL1body_16_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_16_OMP
    SUBROUTINE  GetCIJKL2body_16_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_16_OMP
    SUBROUTINE  GetCIJKL1body_17_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_17_OMP
    SUBROUTINE  GetCIJKL2body_17_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_17_OMP
    SUBROUTINE  GetCIJKL1body_18_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_18_OMP
    SUBROUTINE  GetCIJKL2body_18_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_18_OMP
    SUBROUTINE  GetCIJKL1body_19_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_19_OMP
    SUBROUTINE  GetCIJKL2body_19_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_19_OMP
    SUBROUTINE  GetCIJKL1body_20_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_20_OMP
    SUBROUTINE  GetCIJKL2body_20_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_20_OMP
    SUBROUTINE  GetCIJKL1body_21_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_21_OMP
    SUBROUTINE  GetCIJKL2body_21_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_21_OMP
    SUBROUTINE  GetCIJKL1body_22_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_22_OMP
    SUBROUTINE  GetCIJKL2body_22_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_22_OMP
    SUBROUTINE  GetCIJKL1body_23_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_23_OMP
    SUBROUTINE  GetCIJKL2body_23_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_23_OMP
    SUBROUTINE  GetCIJKL1body_24_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_24_OMP
    SUBROUTINE  GetCIJKL2body_24_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_24_OMP
    SUBROUTINE  GetCIJKL1body_25_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_25_OMP
    SUBROUTINE  GetCIJKL2body_25_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_25_OMP
    SUBROUTINE  GetCIJKL1body_26_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_26_OMP
    SUBROUTINE  GetCIJKL2body_26_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_26_OMP
    SUBROUTINE  GetCIJKL1body_27_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_27_OMP
    SUBROUTINE  GetCIJKL2body_27_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_27_OMP
    SUBROUTINE  GetCIJKL1body_28_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_28_OMP
    SUBROUTINE  GetCIJKL2body_28_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_28_OMP
    SUBROUTINE  GetCIJKL1body_29_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_29_OMP
    SUBROUTINE  GetCIJKL2body_29_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_29_OMP
    SUBROUTINE  GetCIJKL1body_30_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_30_OMP
    SUBROUTINE  GetCIJKL2body_30_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_30_OMP
    SUBROUTINE  GetCIJKL1body_31_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_31_OMP
    SUBROUTINE  GetCIJKL2body_31_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_31_OMP
    SUBROUTINE  GetCIJKL1body_32_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_32_OMP
    SUBROUTINE  GetCIJKL2body_32_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_32_OMP
    SUBROUTINE  GetCIJKL1body_33_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_33_OMP
    SUBROUTINE  GetCIJKL2body_33_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_33_OMP
    SUBROUTINE  GetCIJKL1body_34_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_34_OMP
    SUBROUTINE  GetCIJKL2body_34_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_34_OMP
    SUBROUTINE  GetCIJKL1body_35_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_35_OMP
    SUBROUTINE  GetCIJKL2body_35_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_35_OMP
    SUBROUTINE  GetCIJKL1body_36_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_36_OMP
    SUBROUTINE  GetCIJKL2body_36_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_36_OMP
    SUBROUTINE  GetCIJKL1body_37_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_37_OMP
    SUBROUTINE  GetCIJKL2body_37_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_37_OMP
    SUBROUTINE  GetCIJKL1body_38_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_38_OMP
    SUBROUTINE  GetCIJKL2body_38_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_38_OMP
    SUBROUTINE  GetCIJKL1body_39_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_39_OMP
    SUBROUTINE  GetCIJKL2body_39_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_39_OMP
    SUBROUTINE  GetCIJKL1body_40_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_40_OMP
    SUBROUTINE  GetCIJKL2body_40_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_40_OMP
    SUBROUTINE  GetCIJKL1body_41_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_41_OMP
    SUBROUTINE  GetCIJKL2body_41_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_41_OMP
    SUBROUTINE  GetCIJKL1body_42_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_42_OMP
    SUBROUTINE  GetCIJKL2body_42_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_42_OMP
    SUBROUTINE  GetCIJKL1body_43_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_43_OMP
    SUBROUTINE  GetCIJKL2body_43_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_43_OMP
    SUBROUTINE  GetCIJKL1body_44_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_44_OMP
    SUBROUTINE  GetCIJKL2body_44_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_44_OMP
    SUBROUTINE  GetCIJKL1body_45_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_45_OMP
    SUBROUTINE  GetCIJKL2body_45_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_45_OMP
    SUBROUTINE  GetCIJKL1body_46_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_46_OMP
    SUBROUTINE  GetCIJKL2body_46_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_46_OMP
    SUBROUTINE  GetCIJKL1body_47_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_47_OMP
    SUBROUTINE  GetCIJKL2body_47_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_47_OMP
    SUBROUTINE  GetCIJKL1body_48_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_48_OMP
    SUBROUTINE  GetCIJKL2body_48_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_48_OMP
    SUBROUTINE  GetCIJKL1body_49_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_49_OMP
    SUBROUTINE  GetCIJKL2body_49_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_49_OMP
    SUBROUTINE  GetCIJKL1body_50_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL1body_50_OMP
    SUBROUTINE  GetCIJKL2body_50_OMP(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
    END SUBROUTINE  GetCIJKL2body_50_OMP
    END INTERFACE
END MODULE CI_SUBR
