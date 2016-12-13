      subroutine  GetCIJKL_All_body_Par(MYID,VIN)
      use PASS_ARG
      USE CI_SUBR
      USE SHARED_DIMS
      USE CI_All 
      USE Parallel_CI
      USE rR_hW
      implicit NONE
       INTERFACE 
       SUBROUTINE  GetCIJKL1body_Par(MYID,VIN)
       ! CALL  GetCIJKL1body_Par(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
       END SUBROUTINE  GetCIJKL1body_Par
       SUBROUTINE  GetCIJKL2body_Par(MYID,VIN)
        ! CALL  GetCIJKL2body_Par(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
       END SUBROUTINE  GetCIJKL2body_Par
       END INTERFACE 
!c=========================================================
      INTEGER ::  MYID,Ist,Ifn,TRM
      COMPLEX*16 :: VIN(:)
      COMPLEX*16, ALLOCATABLE  :: VIN1(:)
!      external GetCIJKL1body_Par,GetCIJKL2body_Par
!c=========================================================
!      write(6,*)MYID,"GetCIJKL_All_body_Par 0",SUM(VIN)  

           TRM=MYID_TRM(MYID+1)
           TRM_choice: SELECT CASE (TRM)
                    CASE (1)
                    CALL GetCIJKL1body_Par(MYID,VIN) ! Worker on 1b
                    CASE (2)
                    CALL GetCIJKL2body_Par(MYID,VIN) ! Worker on 2b
                    CASE (3)
        Ist=CI_Proc_From(EXCEPTIONAL_ID+1) 
        Ifn=CI_Proc_Till(EXCEPTIONAL_ID+1) 
        CI_Proc_From(EXCEPTIONAL_ID+1)=Ist
        CI_Proc_Till(EXCEPTIONAL_ID+1)=MaxTrm1b
       ALLOCATE( VIN1(SIZE(VIN)) )
       IF(ALLOCATED(VIN1).eqv..FALSE.) &
      write(6,*)" MEM for VIN1 is NOT ok!",SIZE(VIN1)
        VIN1=VIN
        CALL GetCIJKL1body_Par(MYID,VIN)
        CI_Proc_From(EXCEPTIONAL_ID+1)=1
        CI_Proc_Till(EXCEPTIONAL_ID+1)=Ifn
        CALL GetCIJKL2body_Par(MYID,VIN1)
        CI_Proc_From(EXCEPTIONAL_ID+1)=Ist
        CI_Proc_Till(EXCEPTIONAL_ID+1)=Ifn
        VIN=VIN+VIN1
        DEALLOCATE(VIN1)
                    CASE DEFAULT 
       write(6,*)" Something wrong in GetCIJKL_All_body_Par!!!!!!!!!"
                    END SELECT TRM_choice

!      write(6,*)"GetCIJKL_All_body_Par 1",SUM(VIN)  
      end subroutine GetCIJKL_All_body_Par
