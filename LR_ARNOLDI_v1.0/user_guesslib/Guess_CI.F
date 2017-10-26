C=============== Assistent MPI subroutive -assists in H|psi> multiplication

       subroutine Guess_CI(VIN)
       USE   SHARED_DIMS
c       USE   rR_hW
       USE   CI_All

       implicit NONE
c=================== MPI 
C       INCLUDE 'mpif.h'
C       INTEGER ::  ierr,MYID,numprocs
c==========================================================
        INTEGER :: I
c        COMPLEX*16, DIMENSION(Nconf) :: VIN 
        COMPLEX*16 :: VIN(:)
c====================== For SIL
        COMPLEX*16 :: zrho=ZERO,Z,Z1
c==============================================
      real*4   start,finish, exec_time , cnk_time,time   
      real*8  xnrm

      real*8   FKTR,CNK,DZASUM,DZNRM2,DNRM2,DDOT
      COMPLEX*16 :: ZDOTC
      external FKTR,CNK,GetCIJKL2body_Par
      external DZASUM,DZNRM2,DNRM2,ZDOTC,DDOT,ZSCAL
c=============================================================
       IF(Morb==1) VIN(1)=Zoner
c====================
       VIN=ZERO
       VIN(1)=1.0d0
! Here I check an existance of the configurations from V_W_Psi_string.in file and 
! refill the VIN array ATTENTION VIN is complex while UserCnfW is real*8!!!!!!!!!
        IF(sum(UserCnfN).ne.0) then
        VIN=ZERO
        Do i=1,min(MaxUserCnf,Nconf)
        if (UserCnfN(i).gt.0) VIN(UserCnfN(i))=UserCnfW(i) 
        endDO
        endif
! Reading from ********conf.dat has a priority over V_W_Psi_string.in !!!!
       IF(CIc_Rest) call Guess_Read_CI(VIN)

        xnrm=DZNRM2(Nconf,VIN,1) 
        VIN=VIN/xnrm
        
c==============================================
        write(6,*)"IN CI Vec NORM :",xnrm
      end subroutine Guess_CI
