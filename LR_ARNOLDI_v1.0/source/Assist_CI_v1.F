C=============== Assistent MPI subroutive -assists in H|psi>  multiplication
       subroutine Assist_CI
       USE PASS_ARG
!       USE CI_SUBR
       USE SHARED_DIMS
       USE rR_hW 
       USE CI_ALL
       implicit NONE
       INTERFACE
       SUBROUTINE  GetCIJKL_All_body_Par(MYID,VIN)
       ! CALL  GetCIJKL1body_Par(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
       END SUBROUTINE  GetCIJKL_All_body_Par
       END INTERFACE
c=================== MPI 
       INCLUDE 'mpif.h'
       INTEGER ::  ierr,MYID,numprocs
c==========================================================
       INTEGER ::  n,nn
cMaxTerm,Nterms,FromN,TillN
c==========================================================

      COMPLEX*16, ALLOCATABLE  :: VIN(:)
      COMPLEX*16, ALLOCATABLE  :: VOUT(:)
      Logical SIL
c      external CNK,GetCIJKL2body_Par
c================ MPI ===================================
c         call MPI_INIT(ierr)
         call MPI_COMM_RANK(MPI_COMM_WORLD,MYID,ierr)
         call MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,ierr)

c         write(6,*)"Assistent Process ",MYID," of ", numprocs
c         write(6,*)"Assistent Process ",Nconf
c===================== MEMORY allocation for VIN and VOUT
c         write(6,*)"Asist CI  IN Process ",MYID," JOB -- SIL ",SIL
        ALLOCATE(VIN(Nconf))
        IF(ALLOCATED(VIN).eqv..FALSE.)
     . write(6,*)" MEM for VIN in Proc=",MYID," is NOT ok!",SIZE(VIN)
        ALLOCATE(VOUT(Nconf))
        IF(ALLOCATED(VOUT).eqv..FALSE.)
     . write(6,*)" MEM for VOUT in Proc=",MYID," is NOT ok!",SIZE(VOUT)
c===================== MEMORY allocation for VIN and VOUT is OK !!!!!!!
c===================== Receiving  HIJ and WIJKL
c===================== Redistribute HIJ and WIJKL to all the process ====
c       n=Morb*(Morb+1)/2
c       call MPI_BCAST(HIJ,n,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
c       call MPI_BCAST(MaxTrm1b,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c       call MPI_BCAST(TERM_REQ_1B,MaxTrm1b,MPI_DOUBLE_COMPLEX,0,
c     .                                           MPI_COMM_WORLD,ierr)
c       n=Morb*(Morb+1)/2*(Morb*(Morb+1)/2+1)/2
c       call MPI_BCAST(WIJKL,n,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
c       call MPI_BCAST(MaxTrm2b,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c       call MPI_BCAST(TERM_REQ_2B,MaxTrm2b,MPI_DOUBLE_COMPLEX,0,
c     .                                           MPI_COMM_WORLD,ierr)
         call Share_H_W
c=============================================================================0
1       continue 
      call MPI_BCAST(SIL,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
        IF(SIL.eqv..TRUE.)  THEN 
             GOTO 2 
       ELSE
             GOTO 9999
       ENDIF
c==================== H|psi> PART===============================================
2       CONTINUE
      VOUT=ZERO
      VIN=ZERO
      ZRIJKL =Zero
      ZRIJKL1=Zero
      ZRIJ1 =Zero
      ZRIJ =Zero
      call MPI_BCAST(VIN,Nconf,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
c      CALL GetCIJKL2body_Par(MYID,VIN)
      CALL GetCIJKL_All_body_Par(MYID,VIN)
c================ Gathering together MPI_REDUCE ================================
      CALL MPI_REDUCE(VIN,VOUT,Nconf,MPI_DOUBLE_COMPLEX,MPI_SUM,0,
     .                                          MPI_COMM_WORLD,ierr)
c=============== Reduction of the one- and two-particle densities!!!!!!!!!!!!!
       nn=Rdim*(Rdim+1)/2  
      CALL MPI_REDUCE(ZRIJKL,ZRIJKL1,nn,MPI_DOUBLE_COMPLEX,MPI_SUM,0,
     .                                          MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(ZRIJ,ZRIJ1,Rdim,MPI_DOUBLE_COMPLEX,MPI_SUM,0,
     .                                          MPI_COMM_WORLD,ierr)
c      CALL MPI_Barrier(MPI_COMM_WORLD);
       goto 1
c================= END MPI
9999     continue
       DEALLOCATE(VIN)
       DEALLOCATE(VOUT)
c       write(6,*)"CI Assistent ",MYID," is DONE"
c        return
      end subroutine Assist_CI
