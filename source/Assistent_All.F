C=============== Assistent MPI subroutive -assists in H|psi>  multiplication
       subroutine Assist_ALL
       USE SHARED_DIMS
       USE rR_hW
!       USE CI_ALL
       implicit NONE
c=================== MPI 
       INCLUDE 'mpif.h'
       INTEGER ::  ierr,MYID,numprocs
c==========================================================
c==========================================================
c================ MPI ===================================
c         call MPI_INIT(ierr)
         call MPI_COMM_RANK(MPI_COMM_WORLD,MYID,ierr)
         call MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,ierr)
         CALL MNGR_CI_Part(NUMPROCS)
         CALL MNGR_Orb_Part(NUMPROCS)
c        CALL MNGR_Orb_Part(1)
c         write(6,*)"Assistent Process ",MYID," of ", numprocs
c         write(6,*)"Assistent Process ",Nconf
c=========================================================================
1           continue
c=========================  JOB= 0 - CI      part    =====================
c=========================  JOB= 1 - SCF ORB part    =====================
c=========================  JOB= 3 - CLOSE ALL JOBS  =====================
c         write(6,*)MYID,"Assistent Process 0 CI_SCF=",CI_SCF
         call MPI_BCAST(CI_SCF,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c         write(6,*)MYID,"Assistent Process 1 CI_SCF=",CI_SCF
c         call   MPI_Barrier
           CIvsSCF: SELECT CASE (CI_SCF)
                    CASE (0)
c       write(6,*)MYID,"-ID IN CI of CI_SCF=",CI_SCF
                    CALL Assist_CI
                               GOTO 1
                    CASE (1)
c       write(6,*)MYID,"-ID IN SCF of CI_SCF=",CI_SCF
                    CALL Assist_ORB
                              GOTO 1
                    CASE (2:)
                    GOTO 9999
                    CASE DEFAULT 
       write(6,*)" Something wrong in CI_SCF!!!!!!!!!"
                    GOTO 9999
                    END SELECT CIvsSCF

9999        continue
c===================== MEMORY allocation for VIN and VOUT
       write(6,*)" Assistent ",MYID," is DONE"
        return
      end subroutine Assist_ALL
