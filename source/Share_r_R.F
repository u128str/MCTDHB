       subroutine Share_r_R
       USE SHARED_DIMS
       USE rR_hW
       USE CI_All
       USE Parallel_Orb
       implicit NONE
c=================== MPI
       INCLUDE 'mpif.h'
       INTEGER ::  ierr,MYID,numprocs,i
       integer n
c==========================================================
         call MPI_COMM_RANK(MPI_COMM_WORLD,MYID,ierr)
         call MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,ierr)
c         write(6,*)"Share_r_r Process ",MYID," of ", numprocs
c==========================================================
c         write(6,*)"Share_r_r ZRIJ ",MYID,ZRIJ(1)
c==========================================================================
       n=Morb*(Morb+1)/2
       call MPI_BCAST(ZRIJ,n,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
c       call MPI_BCAST(MaxTrm1b,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c       call MPI_BCAST(TERM_REQ_1B,MaxTrm1b,MPI_DOUBLE_COMPLEX,0,
c     .                                           MPI_COMM_WORLD,ierr)
c         write(6,*)"Share_r_r ZRIJKL ",MYID,ZRIJKL(1)
       n=Morb*(Morb+1)/2*(Morb*(Morb+1)/2+1)/2
       call MPI_BCAST(ZRIJKL,n,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
c       call MPI_BCAST(MaxTrm2b,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c       call MPI_BCAST(TERM_REQ_2B,MaxTrm2b,MPI_DOUBLE_COMPLEX,0,
c     .                                           MPI_COMM_WORLD,ierr)
c         write(6,*)"Share_r_r AllZrij",MYID,AllZRIJ(1,1)
      n=Morb*Morb
      call MPI_BCAST(InvZRIJ,n,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(AllZRIJ,n,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
c==========================================================
      end subroutine Share_r_R

