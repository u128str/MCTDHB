        SUBROUTINE Share_H_W
        USE SHARED_DIMS
        USE rR_hW
        USE CI_ALL
        implicit NONE
c================================== MPI ======================================
       INCLUDE 'mpif.h'
       INTEGER ::  ierr,MYID,numprocs
c=============================================================================
          INTEGER :: n
c=============================================================================
       call MPI_COMM_RANK(MPI_COMM_WORLD,MYID,ierr)
       call MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,ierr)
c       write(6,*)"Get_r_R Process ",MYID," of ", numprocs
c=============================================================================
c========================= MPI PART ======================================
       n=Morb*(Morb+1)/2
       call MPI_BCAST(HIJ,n,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
c       call MPI_BCAST(MaxTrm1b,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c       call MPI_BCAST(TERM_REQ_1B,MaxTrm1b,MPI_DOUBLE_COMPLEX,0,
c     .                                           MPI_COMM_WORLD,ierr)
       n=Morb*(Morb+1)/2*(Morb*(Morb+1)/2+1)/2
       call MPI_BCAST(WIJKL,n,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
c       call MPI_BCAST(MaxTrm2b,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c       call MPI_BCAST(TERM_REQ_2B,MaxTrm2b,MPI_DOUBLE_COMPLEX,0,
c     .                                           MPI_COMM_WORLD,ierr)

        END SUBROUTINE Share_H_W
