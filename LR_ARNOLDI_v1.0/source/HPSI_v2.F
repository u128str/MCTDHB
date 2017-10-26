! V1 - call to CALL GetCIJKL1body_Par(MYID,VOUT1) is installed between MPI_BCAST and MPI_REDUCE 
! this allows to make more flexible menagement of the parallel part - it commes togather with 
! modification of Parallel_MNGR_CI_v1
! V2 -parallization of 1b block
C=============== Short Iterative Lanzcos Vin(o) -> Vout(time)
       subroutine HPSI(VIN,VOUT)
       use PASS_ARG
       USE CI_SUBR
       USE SHARED_DIMS
       USE rR_hW
       USE CI_ALL
       implicit NONE
      INTERFACE 
      SUBROUTINE GetCIJKL_All_body_Par(MYID,VIN)
        ! CALL GetCIJKL_All_body_Par(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
      END SUBROUTINE GetCIJKL_All_body_Par
      END INTERFACE 
c=================== MPI 
       INCLUDE 'mpif.h'
       INTEGER ::  ierr,MYID,numprocs
c==========================================================
       integer i,j,k,l,P,nn

       integer  maxsil,NL,Iter_Total,Nterms
       INTEGER ::  n,FromN,TillN
c=================== F90 F95
c        INTEGER, DIMENSION(:,:), ALLOCATABLE :: MCNK

c        COMPLEX*16, DIMENSION(Nconf) :: VIN 
c        COMPLEX*16, DIMENSION(Nconf) :: VOUT
        COMPLEX*16 :: VIN(:),VOUT(:)
c        COMPLEX*16, DIMENSION(:), ALLOCATABLE :: VOUT1 !,VIN1
c====================== For SIL
        REAL*8 :: time, Error_SIL, E_state

        COMPLEX*16 :: zrho=ZERO,Z,Z1
        COMPLEX*16 :: Escale
c==============================================
      real*4 start,finish,exec_time ,finish1, iter _time       
      REAL*8 :: start_mpi,start_mpi_all
      REAL*8 :: mpi_time_bc,mpi_time_cp,mpi_time_rd

      LOGICAL  CNV,SIL,DNS

c      external GetCIJKL1body_Par,GetCIJKL2body_Par
c      a ,GetCIJKL_All_body_Par
c===========================================================
c            write(6,*)"HPSI -1",SUM(VIN)
            start_mpi_all=MPI_WTIME(ierr)
         call MPI_COMM_RANK(MPI_COMM_WORLD,MYID,ierr)
         call MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,ierr)
c===========================================================
c        write(6,*)"H  =",(ZHJK(TRM_REQ1B(i)),i=1,MaxTrm1b)
c        write(6,*)"W  =",(ZWIJKL(TRM_REQ2B(i)),i=1,MaxTrm2b)
c        write(6,*)"Vin",VIN(407)

c       ALLOCATE(VOUT1(Nconf))
c       ALLOCATE(VIN1(Nconf))
c       IF(ALLOCATED(VOUT1).eqv..FALSE.)
c     . write(6,*) " MEM for VOUT1  is NOT ok!", SIZE(VOUT1)
c==========================================================
c          call cpu_time(start)
             start_mpi=MPI_WTIME(ierr)
c=========================================================
               SIL=.TRUE.
      call MPI_BCAST(SIL,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
c=========================================================

c              VOUT=VIN
c              VOUT1=VIN
c              VIN1=VIN
              VOUT   =ZERO
              ZRIJKL =Zero
              ZRIJKL1=Zero
              ZRIJ   =Zero 
              ZRIJ1  =Zero
c        write(6,*)"TRM_REQ1B  =",(TRM_REQ1B(i),i=1,MaxTrm1b)
!        CALL GetCIJKL1body_I(Npar,Morb,Nconf,VOUT1, ZHJK,ZRJK,TRM_REQ1B,MaxTrm1b)

c       CALL GetCIJKL1body_Par(MYID,VOUT1)
c        write(6,*)" Norm 1b PSI",DZNRM2(Nconf,VOUT1,1)
c        write(6,*)"TRM_REQ1B  =",(TRM_REQ1B(i),i=1,MaxTrm1b)
c==================================================================
C      n=3000
c      call MPI_BCAST(ZWIJKL,n,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
c      call MPI_BCAST(TRM_REQ2B,n,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c      call MPI_BCAST(MaxTrm2b,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c          cal cpu_time(start)
c            start_mpi=MPI_WTIME(ierr)
      call MPI_BCAST(VIN,Nconf,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
c          write(6,*)"SEND VIN",MPI_WTIME(ierr)-start_mpi

c      Nterms=MaxTrm2b/(NUMPROCS)
c      FromN=MYID*Nterms+1
c      TillN=   FromN+Nterms-1
c      IF(TillN.ge.MaxTrm2b) TillN=MaxTrm2b
c      write(6,*)"Master",MYID," works from",FromN," till",TillN,"!"
c==================================================================
c      CALL GetCIJKL1body_Par(MYID,VOUT1)
cfor V1      IF(NUMPROCS.eq.1) CALL GetCIJKL2body_Par(MYID,VIN)
c      CALL GetCIJKL2body_Par(MYID,VIN)
            start_mpi=MPI_WTIME(ierr)
c            write(6,*)"HPSI 0",SUM(VIN)
      CALL GetCIJKL_All_body_Par(MYID,VIN)
c            write(6,*)MYID,":MYID HPSI 1",SUM(VIN)
c          write(6,*)"COMP TRMs",MPI_WTIME(ierr)-start_mpi
c==================================================================
c      write(6,*)"Master Npar,Morb,FromN,TillN",Npar,Morb,FromN,TillN
c================ Gathering togeather MPI_REDUCE
            start_mpi=MPI_WTIME(ierr)
c            write(6,*)"HPSI 1.4",SUM(VIN),SUM(VOUT),Nconf
      CALL MPI_REDUCE(VIN,VOUT,Nconf,MPI_DOUBLE_COMPLEX,MPI_SUM,0,
     .                                          MPI_COMM_WORLD,ierr)
c          write(6,*)"REDUCE VOUT",MPI_WTIME(ierr)-start_mpi
c            write(6,*)"HPSI 1.5"
      nn=Rdim*(Rdim+1)/2
              ZRIJKL1 =ZRIJKL
              ZRIJKL =Zero
      CALL MPI_REDUCE(ZRIJKL1,ZRIJKL,nn,MPI_DOUBLE_COMPLEX,MPI_SUM,0,
     .                                          MPI_COMM_WORLD,ierr)
              ZRIJ1 =ZRIJ
              ZRIJ =Zero
      CALL MPI_REDUCE(ZRIJ1,ZRIJ,Rdim,MPI_DOUBLE_COMPLEX,MPI_SUM,0,
     .                                          MPI_COMM_WORLD,ierr)

c            write(6,*)"HPSI 2",SUM(ZRIJ1)
c      CALL MPI_Barrier(MPI_COMM_WORLD); 
c        VOUT=VOUT+VOUT1
c        VIN=VIN1
c        write(6,*)" Norm 2b PSI",DZNRM2(Nconf,VOUT,1)
c==============================================
c       DEALLOCATE(VOUT1)
c        write(6,*)MYID," HPSI is Done "
c      return
c        write(6,'(a20,15F26.16)')"ZRIJ:",REAL(zrij)
c            write(6,'(a20,15F26.16)')"ZRIJ:",REAL(zrij(1))
c            write(6,'(a20,15F26.16)')"ZRIJKL:",REAL(zrijkl(1))
c          write(6,*)"ALL HPSI",MPI_WTIME(ierr)-start_mpi_all
      end subroutine HPSI
