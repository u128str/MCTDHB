C=============== INtegral Hij Wijkl Evaluation 
       subroutine get_h_W(PSI,time)
c       use PASS_ARG
       USE SHARED_DIMS
       USE rR_hW
       USE CI_All
       USE Parallel_Orb
       implicit NONE
c=================== MPI 
       INCLUDE 'mpif.h'
       INTEGER ::  ierr,MYID,numprocs,i
c==========================================================
c       COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb)  :: PSI
c       COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb)  :: OPSI
c       COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb)  :: OPSI
c        COMPLEX*16 ::  PSI(:,:)
       COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb)  :: PSI
        COMPLEX*16, ALLOCATABLE ::  OPSI(:,:)
c==========================================================
       integer n
       REAL*8 :: time
       COMPLEX*16 :: Z1,ZDOTC
c==========================================================
         allocate(OPSI(NDX*NDY*NDZ,Morb),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in OPSI get_h_w"
c==========================================================
         call MPI_COMM_RANK(MPI_COMM_WORLD,MYID,ierr)
         call MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,ierr)
c         write(6,*)"Get H-W Process ",MYID," of ", numprocs
c==========================================================
        CI_SCF=1 ! CI_SCF=1 -> SCF ORBITAL  PART
        call MPI_BCAST(CI_SCF,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c===========================================================================0
c       n=Morb*(Morb+1)/2
c       call MPI_BCAST(ZRIJ,n,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
c       call MPI_BCAST(MaxTrm1b,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c       call MPI_BCAST(TERM_REQ_1B,MaxTrm1b,MPI_DOUBLE_COMPLEX,0,
c     .                                           MPI_COMM_WORLD,ierr)
c       n=Morb*(Morb+1)/2*(Morb*(Morb+1)/2+1)/2
c       call MPI_BCAST(ZRIJKL,n,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
c       call MPI_BCAST(MaxTrm2b,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c       call MPI_BCAST(TERM_REQ_2B,MaxTrm2b,MPI_DOUBLE_COMPLEX,0,
c     .                                           MPI_COMM_WORLD,ierr)
c      n=Morb*Morb
cc      call MPI_BCAST(InvZRIJ,n,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
c      call MPI_BCAST(AllZRIJ,n,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
c         print * , "Before Share_r_R"
        call Share_r_R
c         print * , "After Share_r_R"
c==========================================================
        ENRG_EVAL=.TRUE.
c        JOB_TYPE=2
c       call MPI_BCAST(JOB_TYPE,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c       write(6,*)" We need just MaxTrm1b=",MaxTrm1b
c       write(6,*)" We need just MaxTrm2b=",MaxTrm2b
c        pause

c        pause
c          write(6,*)"Get _h_W 0 is ok"
c     .,SUM(CONJG(PSI(:,1))*PSI(:,1))

c         print * , "Before Func norm PSI(1)"
c     .,ZONER/DSQRT(ABS(ZDOTC(NDX*NDY*NDZ,PSI(:,1),1,PSI(:,1),1)))
       call FUNCSTR(time,PSI,OPsi)
c         print * , "After Func norm PSI(1)"
c     .,ZONER/DSQRT(ABS(ZDOTC(NDX*NDY*NDZ,PSI(:,1),1,PSI(:,1),1)))
c          write(6,*)"Get _h_W 1 is ok"
c       write(6,*)"For ring  we need just =",MaxTrm2b+MaxTrm1b
c       write(6,'(a10,100F16.10)')
c     . "0:H =",(Real(HIJ(TERM_REQ_1B(i))),i=1,MaxTrm1b)
c       write(6,'(a10,900F16.10)')
c     . "0:W =",(Real(WIJKL(TERM_REQ_2B(i))),i=1,MaxTrm2b)
c       pause

        JOB_TYPE=999
       call MPI_BCAST(JOB_TYPE,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

c========================= MPI PART ======================================
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
c========================= END MPI PART ======================================
       deallocate(OPSI)
      end subroutine get_h_W
