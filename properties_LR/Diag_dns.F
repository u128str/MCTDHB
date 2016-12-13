! Subroutine to provide Op_x, Op_Y,Op,_Z for module DVR_ALL
        SUBROUTINE Diag_DNS
        use PASS_ARG
        USE SHARED_DIMS
        USE rR_hW
        USE CI_ALL
        implicit NONE
       COMPLEX*16, ALLOCATABLE :: AllZRIJ1(:,:),AllZRIJ2(:,:)
       COMPLEX*16, ALLOCATABLE :: AllZRIJ3(:,:),WORK(:)

       REAL*8, ALLOCATABLE :: RWORK(:)

       INTEGER,  ALLOCATABLE :: ipiv(:)

c================================== MPI ======================================
       INTEGER ::  ierr,MYID,numprocs
c=============================================================================
          INTEGER :: I,J,K,icntr,Nadr
          INTEGER :: Iorb,Jorb,iflag,info
          INTEGER :: N,n1,n2,n3,LWORK
c          REAL*8, DIMENSION(Morb) :: Nocc
          REAL*8 :: x,EnrgIN,Errtol,time
          COMPLEX*16 :: ZDOTC
          EXTERNAL DASUM,Nadr,ZGETRF,ZGETRI
          LOGICAL :: SIL,SIMPLE_INV
c=============================================================================
   
c        if(ierr /= 0)write(*,*)"allocation error in Get_r_R  "
        allocate(AllZRIJ1(Morb,Morb),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in Get_r_R"
        allocate(AllZRIJ2(Morb,Morb),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in Get_r_R"
        allocate(AllZRIJ3(Morb,Morb),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in Get_r_R"
        allocate(WORK(Morb*Morb),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in Get_r_R"

        allocate(RWORK(3*Morb-2),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in Get_r_R"
        allocate(ipiv(Morb),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in Get_r_R"

        
c=============================================================================
          AllZRIJ=Zero
          InvZRIJ=Zero
c============= Here I reconstruct AllZRIJ(Morb,Morb)==========================
        icntr=1
        DO Iorb=1,Morb
        DO Jorb=Iorb,Morb

!           print *,iorb,jorb, "ZRhoij",ZRIJ(icntr)
        AllZRIJ(Iorb,Jorb)=ZRIJ(icntr)
        IF(Jorb.ne.Iorb) AllZRIJ(Jorb,Iorb)=Conjg(ZRIJ(icntr)) ! Bug fixed at 27th of November it does not impact in propagation because
        icntr=icntr+1
        EndDO
        EndDO
        DO Iorb=1,Morb
        DO Jorb=1,Morb
!           print *,iorb,jorb, "Rhoij",AllZRIJ(Iorb,Jorb)
        EndDO
        EndDO
c================== Natural Occupations and Natural Vectors ==================
        AllZRIJ1=AllZRIJ
       CALL ZHEEV('V','U',Morb,AllZRIJ1,Morb,NOcc,WORK,Morb*Morb,
     .                                    RWORK,INFO)
       NatVec=AllZRIJ1
!           do n1=1,Morb
!           print *, "Natvec",natvec(:,n1)
!           enddo
         deallocate(AllZRIJ1,AllZRIJ2,AllZRIJ3,WORK,RWORK,ipiv)
         END  SUBROUTINE Diag_DNS
