! Subroutine to provide Op_x, Op_Y,Op,_Z for module DVR_ALL
        SUBROUTINE Get_r_R(VIN,VOUT,time)
!        use PASS_ARG
        USE SHARED_DIMS
        USE rR_hW
        USE CI_ALL
        implicit NONE
      INTERFACE 
      SUBROUTINE HPSI(VIN,VOUT)
       COMPLEX*16 ::  VIN(:),VOUT(:)
      END SUBROUTINE HPSI
      END INTERFACE 
c        COMPLEX*16, DIMENSION(Nconf) :: VIN
c        COMPLEX*16, DIMENSION(Nconf) :: VOUT
        COMPLEX*16 :: VIN(:),VOUT(:)
c        COMPLEX*16, DIMENSION(Nconf) :: VOUT1
        COMPLEX*16, DIMENSION(:), ALLOCATABLE :: VOUT1

c        COMPLEX*16, DIMENSION(Morb,Morb) :: AllZRIJ1,AllZRIJ2,AllZRIJ3
c        COMPLEX*16, DIMENSION(Morb*Morb) :: WORK
       COMPLEX*16, ALLOCATABLE :: AllZRIJ1(:,:),AllZRIJ2(:,:)
       COMPLEX*16, ALLOCATABLE :: AllZRIJ3(:,:),WORK(:)

c        REAL*8, DIMENSION(3*Morb-2) :: RWORK
       REAL*8, ALLOCATABLE :: RWORK(:)

c        INTEGER, DIMENSION(Morb) :: ipiv
       INTEGER,  ALLOCATABLE :: ipiv(:)

c================================== MPI ======================================
       INCLUDE 'mpif.h'
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
        allocate(VOUT1(Nconf),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in Get_r_R"
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
       call MPI_COMM_RANK(MPI_COMM_WORLD,MYID,ierr)
       call MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,ierr)
c       write(6,*)"Get_r_R Process ",MYID," of ", numprocs
c           pause
c=============================================================================
       CI_SCF=0 ! CI_SCF=0 -> CI PART
       call MPI_BCAST(CI_SCF,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c========================= MPI PART ======================================
c       n=Morb*(Morb+1)/2
c       call MPI_BCAST(HIJ,n,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
c       call MPI_BCAST(MaxTrm1b,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c       call MPI_BCAST(TERM_REQ_1B,MaxTrm1b,MPI_DOUBLE_COMPLEX,0,
c     .                                           MPI_COMM_WORLD,ierr)
cc       n=Morb*(Morb+1)/2*(Morb*(Morb+1)/2+1)/2
c       call MPI_BCAST(WIJKL,n,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
c       call MPI_BCAST(MaxTrm2b,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c       call MPI_BCAST(TERM_REQ_2B,MaxTrm2b,MPI_DOUBLE_COMPLEX,0,
c     .                                           MPI_COMM_WORLD,ierr)
c       write(6,*)MYID,"Get_r_R before share h w"
          call Share_H_W
c       write(6,*)MYID,"oK Get_r_R after share h w"
           VOUT1=VIN
c          write(6,*)"Vin Get_r_R",Vin(176751)
c          write(6,*)"Vin",DREAL(ZDOTC(Nconf,VIN,1,VIN,1))
c           write(6,*)"hij",HIJ(1),HIJ(5),HIJ(8),HIJ(10)
c           write(6,'(2F26.16)')SUM(Vout1),SUM(VIN)
c           write(6,*)"Get_r_R 0",size(VOUT1),VOUT1(1),VOUT1(NCONF)
c            pause 
       CALL HPSI(VOUT1,VOUT)
c            write(6,*)"Get_r_R 1",size(VOUT),VOUT(1),VOUT(NCONF)
c            pause 
c          write(6,*)"Vout1 Get_r_R",Vin (176751)
c          write(6,*)"Vout  Get_r_R",Vout(176751)
        EnrgIN=DREAL(ZDOTC(Nconf,VIN,1,VOUT,1))
        Energy=EnrgIN
c          write(6,'(a50,F36.16,a20,F26.16)')
c     ."In Get_r_R Total CI Energy=",EnrgIN," At T=",time
c           write(6,'(a10,4F26.16)')"get_r_r",Real(SUM(Vout1)),SUM(VOUT)
c           write(6,*)"rho1",ZRIJ
c           write(6,*)"rho2",ZRIJKL
c           write(6,*)"VIN",VIN
c           write(6,*)"VOUT",VOUT

c          pause
c=============================================================================
!!!!!!   Here after CALL HPSI(VIN,VOUT) all ZRIJ and ZRIJKL are known  !!!!!!! 
c=============================================================================
          AllZRIJ=Zero
          InvZRIJ=Zero
c============= Here I reconstruct AllZRIJ(Morb,Morb)==========================
        icntr=1
        DO Iorb=1,Morb
        DO Jorb=Iorb,Morb
cOrg        AllZRIJ(Iorb,Jorb)=ZRIJ(icntr)
cOrg        AllZRIJ(Jorb,Iorb)=Conjg(ZRIJ(icntr))

        AllZRIJ(Iorb,Jorb)=ZRIJ(icntr)
        IF(Jorb.ne.Iorb) AllZRIJ(Jorb,Iorb)=Conjg(ZRIJ(icntr)) ! Bug fixed at 27th of November it does not impact in propagation because

c        AllZRIJ(Jorb,Iorb)=ZRIJ(icntr)
c        IF(Jorb.ne.Iorb) AllZRIJ(Iorb,Jorb)=Conjg(ZRIJ(icntr)) ! Bug fixed at 27th of November it does not impact in propagation because
                                                               ! ZHEMM used did not adressed the complex part of diagonal elements
        icntr=icntr+1
        EndDO
        EndDO
c         AllZRIJ3=Transpose(AllZRIJ)
c       AllZRIJ=0.05
c       AllZRIJ(1,1)=20
c       AllZRIJ(2,2)=2
c       AllZRIJ(3,3)=0.5
c       AllZRIJ(4,4)=0.1
c       AllZRIJ(5,5)=0.09
c       AllZRIJ(6,6)=0.07
c       AllZRIJ(7,7)=Npar-AllZRIJ(1,1)-AllZRIJ(2,2)-AllZRIJ(3,3)
c     .                  -AllZRIJ(4,4)-AllZRIJ(5,5)-AllZRIJ(6,6)
c============= Here I (re)construct InvZRIJ(Morb,Morb)==========================
         AllZRIJ1=AllZRIJ
         InvZRij=Zero
         Do I=1,Morb
         InvZRij(I,I)=(1.0d0,0d0)
         ENDDO
c       write(6,*)"Before INverse density",AllZRIJ1
         call ZHESV('U',Morb,Morb,AllZRIJ1,Morb,ipiv,InvZRij,
     .                          Morb,WORK,Morb*Morb,iflag)
c       write(6,*)"After INverse density",iflag
c       write(6,*)"                     density",AllZRij
            
c       write(6,*)"             INverse density",InvZRij

c       write(6,*)"           Unit? ",MATMUL(InvZRij,AllZRij1)
c       write(6,*)"           Unit? ",MATMUL(AllZRij,InvZRij)
c       write(6,*)"           Unit? ",SUM(MATMUL(AllZRij,InvZRij))
c       write(6,*)"           Unit? ",MATMUL(AllZRij,InvZRij)
c====================== Inverse U ===========================
         InvZRij=AllZRIJ
         SIMPLE_INV=.False.
        CALL ZGETRF( Morb, Morb , InvZRij, Morb, IPIV, INFO )  ! LU Factorization
        iF(INFO.ne.0) write(6,*)"Problem in Invertion of U 1",INFO
c        iF(INFO.ne.0) pause
        LWORK=Morb*Morb
        CALL ZGETRI( Morb, InvZRij, Morb, IPIV, WORK, LWORK, INFO )
        iF(INFO.ne.0) write(6,*)"Problem in Invertion of U 2",INFO
        iF(INFO.ne.0) SIMPLE_INV=.True.
        
c       write(6,*)"After INV",iNFO
c================== Natural Occupations and Natural Vectors ==================
        AllZRIJ1=AllZRIJ
       CALL ZHEEV('V','U',Morb,AllZRIJ1,Morb,NOcc,WORK,Morb*Morb,
     .                                    RWORK,INFO)
       NatVec=AllZRIJ1
c       write(6,*)"        Unitary ",SUM(MATMUL(Transpose(Conjg(NatVec)),
c     .MATMUL(AllZRij,NatVec)))
       ErrTOL=1.0d-11
       ErrTOL=Time_Min_Occ_Allowed
c       ErrTOL=epsilon(0.d0)
c       write(6,*)info,"Sorted Natural occupation numbers" ,Nocc
       Do I =1,Morb
       if(Nocc(I).le.ErrTOL) then
       iflag=-I
       goto 11
       Endif
       EndDo
11          continue 
c               IF(SIMPLE_INV) THEN
c                 write(6,*)"SIMPLE_INV",SIMPLE_INV
c                   pause
               IF(iflag.ne.0) Then
       write(6,*)"!!!! PROBLEM !!!!! At least NOCC(",iflag,
     .")=0 so simplified inverse is actived:  "
      write(6,*)"1) All Nocc(i).le.",ErrTOL," are artificially put to 0"
      write(6,*)"2) Rho_INV=0 and Rho_INV(I,I)==1/Nocc(i) are Forced"
      write(6,*)"3) Decrease ErrTOL to 1.0d-12 in Get_r_R.F line 166"
c      ErrTOL=1.0d-10
       InvZRij=Zero
       Do I =1,Morb
c      InvZRij(I,I)=dcmplx(Exp(-Nocc(I)/ErrTOL))
c      IF(Nocc(I).ge.ErrTOL) then
c      InvZRij(I,I)=Zoner/Nocc(I)
         IF(ABS(AllZRIJ(I,I)).ge.ErrTOL) then
       InvZRij(I,I)=Zoner/AllZRIJ(I,I)
         ElSE
       InvZRij(I,I)=Zero
c      InvZRij(I,I)=dcmplx(Exp(-Nocc(I)/ErrTOL))
       write(6,*)"Nocc(",I,")=0d0"
         ENDIF
       EndDO
C       CALL ZGEMM('N','N',Morb,Morb,Morb,Zoner,AllZRIJ1,Morb,
C     .                                 InvZRij,Morb,Zero,AllZRIJ2,Morb)
C       CALL ZGEMM('N','N',Morb,Morb,Morb,Zoner,AllZRIJ2,Morb,
C     .                                AllZRIJ1,Morb,Zero,InvZRij,Morb)

c       write(6,*)"Renormalized INverse density",InvZRij
c       write(6,*)"           Unit? ",MATMUL(InvZRij,AllZRij)
c       write(6,*)"           Unit? ",MATMUL(AllZRij,InvZRij)
c        write(6,*)"How good is inverse SUM D^-1 D=Morb?",
c     . SUM(MATMUL(AllZRij,InvZRij))
                         ENDIF
c                         ENDIF

c        write(6,*)"How good is inverse SUM D^-1 D=Morb?"
         IF(print_level.ge.3) then
        IF(ABS(SUM(MATMUL(AllZRij,InvZRij))-Morb).ge.0.0000001d0)
     .   write(6,*)"How good is inverse SUM D^-1 D=Morb?",
     . SUM(MATMUL(AllZRij,InvZRij))
                ENDIF !PRINT_LEVEL
c     .,InvZRij,AllZRij
c=============================================================================
        SIL=.FALSE.
        call MPI_BCAST(SIL,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)

         deallocate(VOUT1,AllZRIJ1,AllZRIJ2,AllZRIJ3,WORK,RWORK,ipiv)
         END  SUBROUTINE Get_r_R
