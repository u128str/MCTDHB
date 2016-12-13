c============== This subroutine computes Right Hand Side part of MCTDH equations 
c============== to provide them to Dieter's Integrator 
c==============
c================ Func(NextTime,Psi,DtPsi,CData,RData,IData,LData)
C *   Func:          Computes the time derivative DtPsi of Psi at time *
C *                  AbsTime. Called as Func(AbsTime,Psi,DtPsi,CData,  *
C *                  RData,IData,LData).
c        subroutine FUNCSTR(NextTime,dv0,dv,SO,XDVR,NPr,Ldata)
       SUBROUTINE FUNCSTR(AbsTime,Psi,OPsi)
       use PASS_ARG
       USE SHARED_DIMS
       USE DVR_ALL
       USE Parallel_Orb
       USE rR_hW
       USE W_INTERPARTICLE
       IMPLICIT NONE
c=====================================================================
       INCLUDE 'mpif.h'
       INTEGER ::  ierr,MYID,numprocs
c=====================================================================
       REAL*8  :: AbsTime
       REAL*8  :: local_time
c=====================================================================
c=================== END MPI ================================
       COMPLEX*16,  DIMENSION(NDX*NDY*NDZ,Morb) ::PSI,OPsi
       COMPLEX*16,ALLOCATABLE :: temp(:,:)
!        COMPLEX*16,ALLOCATABLE :: lhs(:,:) !JULIAN LR AIS 18JUL2014
       COMPLEX*16 :: W,E,mu,T1,E1
c============================== Timings  =============================
       real*4   start,finish, exec_time , cnk_time       
       real*8 :: tm,tmb
c=====================================================================
       INTEGER :: K,S,L,Q,n,nn,i,iii,jj,ii
       INTEGER ::cI,CJ,cK,cL,i2,WTYPE 
c========================= ===========================================
       INTEGER :: FromN,TillN,Nadr
       INTEGER :: Iorb,IPRC,Jorb,DM,icntr,Korb,Lorb
c======== Data with 2-body density elements ==========================
       COMPLEX*16 :: RhoAll, Z,z1,z2,ZDOTU,ZDOTC
       EXTERNAL Get_Op_PSI,ZSCAL,Nadr, VTRAP_EXT_TD,ZCOPY
c=====================================================================
       Wtype=Wxx_Type
       tm=MPI_WTIME(ierr)
       tmb=MPI_WTIME(ierr)
c=====================================================================
       DM=NDX*NDY*NDZ        
       n=NDX*NDY*NDZ*MORB    
c=====================================================================
       allocate(temp(NDX*NDY*NDZ,Morb),stat=ierr)
       if(ierr /= 0)write(*,*)"allocation error in FUNCSTR"
!       allocate(lhs(NDX*NDY*NDZ,Morb),stat=ierr)
!       if(ierr /= 0)write(*,*)"lhs allocation error in FUNCSTR"
c=====================================================================


       If(ABS(Real(Prefac)).ge.epsilon(0d0)) then
         IF(Morb.gt.1) call  schmidtortho(PSI,NDX*NDY*NDZ,Morb,ierr)
         IF(Morb.eq.1) THEN
           Z=ZONER/ DSQRT( ABS(ZDOTC(DM,PSI(:,1),1,PSI(:,1),1)))
           CALL ZSCAL(NDX*NDY*NDZ,Z,PSI(:,1),1)
         ENDIF
       endIf
c=============================== MPI =================================
        call MPI_COMM_RANK(MPI_COMM_WORLD,MYID,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,ierr)
        IPRC=MYID+1
                                JOB_TYPE=1
       IF(ENRG_EVAL.eqv..TRUE.) JOB_TYPE=2
       call MPI_BCAST(JOB_TYPE,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

c       CALL  MPI_ScatterV(PSI,Orb_Block,Orb_Displ,MPI_DOUBLE_COMPLEX,
c     . temp,Orb_Block(IPRC),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(PSI,n,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
          temp=PSI
c===================== time for an EXTERNAL TRAP ============================
       local_time=AbsTime
       call MPI_BCAST(local_time,1,MPI_DOUBLE_PRECISION,0,
     .               MPI_COMM_WORLD,ierr)
c=====================================================================
        tm=MPI_WTIME(ierr)
        FromN=Proc_WSL_Starts(IPRC)
        TillN=Proc_WSL_Finish(IPRC)


       IF(DIM_MCTDHB.eq.1.and.Wxx_TYPE.eq.0.and.NDX.le.1024) THEN
       CALL Get_OPSI_WSL_disbalanced_OMP (PSI,OPSI,FromN,TillN,
     &                                    local_time)
       ELSE
       CALL Get_OPSI_WSL_balanced_OMP (PSI,OPSI,FromN,TillN,local_time)
       ENDIF
c=====================================================================
          IF(TillN.ne.Morb*(Morb+1)/2) THEN ! Case WSL MPI - paralellized
      CALL MPI_REDUCE(MPI_In_Place,OPSI,n,MPI_DOUBLE_COMPLEX,MPI_SUM,0,
     .                                          MPI_COMM_WORLD,ierr)
        IF(ENRG_EVAL.eqv..TRUE.) THEN
      CALL MPI_REDUCE(MPI_IN_PLACE,WIJKL,Rdim1,
     .MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        ENDIF
          ENDIF
        timeCPU(7)=timeCPU(7)+MPi_WTIME(ierr)-tm
        tm=MPI_WTIME(ierr)
c=====================================================================
       tm=MPI_WTIME(ierr)


       CALL OPERS(temp,local_time)

c=====================================================================
!       CALL  MPI_GatherV(temp,Orb_Block(IPRC),MPI_DOUBLE_COMPLEX,temp,
       CALL  MPI_GatherV(MPI_In_Place,Orb_Block(IPRC),
     . MPI_DOUBLE_COMPLEX,temp,
     . Orb_Block,Orb_Displ,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
       timeCPU(1)=timeCPU(1)+MPI_WTIME(ierr)-tm
       tm=MPI_WTIME(ierr)
c============================   Evaluation of the H_ij integrals ======================== 
        IF(ENRG_EVAL.eqv..TRUE.) THEN
        E=Zero
        icntr=1
        DO Iorb=1,Morb
        DO Jorb=Iorb,Morb
        HIJ(icntr)=ZDOTC(DM,PSI(:,Iorb),1,temp(:,Jorb),1) ! Original version
c        HIJ(icntr)=ZDOTC(DM,PSI(:,Iorb),1,PSI1(:,Jorb),1) ! Original version
        AllHIJ(Iorb,Jorb)=HIJ(icntr)
        AllHIJ(Jorb,Iorb)=Conjg(AllHIJ(Iorb,Jorb))
        icntr=icntr+1
        EndDO
        EndDO
        ENDIF
c================ Total Energy evaluation only for Real wij and rij
      IF(ENRG_EVAL.eqv..TRUE.) THEN
         W=Zero
       E=SUM(AllHIJ*AllZRIJ)
       W=SUM(ZRIJKL*WIJKL*(DREAL(RESCALE_2B)))
     . +SUM(Conjg(ZRIJKL)*Conjg(WIJKL)*(DIMAG(RESCALE_2B)))
        Energy=Real(E+W)
        ENDIF


c      CALL ZGEMM('N','T',DM,Morb,Morb,ZONER,OPSI,DM, 
c     .                           InvZRIJ,Morb,Zero,PSI2,DM)

      CALL ZGEMM('N','T',DM,Morb,Morb,ZONER,OPSI,DM, 
     .                           InvZRIJ,Morb,Zoner,temp,DM)
!         julian AIS 18JUL2014
!      CALL ZGEMM('N','T',DM,Morb,Morb,ZONER,temp,DM,
!     .                           AllZRIJ,Morb,0.0d0,lhs,DM)


c             opsi=temp

Cc      OPSI=PSI2+PSI1 ! Here I add  kinetic and potential terms
C!$OMP PARALLEL DO PRIVATE(ii,k) SHARED (OPSI,PSI2,temp,ndx,ndy,ndz,Morb)
Cc                Do k=1,Morb
C                Do ii=1, ndx*ndy*ndz
C                Do k=1,Morb
Cc        OPSI(ii,k)=PSI2(ii,k)+PSI1(ii,k)
C        OPSI(ii,k)=PSI2(ii,k)+temp(ii,k)
C                 Enddo
C                 Enddo
C!$END OMP PARALLEL DO
c============================   P R O J E C T I O N ======================== 
C========================= FIRST PROJECTION ===============================
C      DO Iorb=1,Morb
C         DO Jorb=1,Morb
Cc            Zmu(Jorb,Iorb)=-ZDOTC(DM,PSI(:,Jorb),1,OPSI(:,Iorb),1)
C            Zmu(Jorb,Iorb)=-ZDOTC(DM,PSI(:,Jorb),1,temp(:,Iorb),1)
C         EndDO
C      EndDO
C      CALL ZGEMM('N','N',DM,Morb,Morb,ZONER,PSI,DM, 
C     .                               ZMU,Morb,ZONER,temp,DM)
Cc     .                               ZMU,Morb,ZONER,OPSI,DM)

C========================= SECOND  PROJECTION ===============================
       DO Iorb=1,Morb
          DO Jorb=1,Morb
c              Zmu(Jorb,Iorb)=-ZDOTC(DM,PSI(:,Jorb),1,OPSI(:,Iorb),1)
!                ZmuR(Jorb,Iorb)=ZDOTC(DM,PSI(:,Jorb),1,lhs(:,Iorb),1)
              Zmu(Jorb,Iorb)=-ZDOTC(DM,PSI(:,Jorb),1,temp(:,Iorb),1)
          EndDO
       EndDO
       CALL ZGEMM('N','N',DM,Morb,Morb,ZONER,PSI,DM,
     .                               ZMU,Morb,ZONER,temp,DM)
c     .                               ZMU,Morb,ZONER,OPSI,DM)

c      OPSI=Prefac*OPSI !Propagation
!$OMP PARALLEL DO PRIVATE(ii,k) SHARED (OPSI,Prefac,ndx,ndy,ndz,Morb)
c                Do k=1,Morb
                Do ii=1, ndx*ndy*ndz
                Do k=1,Morb
c        OPSI(ii,k)=OPSI(ii,k)*Prefac
        OPSI(ii,k)=temp(ii,k)*Prefac
                 Enddo
                 Enddo
!$END OMP PARALLEL DO
c================================================================
       ABMiter=ABMiter+1
c       deallocate(PSI1,PSI2,temp,tempX,tempY,tempZ,tempIN,tempOUT)
       deallocate(temp)
          timeCPU(2)=timeCPU(2)+MPI_WTIME(ierr)-tm
          timeCPU(3)=timeCPU(3)+MPI_WTIME(ierr)-tmb
      end subroutine  FUNCSTR

