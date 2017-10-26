MODULE LINEAR_RESPONSE
USE Prop_MB
USE W_INTERPARTICLE
USE CORRELATIONFUNCTIONS
USE PASS_ARG
USE SHARED_DIMS
USE rR_hW

IMPLICIT NONE
         SAVE                                                                           
         INTEGER, public  :: dimL_orb, dimL, ND, numeig
         INTEGER, DIMENSION(100), public  :: c1_i,c1_j
         INTEGER, DIMENSION(10000), public  :: c2_i,c2_j,c2_k,c2_l
         COMPLEX*16, ALLOCATABLE, public :: Rho_ij(:,:),Rho_ijkl(:,:,:,:), H_ij(:,:), W_ijkl(:,:,:,:)
         REAL*8, public :: dx
!         LOGICAL :: CI_PRD=.FALSE.

CONTAINS

!==============create SIL matrix
!=============== Diagonalize CI-Matrix in Krylov subspace
        SUBROUTINE diag_CI(VIN,VOUT,EV,EVALS,Nc,Error_SIL,stateCI)

        USE CI_ALL
        USE Parallel_CI 
        IMPLICIT NONE
        COMPLEX*16 :: VIN(Nc), VOUT(Nc)
!=================== MPI ==================================
        INCLUDE 'mpif.h'
        INTEGER ::  ierr,MYID,numprocs
!==========================================================
        INTEGER :: i,j,k,l,P,Nc
        INTEGER ::  n,FromN,TillN
!====================== For SIL ===========================
        INTEGER ::  maxsil,NL,Iter_Total,maxsil1,minsil
        COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: VSIL(:,:)
        REAL*8   , DIMENSION(125) :: SILALPHA,SILBETA,SILPROJ
        REAL*8 :: time,Error_SIL, E_state,xerr,t,t1,Error_SIL_WORK
        COMPLEX*16 :: Z,Z1
        COMPLEX*16, DIMENSION(125) :: SCVECP,SCVECM
!==============================================

        REAL*4 start,finish,exec_time ,finish1,iter_time       
        REAL*8 time_start,time_finish
        REAL*8 xnrm

        CHARACTER*20 message
        CHARACTER*16 lname 
        CHARACTER*26 lnamenew,lnamerest  
        LOGICAL  CNV,SIL

        REAL*8 ::   FKTR,CNK,DZASUM,DZNRM2,DNRM2,DDOT
        COMPLEX*16 :: ZDOTC

        REAL*8   , DIMENSION(LR_maxsil-1)     :: EVALS
        REAL*8   , DIMENSION(LR_maxsil-1,LR_maxsil-1)  :: EV
        INTEGER :: stateCI

        EXTERNAL FKTR,CNK
        EXTERNAL DZASUM,DZNRM2,DNRM2,ZDOTC,DDOT,ZSCAL

!========================================================
!===========================================================
        CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,ierr)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,ierr)
        CALL MNGR_CI_Part_julian(numprocs) 

!===========================================================
!        WRITE(6,*) maxtrm1b, maxtrm2b, MYID, numprocs
        Error_SIL_WORK=Error_SIL/1.0d0
        xerr=1.0d0
!============================================
        maxsil=LR_maxsil
        minsil=Time_minsil
        MaxSil1=MaxSil
        ALLOCATE(VSIL(Nc,maxsil))
        IF(ALLOCATED(VSIL).eqv..FALSE.) THEN
           WRITE(6,*) " MEM for SIL", & 
          maxsil," is NOT ok!", SIZE(VSIL)
        END IF

!===================== S I L  =============================================
        SIL=.TRUE.
        Iter_Total=0
!==========================================================================
        VSIL=ZERO
        CNV=.FALSE.
!================== On-line SIL integrator/diagonalizer
        VIN=VIN/DZNRM2(Nc,VIN,1) ! norm
!==================== First Iteration out of DO-LOOP
        time_start=MPI_WTIME(ierr)
        VSIL(:,1)=VIN
        CALL HPSI_LR(VIN,VOUT,Nc)
        VSIL(:,2)=VOUT

!========================== GP case 
        IF(Morb.eq.1) THEN
           SIL=.FALSE.
           CALL MPI_BCAST(SIL,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
           VIN(1) =Zoner
           VOUT(1)=Zoner
           GOTO 102 
        ENDIF

!========================== initialize SIL 
        SILALPHA(1)=DREAL(ZDOTC(Nc,VSIL(:,1),1,VSIL(:,2),1))
        VSIL(:,2)=VSIL(:,2)-VSIL(:,1)*SILALPHA(1)
        SILBETA(1)=DZNRM2(Nc,VSIL(:,2),1)
        VSIL(:,2)=VSIL(:,2)/SILBETA(1)
        E_state=SILALPHA(1)
        Iter_Total=Iter_Total+1

        xerr=SILBETA(1)*t
!        WRITE(6,*) 'SILBETA', SILBETA, 'SILALPHA', SILALPHA

        Do J=2,maxsil-1
!==========================  VSIL(:,I-1)=|PSI>
            VIN=VSIL(:,J)
            VOUT=ZERO
            CALL HPSI_LR(VIN,VOUT,Nc)
            VSIL(:,J+1)=VOUT
!==========================  VSIL(:,I)=H|PSI> r=AV_j
            VSIL(:,J+1)=VSIL(:,J+1)-VSIL(:,J-1)*SILBETA(J-1) 
!==========================  r=r-V_{J-1} Beta_{j-1}
            SILALPHA(J)=DREAL(ZDOTC(Nc,VSIL(:,J),1,VSIL(:,J+1),1))
        SILALPHA(J)=0
        DO i=1,Nc
           SILALPHA(J)=SILALPHA(J)+DREAL(CONJG(VSIL(i,J))*VSIL(i,J+1))
        END DO
            
            VSIL(:,J+1)=VSIL(:,J+1)-VSIL(:,J)*SILALPHA(J) 
!==========================  r=r-V_{J} Alpha_{j}
            SILBETA(J)=DZNRM2(Nc,VSIL(:,J+1),1)
            VSIL(:,J+1)=VSIL(:,J+1)/SILBETA(J)

            NL=J
            SILPROJ=SILALPHA

        ENDDO


        SIL=.FALSE.
        CALL MPI_BCAST(SIL,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
!================= Solve Eigenproblem for Krylov subspace DIM==maxsil
        CALL &
           Get_Eigen_Krylov(NL,SILPROJ,SILBETA,E_state,Error_SIL_WORK, &
           CNV,maxsil,t,EV,EVALS)
! LOWEST EVA!!!
        VOUT=ZERO
        DO J=1,NL
!           DO K=2,2
              VOUT=VOUT+EV(J,stateCI)*VSIL(:,J)
!           ENDDO
        ENDDO

102     CONTINUE
        xnrm=DZNRM2(Nc,VOUT,1)
        VOUT=VOUT/DZNRM2(Nc,VOUT,1)
        WRITE(6,*)"Sil Relax CI vectors are renormalized",xnrm

!==============================================
        DEALLOCATE(VSIL)
        END SUBROUTINE diag_CI  

!========================================================================================================
!================= Diagonalization of SIL subspace Tridiagonal matrix is 
!================= assumed to be constructed out off DIAG0,DIAG1
       SUBROUTINE Get_Eigen_Krylov(NL,D0,D1,E_state,Error_SIL,CNV,MaxSil &
    &             ,t,EV,EVALS)
       IMPLICIT NONE
       REAL*8 :: Error_SIL, E_state,t,ABSTOL,VL,VU
       REAL*8   , DIMENSION(125) :: D0,D1 
       LOGICAL CNV
       INTEGER NL,INFO,MaxSIl,I,J,K,IL,IU,M
       REAL*8   , DIMENSION(NL)     :: EVALS
       REAL*8   , DIMENSION(NL-1)   :: OffDg2
       REAL*8   , DIMENSION(5*NL) :: WORK     
       INTEGER   , DIMENSION(5*NL) :: IWORK     
       INTEGER   , DIMENSION(NL) :: IFAIL     
       REAL*8   , DIMENSION(NL,NL)  :: EV

!========================== GP case 
       IF(Morb.eq.1) THEN
          RETURN 
       ENDIF

!========================== Initialization
       EV=0.0d0
       DO i=1,NL
          EVALS(I)=D0(I)
          EV(I,I)=1.0d0
       ENDDO
       DO i=1,NL-1
          OffDg2(I)=D1(I)
       ENDDO
       ABSTOL=1.0d-14

!========================== Eigensystem solver
       CALL DSTEVX('V','A',NL, D0,OffDg2,VL,VU, IL, IU, ABSTOL, &
    &                        M, EVALS, EV, NL, WORK, IWORK, IFAIL, INFO )

!========================== Tests
       IF(M.ne.NL) WRITE(6,*) "Problems in SIL diag ierr", INFO
       IF(INFO.ne.0) WRITE(6,*) "Problems in SIL diag ierr", INFO
       IF(DABS(E_state-EVALS(STATE)).le.Error_SIL) CNV=.TRUE. 
            
       END SUBROUTINE Get_Eigen_Krylov

!=========================================
!=========================================

       subroutine HPSI_LR(VIN,VOUT,Nc)
       use PASS_ARG
       USE CI_SUBR
       USE CI_ALL
        USE CI_prod 
       implicit NONE
      INTERFACE 
      SUBROUTINE GetCIJKL_All_body_Par(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
      END SUBROUTINE GetCIJKL_All_body_Par
      END INTERFACE 
!=================== MPI 
       INCLUDE 'mpif.h'
       INTEGER ::  ierr,MYID,numprocs,Nc
!==========================================================
       integer i,j,k,l,P,nn

       integer  maxsil,NL,Iter_Total,Nterms
       INTEGER ::  n,FromN,TillN
!=================== F90 F95

        COMPLEX*16, INTENT(IN) :: VIN(Nc)
        COMPLEX*16, INTENT(OUT) :: VOUT(Nc)
!====================== For SIL
        REAL*8 :: time, Error_SIL, E_state

        COMPLEX*16 :: zrho=ZERO,Z,Z1
        COMPLEX*16 :: Escale
!==============================================
      real*4 start,finish,exec_time ,finish1, iter_time       
      REAL*8 :: start_mpi,start_mpi_all
      REAL*8 :: mpi_time_bc,mpi_time_cp,mpi_time_rd

      LOGICAL  CNV,SIL,DNS

!===========================================================
            start_mpi_all=MPI_WTIME(ierr)
         call MPI_COMM_RANK(MPI_COMM_WORLD,MYID,ierr)
         call MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,ierr)
!===========================================================

             start_mpi=MPI_WTIME(ierr)
!=========================================================
               SIL=.TRUE.
      call MPI_BCAST(SIL,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
!=========================================================

              VOUT   =ZERO
              ZRIJKL =Zero
              ZRIJKL1=Zero
              ZRIJ   =Zero 
              ZRIJ1  =Zero
!==================================================================
      call MPI_BCAST(VIN,Nconf,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)

            start_mpi=MPI_WTIME(ierr)
      CALL GetCIJKL_All_body_Par(MYID,VIN)

            start_mpi=MPI_WTIME(ierr)
      CALL MPI_REDUCE(VIN,VOUT,Nconf,MPI_DOUBLE_COMPLEX,MPI_SUM,0, &
    &                                           MPI_COMM_WORLD,ierr)
      nn=Rdim*(Rdim+1)/2
              ZRIJKL1 =ZRIJKL
              ZRIJKL =Zero
      CALL MPI_REDUCE(ZRIJKL1,ZRIJKL,nn,MPI_DOUBLE_COMPLEX,MPI_SUM,0, &
    &                                           MPI_COMM_WORLD,ierr)
              ZRIJ1 =ZRIJ
              ZRIJ =Zero
      CALL MPI_REDUCE(ZRIJ1,ZRIJ,Rdim,MPI_DOUBLE_COMPLEX,MPI_SUM,0, &
    &                                           MPI_COMM_WORLD,ierr)

      end subroutine HPSI_LR

!=========================================
!=========================================

      subroutine READ_matrix_elements()
       
        USE CI_ALL
        USE Parallel_CI 
        USE rR_hW
        USE CI_prod 
        USE CI_ALL 

        IMPLICIT NONE

        INTEGER :: i,j
        LOGICAL :: exist_Rij, exist_muij
        REAL*8, DIMENSION(MaxTrm1b) :: HIJr,HIJi,ZRIJr,ZRIJi
        REAL*8, DIMENSION(MaxTrm2b) :: WIJKLr,WIJKLi,ZRIJKLr,ZRIJKLi
        REAL*8, DIMENSION(Morb,Morb) :: AllZRIJr,AllZRIJi,InvZRIJr,InvZRIJi
        REAL*8, DIMENSION(Morb*Morb) :: ZmuRr, ZmuRi
!=================== MPI ==================================
        INCLUDE 'mpif.h'
        INTEGER ::  ierr,MYID,numprocs

        CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,ierr)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,ierr)
!        CALL MNGR_CI_Part_julian(numprocs) 
        CALL MNGR_CI_Part(numprocs) 
          write(6,*) 'axel', size(ind_ci_1b), size(ind_ci_2b)
         ind_ci_1b=-1 
         ind_ci_2b=-1 
!       READ Rij.out and put together complex vectors

        INQUIRE(file='Rij.out',exist=exist_Rij) 

        IF(exist_Rij.eqv..TRUE.) THEN

          Open (unit=28, File='Rij.out', status='old', action='read', &
    &              form='formatted')

          DO i=1,MaxTrm1b
                 read(28,'(I4,I4,F26.16,F26.16,F26.16,F26.16,F26.16,F26.16)',err=11) c1_i(i),c1_j(i),HIJr(i),HIJi(i), &
    &                  ZRIJr(i),ZRIJi(i)
!                 WRITE(6,*) 'Htest:', HIJr(i),HIJi(i),ZRIJr(i),ZRIJi(i)
                 HIJ(i)=DCMPLX(HIJr(i),HIJi(i))
!                 ZRIJ(i)=DCMPLX(ZRIJr(i),ZRIJi(i))
          END DO

          DO i=1,Morb     
             DO j=1,Morb  
                read(28,*) &!'(F26.16,F26.16,F26.16,F26.16)')  &
       &        AllZRIJr(i,j),AllZRIJi(i,j), InvZRIJr(i,j),InvZRIJi(i,j)
                 
                AllZRIJ(i,j)=DCMPLX(AllZRIJr(i,j),AllZRIJi(i,j))
                InvZRIJ(i,j)=DCMPLX(InvZRIJr(i,j),InvZRIJi(i,j))
             END DO     
          END DO         
     
          DO i=1,MaxTrm2b
                 read(28,'(I4,I4,I4,I4,F26.16,F26.16,F26.16,F26.16)',err=11) c2_i(i),c2_j(i),c2_k(i),c2_l(i),WIJKLr(i),WIJKLi(i), &
    &                  ZRIJKLr(i),ZRIJKLi(i)
!                 WRITE(6,*) 'Wtest:', WIJKLr(i),WIJKLi(i),ZRIJKLr(i),ZRIJKLi(i)
                 WIJKL(i)=DCMPLX(WIJKLr(i),WIJKLi(i))
!                 ZRIJKL(i)=DCMPLX(ZRIJKLr(i),ZRIJKLi(i))
          END DO

          GOTO 12          
 11       WRITE(6,*) 'error in reading rij'
 12       CONTINUE
          CLOSE(28)

        ELSE

           write(6,*)  "File muij is corrupted or empty"

        ENDIF

!       READ muij.out and put together complex vectors

        INQUIRE(file='muij.out',exist=exist_muij) 

        IF(exist_muij.eqv..TRUE.) THEN

          Open (unit=117, File='muij.out', status='old', action='read', &
    &              form='formatted')

          DO i=1,Morb*Morb
                 read(117,'(I4,I4,F16.8,F16.6)',err=13) c1_i(i),c1_j(i),ZmuRr(i),ZmuRi(i)
                 ZmuR(c1_i(i),c1_j(i))=DCMPLX(ZmuRr(i),ZmuRi(i))
          END DO

        ELSE

           write(6,*)  "File muij is corrupted or empty"

        ENDIF

          GOTO 14          
 13       WRITE(6,*) 'error in reading muij'
 14       CONTINUE
          CLOSE(117)

      end subroutine READ_matrix_elements


!==============================================================================
      SUBROUTINE  MNGR_CI_Part_julian(NprocORG)
          USE Parallel_CI
          USE CI_prod
          USE CI_ALL
          IMPLICIT NONE
          INCLUDE 'mpif.h'
!==============================================================================
          INTEGER :: MYID
          INTEGER :: Nterms,FromN,TillN,Master_From, Master_Till,MaxTrm
          INTEGER :: I,J,K,icntr,Nproc,M,NprocORG,ia,ib
          INTEGER ::  ierr,IPRC,TRM
       WRITE(6,*)"CI Parallel Menager CI Part ",NprocORG," processors"
       WRITE(6,*)"CI MaxTrm1b",MaxTrm1b," MaxTrm2b",MaxTrm2b

        CI_Proc_From=0
        CI_Proc_Till=0


        IF(NprocORG.eq.1) THEN
        Nproc=NprocORG
        CI_Proc_From(1)=1
        CI_Proc_Till(1)=MaxTrm2b
        write(6,*)"CI Master does all 2b terms from",CI_Proc_From(1), &
     & " till",CI_Proc_Till(1),"!"
        ENDIF

        Nproc=NprocORG

        MaxTrm=MaxTrm2b+MaxTrm1b
        Master_From=1
        Master_Till=1
         
       WRITE(6,*)"CI MaxTrm",MaxTrm

       IF(NPROC.GT.MaxTrm) stop "NPROC>MaxTrm: reduce nproc!!!"

       DO MYID=0,Nproc-1
       Nterms=MaxTrm/(NPROC)
       FromN=MYID*Nterms+1
       CI_Proc_From(MYID+1)=FromN 
       IF(MYID+1.gt.NPROC-(MaxTrm-Nterms*NPROC)) THEN
       FromN=FromN+(MYID-(NPROC-(MaxTrm-Nterms*NPROC)))
       CI_Proc_From(MYID+1)=FromN
       Nterms=MaxTrm/(NPROC)+1
       ENDIF
       TillN=   FromN+Nterms-1
       CI_Proc_Till(MYID+1)=TillN
       IF(TillN.ge.MaxTrm) THEN
       TillN=   MaxTrm
       CI_Proc_Till(MYID+1)=TillN
       ENDIF
      write(6,*)"CI oRG Proc ",MYID+1," works from",FromN," till",TillN,"!"
       ENDDO


!       EXCEPTIONAL_ID=-1 ! negative means all jobs are nicely redistributed between PROCS
        MYID_TRM=0
       DO I=1,Nproc
      IF(CI_Proc_From(I).le.MaxTrm1b.and.CI_Proc_Till(I).gt.MaxTrm1b) &
   &   EXCEPTIONAL_ID=I-1
      IF(CI_Proc_From(I).le.MaxTrm1b.and.CI_Proc_Till(I).gt.MaxTrm1b) &
   &   MYID_TRM(I)=1

      IF(CI_Proc_Till(I).le.MaxTrm1b) THEN
      CI_Proc_From(I)=CI_Proc_From(I)
      CI_Proc_Till(I)=CI_Proc_Till(I)
      MYID_TRM(I)=MYID_TRM(I)+1
      ELSE
      MYID_TRM(I)=MYID_TRM(I)+2
      CI_Proc_Till(I)=CI_Proc_Till(I)-MaxTrm1b
      ENDIF
      IF(CI_Proc_From(I).ge.MaxTrm1b+1) THEN 
      CI_Proc_From(I)=CI_Proc_From(I)-MaxTrm1b
      ENDIF
       ENDDO
! Proc            1  works from           1  till          42 ! tRM:           3
! Proc            2  works from          43  till         105 ! tRM:           2
! Proc            3  works from         106  till         168 ! tRM:           2
! Proc            4  works from         169  till         231 ! tRM:           2

      write(6,*)"CI Exceptional Proc:",EXCEPTIONAL_ID+1,"works 1b [", &
   &  CI_Proc_From(EXCEPTIONAL_ID+1),":",MaxTrm1b,"] 2b [1", &
   &  ":",CI_Proc_Till(EXCEPTIONAL_ID+1),"]" 
       DO I=1,NprocORG
       write(6,*)"CI Proc",I," works from",CI_Proc_From(I), &
   &  " till",CI_Proc_Till(I),"!"," tRM:",MYID_TRM(I)
       ENDDO

!==================== Allocation of memory for Prefactors_1b,Prefactors_2b,Ind_CI_1b,Ind_CI_2b arrays =====================
!==================== These arrays, different on each node will be used after CI_Production EQ TRUE =======================
            IF(Morb.eq.1) CI_PRD=.FALSE. !For GP there is no need for parallelization

      write(6,*)"CI_Prd ",CI_Prd
           IF(CI_Prd.eqv..TRUE.) THEN
           CI_Production=.FALSE.
           CI_Production_1b=.FALSE.
           CI_Production_2b=.FALSE.
           IF(CI_Production.eqv..FALSE.) THEN
        call MPI_COMM_RANK(MPI_COMM_WORLD,MYID,ierr) ! to get MYID
        IPRC=MYID+1
        FromN=CI_Proc_From(IPRC)
        TillN=CI_Proc_Till(IPRC)

          nXdim=Nconf ! Dimension of the CI vector in use
!         nYdim= ! number of the ci vectors needed work on each node - depends on MYID


           TRM=MYID_TRM(MYID+1)
           TRM_choice: SELECT CASE (TRM)
                    CASE (1) ! Worker on 1b
            nYdim=TillN-FromN+1
        ALLOCATE( Prefactors_1b(nXdim,nYdim) )
        IF(ALLOCATED(Prefactors_1b).eqv..FALSE.) &
     & write(6,*)" MEM for Prefactors_1b is NOT ok!",SIZE(Prefactors_1b)
        ALLOCATE( Ind_CI_1b(nXdim,nYdim) )
        IF(ALLOCATED(Ind_CI_1b).eqv..FALSE.) &
     & write(6,*)" MEM for Ind_CI_1b is NOT ok!",SIZE(Ind_CI_1b)
        IF(ALLOCATED(Prefactors_1b).eqv..FALSE.) stop "Prefactors_1b"
        IF(ALLOCATED(Ind_CI_1b).eqv..FALSE.)     stop "Ind_CI_1b"
            Ind_CI_1b=1
           Prefactors_1b=0d0
                    CASE (2) ! Worker on 2b
            nYdim=TillN-FromN+1
        ALLOCATE( Prefactors_2b(nXdim,nYdim) )
        IF(ALLOCATED(Prefactors_2b).eqv..FALSE. ) &
     & write(6,*)" MEM for Prefactors_2b is NOT ok!",SIZE(Prefactors_2b)
        ALLOCATE( Ind_CI_2b(nXdim,nYdim) )
        IF(ALLOCATED(Ind_CI_2b).eqv..FALSE.) &
     & write(6,*)" MEM for Ind_CI_2b is NOT ok!",SIZE(Ind_CI_2b)
        IF(ALLOCATED(Prefactors_2b).eqv..FALSE.) stop "Prefactors_2b"
        IF(ALLOCATED(Ind_CI_2b).eqv..FALSE.)     stop "Ind_CI_2b"
            Ind_CI_2b=1
           Prefactors_2b=0d0
                    CASE (3) ! EXCEPTIONAL_ID : both 1b and 2b terms on one node
        FromN=CI_Proc_From(EXCEPTIONAL_ID+1)
        TillN=MaxTrm1b
            nYdim=TillN-FromN+1
            write(6,*)"MYID:",MYID,"Ind_CI_1b",nYdim
        ALLOCATE( Prefactors_1b(nXdim,nYdim) )
        IF(ALLOCATED(Prefactors_1b).eqv..FALSE.) &
     & write(6,*)" MEM for Prefactors_1b is NOT ok!",SIZE(Prefactors_1b)
        ALLOCATE( Ind_CI_1b(nXdim,nYdim) )
        IF(ALLOCATED(Ind_CI_1b).eqv..FALSE.) &
     & write(6,*)" MEM for Ind_CI_1b is NOT ok!",SIZE(Ind_CI_1b)
        IF(ALLOCATED(Prefactors_1b).eqv..FALSE.) stop "Prefactors_1b"
        IF(ALLOCATED(Ind_CI_1b).eqv..FALSE.)     stop "Ind_CI_1b"
            Ind_CI_1b=1
           Prefactors_1b=0d0
!       write(6,*)"MYID:",MYID,"Ind_CI_1b",nYdim,"fr:",FromN," tl:",TillN
!            pause
        FromN=1
        TillN=CI_Proc_Till(EXCEPTIONAL_ID+1)
            nYdim=TillN-FromN+1
        ALLOCATE( Prefactors_2b(nXdim,nYdim) )
        IF(ALLOCATED(Prefactors_2b).eqv..FALSE.) &
     & write(6,*)" MEM for Prefactors_2b is NOT ok!",SIZE(Prefactors_2b)
        ALLOCATE( Ind_CI_2b(nXdim,nYdim) )
        IF(ALLOCATED(Ind_CI_2b).eqv..FALSE.) &
     & write(6,*)" MEM for Ind_CI_2b is NOT ok!",SIZE(Ind_CI_2b)
        IF(ALLOCATED(Prefactors_2b).eqv..FALSE.) stop "Prefactors_2b"
        IF(ALLOCATED(Ind_CI_2b).eqv..FALSE.)     stop "Ind_CI_2b"
            Ind_CI_2b=1
           Prefactors_2b=0d0
                    CASE DEFAULT
       write(6,*)" Something wrong in Alloc of MNGR_CI_Part!!!!!!!!!"
                    END SELECT TRM_choice
           ENDIF ! end of CI_Production.eqv.FALSE
           ELSE
           CI_Production_1b=.False.
           CI_Production_2b=.False.
        ALLOCATE( Prefactors_1b(1,1) )
        ALLOCATE( Prefactors_2b(1,1) )
        ALLOCATE( Ind_CI_1b(1,1) )
        ALLOCATE( Ind_CI_2b(1,1) )
        IF(ALLOCATED(Prefactors_1b).eqv..FALSE.) stop "Prefactors_1b"
        IF(ALLOCATED(Ind_CI_1b).eqv..FALSE.)     stop "Ind_CI_1b"
        IF(ALLOCATED(Prefactors_2b).eqv..FALSE.) stop "Prefactors_2b"
        IF(ALLOCATED(Ind_CI_2b).eqv..FALSE.)     stop "Ind_CI_2b"
           ENDIF ! end of CI_Prd
       END  SUBROUTINE MNGR_CI_Part_julian

!===CONSTRUCT CI MATRIX==============================================================

       SUBROUTINE construct_CI_matrix(CImat,Nc)

          USE CI_ALL
          USE Parallel_CI 
          IMPLICIT NONE
          COMPLEX*16 :: Vhelp1(Nc), Vhelp2(Nc)
          COMPLEX*16 :: CImat(Nc,Nc)
!=================== MPI ==================================
          INCLUDE 'mpif.h'
          INTEGER ::  ierr,MYID,numprocs
!==========================================================
          INTEGER :: i,j,k,l,P,Nc

!========================================================
!===========================================================
!          CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,ierr)
!          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,ierr)
!          CALL MNGR_CI_Part_julian(numprocs) 

!===========================================================

          DO i=1,Nc

             Vhelp1=0.0d0
             Vhelp1(i)=1.0d0
             
             CALL HPSI_LR(Vhelp1,Vhelp2,Nc)

             CImat(:,i)=Vhelp2
             CImat(i,i)=CImat(i,i)-Energy
            
          END DO
!STR       
!            CImat=DREAL(CImat)

       END SUBROUTINE construct_CI_matrix

!===DIAGONALIZE EXACT CI MATRIX======================================================

       SUBROUTINE diag_CI_ex(EV_ex,EVALS_ex,Nc,pathLR)

          USE CI_ALL
          USE Parallel_CI 
          IMPLICIT NONE
          COMPLEX*16 :: CImat(Nc,Nc)
!=================== MPI ==================================
          INCLUDE 'mpif.h'
          INTEGER ::  ierr,MYID,numprocs
!==========================================================
          INTEGER :: i,j,k,l,P,Nc
          REAL*8   , DIMENSION(Nc)     :: EVALS_ex
          REAL*8   , DIMENSION(Nc,Nc)  :: EV_ex

          INTEGER :: info
          INTEGER, PARAMETER :: lwork=20000, lrwork=50000, liwork=10000
          COMPLEX*16, DIMENSION(lwork) :: work
          REAL*8, DIMENSION(lrwork) :: rwork
          INTEGER, DIMENSION(liwork) :: iwork
          character*10 pathLR
!========================================================
!===========================================================
!          CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,ierr)
!          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,ierr)
!          CALL MNGR_CI_Part_julian(numprocs) 
!===========================================================

          CALL construct_CI_matrix(CImat,Nconf)

!            WRITE(6,*) Dreal(CImat)
          CALL ZHEEVD('N', 'U', Nc, CImat, Nc, EVALS_ex, work, lwork, rwork, lrwork, iwork,liwork, info)

          EV_ex=CImat

          if(info.ne.0) WRITE(6,*) 'info ZHEEVD:', info, 'work(1): ', & 
     &              work(1), 'lwork: ', lwork, 'rwork(1): ', rwork(1), &
     &             'lrwork: ', lrwork, 'iwork(1): ', iwork(1), &
     &             'liwork: ', liwork

            WRITE(6,*)"Original MCTDHB energy is=",Energy
               OPEN(66,FILE=trim(pathLR)//path_sep//'Pure_MCHB_CI_spectrum.out')
               WRITE(66,*)"Original MCTDHB energy is=",Energy
               Do I=1,Nconf
!               write(6,*)"CI(",I,")=",(dreal(CImat(I,k)),k=1,5)
!               write(6,'(i6,1x,F26.16,1X,F26.16)') i,EVALS_ex(i),EVALS_ex(i)+Energy
               write(66,'(i6,1x,F26.16,1X,F26.16)') i,EVALS_ex(i),EVALS_ex(i)+Energy
               ENDDO
!              pause
                      CLOSE(66)
       END SUBROUTINE diag_CI_ex    

!===DIAGONALIZE EXACT, doubled CI MATRIX===============================================

       SUBROUTINE diag_CI_double_ex(EV_ex,EVALS_ex,Nc)

          USE CI_ALL
          USE Parallel_CI 
          IMPLICIT NONE
          COMPLEX*16 :: CImat(Nc,Nc)
          COMPLEX*16 :: CImat_d(2*Nc,2*Nc)
!=================== MPI ==================================
          INCLUDE 'mpif.h'
          INTEGER ::  ierr,MYID,numprocs
!==========================================================
          INTEGER :: i,j,k,l,P,Nc
          REAL*8   , DIMENSION(Nc)     :: EVALS_ex
          REAL*8   , DIMENSION(Nc,Nc)  :: EV_ex

          INTEGER :: info
          INTEGER, PARAMETER :: lwork=20000, lrwork=50000, liwork=10000
          COMPLEX*16, DIMENSION(lwork) :: work
          REAL*8, DIMENSION(lrwork) :: rwork
          INTEGER, DIMENSION(liwork) :: iwork
!========================================================
!===========================================================
!          CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,ierr)
!          CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,ierr)
!          CALL MNGR_CI_Part_julian(numprocs) 
!===========================================================

          CALL construct_CI_matrix(CImat,Nconf)

          CImat_d=0.0d0
          CImat_d(1:Nconf,1:Nconf)=CImat
          CImat_d(Nconf+1:2*Nconf,Nconf+1:2*Nconf)=-CONJG(CImat)

          CALL ZHEEVD('V', 'U', 2*Nc, CImat_d, 2*Nc, EVALS_ex, work, lwork, rwork, lrwork, iwork,liwork, info)

          EV_ex=CImat_d

          WRITE(6,*) 'info ZHEEVD:', info, 'work(1): ', & 
     &              work(1), 'lwork: ', lwork, 'rwork(1): ', rwork(1), &
     &             'lrwork: ', lrwork, 'iwork(1): ', iwork(1), &
     &             'liwork: ', liwork

          WRITE(6,*)

       END SUBROUTINE diag_CI_double_ex    

!==========h2===========================================

       SUBROUTINE h2_help(h2_out)

          USE CI_ALL
          USE Parallel_CI 
          USE W_INTERPARTICLE !ADD BY STR
          USE DVR_ALL
          USE rR_hW

          IMPLICIT NONE         

          INTEGER :: i, j, k, q, s, ll, mm, kk 

          COMPLEX*16, DIMENSION(dimL_orb/2,dimL_orb/2), INTENT(OUT) :: h2_out
         
          COMPLEX*16, DIMENSION(ND,ND) :: h2s,Tkin
 
!     Non-interacting
         h2_out=0.0d0;
           IF(Time_DVRMETHODX==4) THEN ! FFT - explisit construction of the FFT Kinetik energy - very inefficient :-)
!STR==============
!AIS 21Jul2014
          call T_FFT_1D(Tkin)
          ELSE 
            Tkin=Op_X
          ENDIF 
!              Tkin=Op_X !STR October 2012 modifications
              DO j=1,ND
              Tkin(j,j)=Tkin(j,j)+REAL(VTRAP_EXT(j))
              END DO
         DO k=1,Morb
             h2_out((k-1)*ND+1:k*ND,(k-1)*ND+1:k*ND)=Tkin(:,:)
         END DO

       END SUBROUTINE h2_help


!==========orbital matrix===========================================

       SUBROUTINE construct_orb_matrix(PSI,L_orb,h2_out)

          USE CI_ALL
          USE Parallel_CI 
          USE W_INTERPARTICLE !ADD BY STR
          USE DVR_ALL
          USE rR_hW

          IMPLICIT NONE         

          INTEGER :: i, j, k, q, s, ll, mm, kk 

          COMPLEX*16, DIMENSION(ND,Morb), INTENT(IN) :: PSI
          COMPLEX*16, DIMENSION(dimL_orb,dimL_orb), INTENT(OUT) :: L_orb
          COMPLEX*16, DIMENSION(dimL_orb/2,dimL_orb/2), INTENT(OUT) :: h2_out
         
          COMPLEX*16, DIMENSION(dimL_orb/2,dimL_orb/2) :: h2
          COMPLEX*16, DIMENSION(dimL_orb/2,dimL_orb/2) :: A, B, mumat
          COMPLEX*16, DIMENSION(ND,ND) :: h2s,Tkin
          COMPLEX*16, DIMENSION(ND,ND) :: KSL,KSL1 ! ADDED BY STR to KSL - non-local potential
          COMPLEX*16, DIMENSION(dimL_orb,dimL_orb) :: Pb, Pbb, PL, PLL, L_orb_help1, L_orb_help2
          COMPLEX*16, DIMENSION(dimL_orb/2,dimL_orb/2) :: P
          COMPLEX*16, DIMENSION(Morb*Nd,Morb*Nd) :: TST
          COMPLEX*16, DIMENSION(Morb,Morb) :: muTST
          COMPLEX*16, DIMENSION(ND*Morb)  :: PSITST,PSITST1
          INTEGER, DIMENSION(2) :: ind11, ind12, ind21, ind22
          INTEGER, DIMENSION(2) :: ind11a, ind12a, ind21a, ind22a
          COMPLEX*16, DIMENSION(Morb,Morb) :: mu_PSI, RhoInvSQR, rmu
          COMPLEX*16, DIMENSION(ND) :: psihelp
          COMPLEX*16, DIMENSION(NDX*NDY*NDZ)  :: WSL,WSL1 !ADDED BY STR
          COMPLEX*16, DIMENSION(ND,Morb) :: PSI_help
          REAL*8 :: exp_val
 
           xlambda0=xlambda_0
!          OPEN(55,File='MC_r2mu.out')

!     Build matrix

!     Non-interacting

         A=0.0d0; B=0.0d0; L_orb=0.0d0;Pb=0.0d0;
         mumat=0.0d0; h2=0.0d0; h2s=0.0d0; h2_out=0.0d0;
 
                h2=0d0
          IF(Time_DVRMETHODX==4) THEN ! FFT - explicit construction Tkin matrix
          call T_FFT_1D(Tkin)
          ELSE 
          Tkin=Op_X
          ENDIF

         DO j=1,ND
              Tkin(j,j)=Tkin(j,j)+REAL(VTRAP_EXT(j))
         END DO
         DO k=1,Morb
            DO q=1,Morb
             h2((k-1)*ND+1:k*ND,(q-1)*ND+1:q*ND)=Rho_ij(k,q)*Tkin(:,:)
            END DO
         END DO
!              Checked with HO,SIN,EXP-DVR 1D  
!=============================================================================

         CALL Get_Exp_Val_1bOp(PSI,Tkin,Rho_ij,exp_Val)
         write(6,'(a25,2F20.16)')"<h>_old,<T+V>_recomp:",exp_Val,Real(SUM(Rho_ij*H_IJ))

!     get Lagrange multipliers
!         CALL get_mu(mu_PSI)
              ZmuR= MAtmul(Zmu,Transpose(-Rho_ij)) !AIS 18JUL14 Zmu is already available after execution  Func.F
!         stop
!       CALL Get_h_W(PSI,0d0)
!                mu_PSI=ZmuR
!                ZmuR=mu_PSI
!         write(6,*)"MuR",ZmuR
!         write(6,*)"Mu * Rho",MAtmul(Zmu,Transpose(-Rho_ij))
!     nonlinearities

            DO k=1,Morb
               DO q=1,Morb
                  DO i=1,ND
                     mumat((k-1)*ND+i,(q-1)*ND+i)= &
     &                   ZmuR(k,q)
                  END DO
               END DO
            END DO  

            DO k=1,Morb
               DO q=1,Morb
            
                     DO s=1,Morb
                        DO ll=1,Morb
                           CALL  Get_KSL(KSL,psi(:,s),psi(:,ll))  
                           CALL  Get_KSL(KSL1,CONJG(psi(:,s)),psi(:,ll)) 
                           CALL  Get_WSL_omp(WSL,psi(:,s),psi(:,ll))  !STR corrected back
                  DO i=1,ND
                           A((k-1)*ND+i,(q-1)*ND+i)= A((k-1)*ND+i,(q-1)*ND+i) +Rho_ijkl(k,s,ll,q)*WSL(i) 
                  DO j=1,ND
 
                           A((k-1)*ND+i,(q-1)*ND+j)= &
     &                        A((k-1)*ND+i,(q-1)*ND+j) &
     &                        +Rho_ijkl(k,s,ll,q)*KSL(i,j) 


                           B((k-1)*ND+i,(q-1)*ND+j)= &
     &                        B((k-1)*ND+i,(q-1)*ND+j) &
     &                        +Rho_ijkl(k,q,s,ll)*KSL1(i,j)  
                  END DO
                  END DO
                        
                        END DO
                     END DO

               END DO
            END DO  

!     Assemble whole matrix

            DO k=1,Morb
               DO j=1,Morb
                  DO i=1,ND
                     DO q=1,ND

                     ! for the blocks
                     ind11=(/(k-1)*ND+i,(j-1)*ND+q/)
                     ind12=(/(k-1)*ND+i,Morb*ND+(j-1)*ND+q/)
                     ind21=(/Morb*ND+(k-1)*ND+i,(j-1)*ND+q/)
                     ind22=(/Morb*ND+(k-1)*ND+i,Morb*ND+(j-1)*ND+q/)

                     L_orb(ind11(1),ind11(2))=h2(ind11(1),ind11(2)) &
     &                                -mumat(ind11(1),ind11(2)) &
     &                  +1.0d0*A(ind11(1),ind11(2)) ! Prefactor 2 is only for contact because Wsl Ksl
            TST(ind11(1),ind11(2))=h2(ind11(1),ind11(2)) +0.5d0*A(ind11(1),ind11(2)) !STR added to recompute mu_ij and compare with Func.F results
!     &                  +2.0d0*A(ind11(1),ind11(2))
!     &                  +2.0d0*xlambda_0*A(ind11(1),ind11(2))
!STR     &                  +2.0d0*A(ind11(1),ind11(2))
!ORG     &                  +2.0d0*xlambda_0*A(ind11(1),ind11(2))

                     L_orb(ind12(1),ind12(2))= &
     &                          B(ind11(1),ind11(2))
!      &                          xlambda_0*B(ind11(1),ind11(2))
!STR     &                          B(ind11(1),ind11(2))
!ORG      &                          xlambda_0*B(ind11(1),ind11(2))

!ORG                  L_orb(ind21(1),ind21(2))=-xlambda_0* &
                     L_orb(ind21(1),ind21(2))=- &
     &                          CONJG(B(ind11(1),ind11(2)))

                     L_orb(ind22(1),ind22(2))=- &
     &                   CONJG(L_orb(ind11(1),ind11(2))) 

                     END DO
                  END DO
               END DO
            END DO  
!STR Added to ensure that mu_ij are reproduced!!! YES -they are!!
          Do i=1,Morb
                Do j=1,ND
          PSITST((i-1)*ND+j)=PSI(j,i)
                EndDO
          EndDO
           PSITST1=Matmul(TST,PSITST)
           muTST=0d0
          Do i=1,Morb
          Do j=1,Morb
              Do k=1,Nd
              muTST(i,j)=muTST(i,j)+PSITST1((i-1)*Nd+k)*Conjg(PSI(k,j))
              ENDdo
!            write(6,*)"mu_old-mu_recomut",i,j,real(muTST(i,j)-ZmuR(i,j))
          EndDO
          EndDO
            write(6,*)"mu_old - for Delta only ! mu_recomputed==0",real(SUM(muTST-ZmuR))
!                 pause

!!     nonlinearities, test
!
!         DO k=1,Morb
!         DO kk=1,Morb
!            DO q=1,Morb
!               DO i=1,ND
!
!                     mumat((k-1)*ND+i,(q-1)*ND+i)= &
!     &                mumat((k-1)*ND+i,(q-1)*ND+i) +    InvZRIJ(k,kk)*ZmuR(kk,q)
!                  DO s=1,Morb
!
!
!!                     DO ll=1,Morb
! 
!                           A((k-1)*ND+i,(q-1)*ND+i)= &
!     &                        A((k-1)*ND+i,(q-1)*ND+i) &
!     &                        +InvZRIJ(k,kk)*Rho_ijkl(kk,s,q,ll)*CONJG(PSI(i,s))*PSI(i,ll) 
!                       
!                           B((k-1)*ND+i,(q-1)*ND+i)= &
!     &                        B((k-1)*ND+i,(q-1)*ND+i) &
!     &                        +InvZRIJ(k,kk)*Rho_ijkl(kk,s,q,ll)*PSI(i,s)*PSI(i,ll) 
!                  
!                     END DO
!                  END DO
!
!               END DO
!            END DO
!         END DO  
!         END DO  

!     Calculate projected matrix: PL=Pb*L

!     Get projector matrix

!         IF(Time_DVRMETHODX==4) THEN
!            CALL get_proj(PSI*SQRT(dx),Pb,P)
!         ELSE
            CALL get_proj(PSI,Pb,P)
!         END IF

!     Multiply proj. matrix (double)

         CALL ZGEMM('N', 'N', dimL_orb, dimL_orb, dimL_orb, DCMPLX(1.0d0), Pb, &
     &             dimL_orb, Pb, dimL_orb, DCMPLX(0.0d0), Pbb, dimL_orb)

         CALL ZGEMM('N', 'N', dimL_orb, dimL_orb, dimL_orb, DCMPLX(1.0d0), Pbb, &
     &             dimL_orb, L_orb, dimL_orb, DCMPLX(0.0d0), PL, dimL_orb)
!STR          Pbb=Pb
!STR         CALL ZGEMM('N', 'N', dimL_orb, dimL_orb, dimL_orb, DCMPLX(1.0d0), Pb, &
!STR     &             dimL_orb, L_orb, dimL_orb, DCMPLX(0.0d0), PL, dimL_orb)

         CALL ZGEMM('N', 'N', dimL_orb, dimL_orb, dimL_orb, DCMPLX(1.0d0), PL, &
     & dimL_orb, Pbb, dimL_orb, DCMPLX(0.0d0), PLL, dimL_orb)

         L_orb=PLL
!         goto 2
!      Metric
!      get squareroot of inverse of 

         CALL squarerootRInv(RhoInvSQR)
         

!      before
         L_orb_help1=0.0d0; L_orb_help2=0.0d0
         DO k=1,Morb
            DO s=1,Morb
               DO i=1,ND
                  DO q=1,ND
                     DO j=1,Morb

                  ! for the blocks
                  ind11=(/(j-1)*ND+i,(s-1)*ND+q/)
                  ind12=(/(j-1)*ND+i,Morb*ND+(s-1)*ND+q/)
                  ind21=(/Morb*ND+(j-1)*ND+i,(s-1)*ND+q/)
                  ind22=(/Morb*ND+(j-1)*ND+i,Morb*ND+(s-1)*ND+q/)
                  ind11a=(/(k-1)*ND+i,(s-1)*ND+q/)
                  ind12a=(/(k-1)*ND+i,Morb*ND+(s-1)*ND+q/)
                  ind21a=(/Morb*ND+(k-1)*ND+i,(s-1)*ND+q/)
                  ind22a=(/Morb*ND+(k-1)*ND+i,Morb*ND+(s-1)*ND+q/)

                  L_orb_help2(ind11a(1),ind11a(2))=L_orb_help2(ind11a(1),ind11a(2))+RhoInvSQR(k,j)*L_orb(ind11(1),ind11(2))

                  L_orb_help2(ind12a(1),ind12a(2))=L_orb_help2(ind12a(1),ind12a(2))+RhoInvSQR(k,j)*L_orb(ind12(1),ind12(2))

                  L_orb_help2(ind22a(1),ind22a(2))=L_orb_help2(ind22a(1),ind22a(2))+CONJG(RhoInvSQR(k,j))*L_orb(ind22(1),ind22(2))

                  L_orb_help2(ind21a(1),ind21a(2))=L_orb_help2(ind21a(1),ind21a(2))+CONJG(RhoInvSQR(k,j))*L_orb(ind21(1),ind21(2))

                     END DO
                  END DO
               END DO
            END DO
         END DO

!      after
         DO k=1,Morb
            DO s=1,Morb
               DO i=1,ND
                  DO q=1,ND
                     DO j=1,Morb

                  ! for the blocks
                  ind11=(/(k-1)*ND+i,(j-1)*ND+q/)
                  ind12=(/(k-1)*ND+i,Morb*ND+(j-1)*ND+q/)
                  ind21=(/Morb*ND+(k-1)*ND+i,(j-1)*ND+q/)
                  ind22=(/Morb*ND+(k-1)*ND+i,Morb*ND+(j-1)*ND+q/)
                  ind11a=(/(k-1)*ND+i,(s-1)*ND+q/)
                  ind12a=(/(k-1)*ND+i,Morb*ND+(s-1)*ND+q/)
                  ind21a=(/Morb*ND+(k-1)*ND+i,(s-1)*ND+q/)
                  ind22a=(/Morb*ND+(k-1)*ND+i,Morb*ND+(s-1)*ND+q/)


                  L_orb_help1(ind11a(1),ind11a(2))=L_orb_help1(ind11a(1),ind11a(2))+&
                  L_orb_help2(ind11(1),ind11(2))*RhoInvSQR(j,s)

                  L_orb_help1(ind12a(1),ind12a(2))=L_orb_help1(ind12a(1),ind12a(2))+&
                  L_orb_help2(ind12(1),ind12(2))*CONJG(RhoInvSQR(j,s))

                  L_orb_help1(ind22a(1),ind22a(2))=L_orb_help1(ind22a(1),ind22a(2))+&
                  L_orb_help2(ind22(1),ind22(2))*CONJG(RhoInvSQR(j,s))

                  L_orb_help1(ind21a(1),ind21a(2))=L_orb_help1(ind21a(1),ind21a(2))+&
                  L_orb_help2(ind21(1),ind21(2))*RhoInvSQR(j,s)

                     END DO
                  END DO
               END DO
            END DO
         END DO

         L_orb=L_orb_help1
!2        continue

!STR
!STR         DO i=1,dimL_orb
!STR            DO j=1,dimL_orb
!STR               WRITE(56,'(I8,1X,I8,1X,2E16.8)') i, j, &
!STR     &                          L_orb(i,j)
!STR            END DO
!STR         END DO         
!STR         CLOSE(56)


!STR          CLOSE(55)

      END SUBROUTINE construct_orb_matrix

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     calculate projectors
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE get_proj(PSI,Pb,P)

          USE CI_ALL
          USE Parallel_CI 
          USE DVR_ALL

          IMPLICIT NONE         

          INTEGER :: i, j, k, q, s, ll, mm 

          COMPLEX*16, DIMENSION(ND,Morb), INTENT(IN) :: PSI
          COMPLEX*16, DIMENSION(dimL_orb,dimL_orb), INTENT(OUT) :: Pb
          COMPLEX*16, DIMENSION(dimL_orb/2,dimL_orb/2), INTENT(OUT) :: P

          INTEGER, DIMENSION(2) :: ind11, ind12, ind21, ind22
          INTEGER, DIMENSION(2) :: ind11a, ind12a, ind21a, ind22a

          P=0.0d0
     
          DO i=1,ND
             DO j=1,ND
                DO k=1,Morb

                   P(i+(k-1)*ND,j+(k-1)*ND)=-SUM(PSI(i,:) &
     &                 *CONJG(PSI(j,:)))
!ORG     &                 *CONJG(PSI(j,:)))*dx

                END DO

             END DO

            DO k=1,Morb
               P(i+(k-1)*ND,i+(k-1)*ND)=P(i+(k-1)*ND,i+(k-1)*ND) &
     &              +1.0d0
            END DO

         END DO

!     Big projector matrix

         DO k=1,Morb
            DO j=1,Morb
               DO i=1,ND
                  DO q=1,ND

                  ! for the blocks
                  ind11=(/(k-1)*ND+i,(j-1)*ND+q/)
                  ind12=(/(k-1)*ND+i,Morb*ND+(j-1)*ND+q/)
                  ind21=(/Morb*ND+(k-1)*ND+i,(j-1)*ND+q/)
                  ind22=(/Morb*ND+(k-1)*ND+i,Morb*ND+(j-1)*ND+q/)

                  Pb(ind11(1),ind11(2))=P(ind11(1),ind11(2))

                  Pb(ind12(1),ind12(2))=0.0d0

                  Pb(ind21(1),ind21(2))=0.0d0

                  Pb(ind22(1),ind22(2))=CONJG(P(ind11(1),ind11(2)))

                  END DO
               END DO
            END DO
         END DO  

      END SUBROUTINE get_proj

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     calculate chemical potentials
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE get_mu(mu_PSI)

       USE CI_ALL
       USE Parallel_CI 
       USE DVR_ALL


          IMPLICIT NONE         

          INTEGER :: i, j, k, q, s, ll, mm 

          COMPLEX*16, DIMENSION(Morb,Morb), INTENT(OUT) :: mu_PSI
          COMPLEX*16, DIMENSION(Morb,Morb) :: rmu
          COMPLEX*16 :: en_PSI

          OPEN(72,FILE='MC_mu.out')
          OPEN(73,FILE='MC_rmu.out')
          OPEN(74,FILE='MC_en.out')
          OPEN(75,FILE='MC_Hij.out')
          OPEN(76,FILE='MC_Rij.out')
          OPEN(77,FILE='MC_Wijkl.out')
          OPEN(78,FILE='MC_Rijkl.out')

          mu_PSI=0.0d0; rmu=0.0d0; en_PSI=0.0d0

          DO k=1,Morb
             DO j=1,Morb

                DO q=1,Morb

                   mu_PSI(k,j)=mu_PSI(k,j) + Rho_ij(k,q)*H_IJ(j,q) 
                END DO


                DO q=1,Morb
                   DO s=1,Morb
                      DO ll=1,Morb

                         mu_PSI(k,j)=mu_PSI(k,j) + Rho_ijkl(k,s,q,ll)*W_ijkl(j,s,q,ll)       
                      END DO
                   END DO
                END DO

                WRITE(72,'(I4,3X,I4,3X,2E16.8)') k,j,mu_PSI(k,j)

             END DO
          END DO

          DO k=1,Morb
             DO j=1,Morb

                   en_PSI=en_PSI + Rho_ij(k,j)*H_IJ(k,j) 

                         WRITE(75,'(I4,2X,I4,2X,2E16.8)') k,j,H_IJ(k,j)
                         WRITE(76,'(I4,2X,I4,2X,2E16.8)') k,j,Rho_ij(k,j)
                DO q=1,Morb
                   DO s=1,Morb

                         en_PSI=en_PSI + 0.5d0*Rho_ijkl(k,j,q,s)*W_ijkl(k,j,q,s)       
                         WRITE(77,'(I4,2X,I4,2X,I4,2X,I4,2X,2E16.8)') &
      &                   k,j,q,s,W_ijkl(k,j,q,s)
                         WRITE(78,'(I4,2X,I4,2X,I4,2X,I4,2X,2E16.8)') &
      &                   k,j,q,s,Rho_ijkl(k,j,q,s)
                   END DO
                END DO

             END DO
          END DO
          WRITE(74,'(F26.16, F26.16)') Real(en_PSI), Real(SUM(Rho_ij*H_IJ)+0.5*SUM(Rho_ijkl*W_ijkl)) 
!Here there is a problem RESCALE_2B are forgotten!!! 24/10 2012 
!-fixed By removing 0.5 in rewriting the matrix elements W_ijkl

          
         ! Test rhoinv * mu
          DO k=1,Morb
             DO j=1,Morb
                DO q=1,Morb
                   rmu(k,j)=rmu(k,j)+InvZRIJ(k,q)*mu_PSI(q,j)
                END DO
                WRITE(73,'(I4,3X,I4,3X,2E16.8)') k,j,rmu(k,j)
             END DO
          END DO

!          CLOSE(72-77)
          CLOSE(72)
          CLOSE(73)
          CLOSE(74)
          CLOSE(75)
          CLOSE(76)
          CLOSE(77)
          CLOSE(78)

      END SUBROUTINE get_mu

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     calculate chemical potentials test
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE get_mu_test(mu_PSI)

          USE CI_ALL
          USE Parallel_CI 
          USE DVR_ALL

          IMPLICIT NONE         

          INTEGER :: i, j, k, q, s, ll, mm 

          COMPLEX*16, DIMENSION(Morb,Morb), INTENT(OUT) :: mu_PSI
          COMPLEX*16, DIMENSION(Morb,Morb) :: rmu
          COMPLEX*16 :: en_PSI

          OPEN(72,FILE='MC_mu.out')
          OPEN(73,FILE='MC_rmu.out')
          OPEN(74,FILE='MC_en.out')
          OPEN(75,FILE='MC_Hij.out')
          OPEN(76,FILE='MC_Rij.out')
          OPEN(77,FILE='MC_Wijkl.out')
          OPEN(78,FILE='MC_Rijkl.out')

          mu_PSI=0.0d0; rmu=0.0d0; en_PSI=0.0d0

          DO k=1,Morb
             DO j=1,Morb

                DO q=1,Morb

                   mu_PSI(k,j)=mu_PSI(k,j) + Rho_ij(k,q)*H_IJ(j,q) 
                END DO


                DO q=1,Morb
                   DO s=1,Morb
                      DO ll=1,Morb

                         mu_PSI(k,j)=mu_PSI(k,j) + Rho_ijkl(k,s,q,ll)*W_ijkl(j,s,q,ll)       
                      END DO
                   END DO
                END DO

                WRITE(72,'(I4,3X,I4,3X,2E16.8)') k,j,mu_PSI(k,j)

             END DO
          END DO

          DO k=1,Morb
             DO j=1,Morb

                   en_PSI=en_PSI + Rho_ij(k,j)*H_IJ(k,j) 

                         WRITE(75,'(I4,2X,I4,2X,2E16.8)') k,j,H_IJ(k,j)
                         WRITE(76,'(I4,2X,I4,2X,2E16.8)') k,j,Rho_ij(k,j)
                DO q=1,Morb
                   DO s=1,Morb

                         en_PSI=en_PSI + 0.5d0*Rho_ijkl(k,j,q,s)*W_ijkl(k,j,q,s)       
                         WRITE(77,'(I4,2X,I4,2X,I4,2X,I4,2X,2E16.8)') &
      &                   k,j,q,s,W_ijkl(k,j,q,s)
                         WRITE(78,'(I4,2X,I4,2X,I4,2X,I4,2X,2E16.8)') &
      &                   k,j,q,s,Rho_ijkl(k,j,q,s)
                   END DO
                END DO

             END DO
          END DO
                WRITE(74,'(2E16.8)') en_PSI

          

         ! Test rhoinv * mu
          DO k=1,Morb
             DO j=1,Morb
                   rmu(k,j)=mu_PSI(k,j)
                WRITE(73,'(I4,3X,I4,3X,2E16.8)') k,j,rmu(k,j)
             END DO
          END DO

!          CLOSE(72-77)
          CLOSE(72)
          CLOSE(73)
          CLOSE(74)
          CLOSE(75)
          CLOSE(76)
          CLOSE(77)
          CLOSE(78)

      END SUBROUTINE get_mu_test

!=================================================
!======Diagonalize orbital part====================
!=================================================

      SUBROUTINE diag_orb_matrix(PSI,L_orb,wo,u,v)

          USE CI_ALL
          USE Parallel_CI 
          USE DVR_ALL

          IMPLICIT NONE         

          INTEGER :: i, j, k, q, s, ll, mm 

          COMPLEX*16, DIMENSION(ND,Morb), INTENT(IN) :: PSI
          COMPLEX*16, DIMENSION(dimL_orb,dimL_orb), INTENT(IN) :: L_orb
         
          COMPLEX*16, DIMENSION(dimL_orb), INTENT(OUT) :: wo 
          COMPLEX*16, DIMENSION(ND,numeig,Morb), INTENT(OUT) :: u, v

          COMPLEX*16, DIMENSION(dimL_orb,dimL_orb) :: Lh

          REAL*8, DIMENSION(dimL_orb) :: wh
          COMPLEX*16, DIMENSION(dimL_orb) :: w

          COMPLEX*16, DIMENSION(dimL_orb/2,dimL_orb/2) :: A, B, h2, mumat, P
          COMPLEX*16, DIMENSION(dimL_orb,dimL_orb) :: Pb, Pbb, PL, PLL
         INTEGER, DIMENSION(2) :: ind11, ind12, ind21, ind22
          COMPLEX*16, DIMENSION(Morb,Morb) :: mu_PSI

 !        INTEGER, PARAMETER :: lwork=245760, lwork2=32010, lwork4=9600
         INTEGER, PARAMETER :: lwork=2045760, lwork2=302010, lwork4=90600
         INTEGER :: info
         COMPLEX*16 :: z
         COMPLEX*16, DIMENSION(dimL_orb-1) :: tau
         INTEGER, DIMENSION(dimL_orb) :: ind
         COMPLEX*16, DIMENSION(lwork) :: work
         COMPLEX*16, DIMENSION(lwork2) :: work2
         COMPLEX*16, DIMENSION(lwork4) :: work4
         REAL*8, DIMENSION(dimL_orb) :: rwork3
         COMPLEX*16, DIMENSION(dimL_orb*dimL_orb) :: work3
 
         LOGICAL, DIMENSION(dimL_orb) :: select
         INTEGER, DIMENSION(numeig) :: ifaill, ifailr
         COMPLEX*16 :: vl
         COMPLEX*16, DIMENSION(dimL_orb,numeig) :: vr
         INTEGER :: mn

         COMPLEX*16, DIMENSION(numeig,numeig) :: norm_mat
         COMPLEX*16, DIMENSION(numeig,numeig) :: orth_mat
         INTEGER, DIMENSION(numeig) :: indvec, indorder, indorder2
         REAL*8, DIMENSION(numeig) :: indvecr
         INTEGER :: nz
        
         OPEN (20,FILE='MC_w_test.out')
         OPEN (21,FILE='BMF_wtest.out')
         OPEN (22,FILE='BMF_ev_u.out')
         OPEN (23,FILE='BMF_ev_v.out')
         OPEN (28,FILE='BMF_norm.out')

!     Diagonalize matrix (eigenvalues)

         CALL ZGEHRD(dimL_orb, 1, dimL_orb, L_orb, dimL_orb, tau, work, lwork, info)

         Lh=L_orb  ! save projected matrix for later (eigenvector calculation) (upper Hessenberg form)

         WRITE(6,*)
         WRITE(6,*) 'info_lr_diag: ', info, 'work(1): ', work(1), & 
     &             'lwork: ', lwork
         WRITE(6,*)

         CALL ZHSEQR('E', 'N', dimL_orb, 1, dimL_orb, L_orb, dimL_orb, w, z, 1, &
     &                work2, lwork2, info)
         
         WRITE(6,*)
         WRITE(6,*) 'info_lr_diag_2": ', info, 'work(1): ', work2(1), &
     &             'lwork: ', lwork2
         WRITE(6,*)

!     Sort eigenvalues

         DO i=1,dimL_orb
            WRITE(21,'(I8,1X,2E16.8)') i, w(i)
         END DO

         wh=REAL(w)

         DO i=1,dimL_orb
            ind(i)=i
         END DO

         CALL SHELL(dimL_orb,wh,ind) ! "ind" gives location of i-th value
        
         DO i=dimL_orb/2+1,dimL_orb
         WRITE(20,'(I8,1X,I8,1X,2E16.8,1x,E16.8)') i-dimL_orb/2, ind(i), w(ind(i)),abs(w(ind(i)))
            wo(i-dimL_orb/2)=w(ind(i))
         END DO

!     Compute selected eigenvectors

!     Define subset of eigenvectors to be computed

         indvec=ind(dimL_orb/2+1:dimL_orb/2+numeig) ! positive ones
!         indvec=ind(1:mv) ! first ones
        
!     Compute order of subset: 

         DO i=1,numeig
            indorder(i)=i
            indorder2(i)=i
         END DO

         indvecr=REAL(indvec) ! subset of ind 731 732 730
         CALL SHELL(numeig,indvecr,indorder) ! gives the places of each rank 3 1 2  
         indvecr=REAL(indorder)  
         CALL SHELL(numeig,indvecr,indorder2) ! gives the rank of each place 2 3 1        
        
!     Calculate number of zero vectors
 
         nz=0;
         DO i=1,numeig
            IF(ABS(wo(i))<1.0d-1) THEN
               nz=nz+1;
            ELSE
               EXIT
            END IF
               WRITE(6,*) 'nz', nz, ABS(wo(i))
         END DO

!     Eigenvector calculation
  
         vl=0.0d0; vr=0.0d0
         select=.FALSE.
         select(indvec)=.TRUE.
         mn=1
         work3=0.0d0; rwork3=0.0d0

         CALL ZHSEIN('R', 'Q', 'N', select, dimL_orb, Lh, dimL_orb, w, vl, & 
     &                  1, vr, dimL_orb, numeig, mn, work3, rwork3, ifaill, &
     &                  ifailr, info)

         WRITE(6,*)
         WRITE(6,*) 'info_lr_ev: ', info
         WRITE(6,*)

         CALL ZUNMHR('L', 'N', dimL_orb, numeig, 1, dimL_orb, Lh, dimL_orb, tau, vr, &
     &                  dimL_orb, work4, lwork4, info)

         WRITE(6,*)
         WRITE(6,*) 'info_lr_ev_2: ', info, 'work(1): ', work4(1), &
     &             'lwork: ', lwork4
         WRITE(6,*)

!     Check norm and orthogonalization

         CALL norm_orb(PSI,vr,norm_mat,orth_mat,indorder2, &
     &                   u,v)

!     Write

         DO j=1,numeig
            DO q=1,Morb
               DO i=1,ND
                  WRITE(22,'(E16.8,1X,I4,1X,I4,1X,2E16.8,1X,I2)') &
     &                     ort_x(i), j, q, u(i,j,q)
                  WRITE(23,'(E16.8,1X,I4,1X,I4,1X,2E16.8,1X,I2)') &
     &                     ort_x(i), j, q, v(i,j,q)
               END DO
            END DO
         END DO
        
         DO i=1,numeig
            DO j=1,numeig
               WRITE(28,'(I8,1X,I8,1X,2E16.8,1X,2E16.8)') i, j, &
     &                          norm_mat(i,j), orth_mat(i,j)
            END DO
         END DO

!        CLOSE(20-23); 
        CLOSE(20); 
        CLOSE(21); 
        CLOSE(22); 
        CLOSE(23); 
        CLOSE(28);

      END SUBROUTINE diag_orb_matrix

!ccccccccccccccccccccccccccccccccccc
!     Sort eigenvalues
!ccccccccccccccccccccccccccccccccccc

      SUBROUTINE SHELL(N,ARR,ind)
        
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: N
         INTEGER :: LOGNB2, m, nn, j, k, l, i
         REAL*8 :: t
         INTEGER :: ti
         REAL*8, parameter :: ALN2I=1./0.69314718,TINY=1.E-5
         REAL*8, DIMENSION(N), INTENT(INOUT) :: ARR 
         INTEGER, DIMENSION(N), INTENT(INOUT) :: ind 
         LOGNB2=INT(ALOG(FLOAT(N))*ALN2I+TINY)
         m=n
         do nn=1,LOGNB2
          m=m/2; k=n-m
           do j=1,k
           i=j
10    continue
       l=i+m
       if(ARR(l).LT.ARR(i)) then
        t=ARR(i); ti=ind(i)
         ARR(i)=ARR(l); ind(i)=ind(l)
         ARR(l)=t; ind(l)=ti
         i=i-m
         if(i.GE.1) GOTO 10
      end if
      end do
      end do
      return
      END SUBROUTINE SHELL

!=========================================
!===Produce C_i_j and C_i_j_k_l
!=========================================

       subroutine Produce_Cij(VIN,VOUT,i,j,k,l,Nc)
       USE CI_All
       USE Parallel_CI
       USE omp_lib
       USE CI_prod
       implicit NONE
!=================== MPI 
       INCLUDE 'mpif.h'
       INTEGER ::  ierr,MYID,IPRC,numprocs
!==========================================================
       COMPLEX*16, INTENT(IN) :: VIN(Nc)
       COMPLEX*16, INTENT(OUT) :: VOUT(Nc)
!==========================================================
       integer :: nn,Iorb,Jorb,ihelp,ihelphelp,jhelp,jhelphelp,khelp,lhelp,conj,suml,sumr
       integer :: maxsil,NL,Iter_Total,Nterms
       INTEGER :: II, KK, PP, cI, cJ, cK, cL, I_current_term
       INTEGER :: myII, myPP, mycI, mycJ, mycK, mycL
       real*8 :: xbprefac
       INTEGER, INTENT(IN) :: i, j, k, l, Nc

!=====1-body===============================================
       IF(k==0.OR.l==0) THEN

!         discriminate ij and ji

          conj=0
          ihelp=i; jhelp=j;
          IF(i>j) THEN
             ihelp=j; jhelp=i;
             conj=1
          END IF

!         create index for mapping

!          II=1
!          DO Iorb=1,Morb
!             DO Jorb=Iorb,Morb
!                IF(ihelp==Iorb.and.jhelp==Jorb) GOTO 119
!                II=II+1
!             EndDO      
!          EndDO    

!         Find index II from ihelp and jhelp
                VOUT=0d0 !Added by str
!              write(6,*)"Im in Produce _CI"
          DO II=1,MaxTrm1b
             I_current_term= TERM_REQ_1B(II)
             PP=TERM_INDEX_1B(I_current_term)
             cJ= INT(PP/100); cI=PP-cJ*100; 
             IF(cI==ihelp.AND.cJ==jhelp) THEN
              myII=II; MyPP=PP; mycI=cI; mycJ=cJ;
             END IF
          END DO
        II=myII; PP=MyPP; cI=mycI; cJ=mycJ;
!              write(6,*)"Im in Produce II, PP ",II,PP
!              write(6,*)"Im in Produce I,J",cI,CJ
!              write(6,*)"Prefactors CI 1 1",Prefactors_1b(:,II)
!              write(6,*)"Im in Produce conj",conj
!STRSTR       WRITE(6,*) 'II, ij', II, 'cI:', cI, 'cJ:', cJ, 'I_curr', I_current_term

!         loop over configurations

          IF(conj==0) THEN
             Do KK=1,Nc 
                PP=Ind_CI_1b(KK,II)
                xbprefac=Prefactors_1b(KK,II)
                VOUT(PP)=VOUT(PP)+VIN(KK)* DSQRT(xbprefac) ! STR corrected
!                write(6,*)" qq", KK, "v",PP," times ",xbprefac
!                write(6,*)"Vin", VIN(KK), " Vout", VOUT(PP)
!                pause
             EndDO       
          ELSE 
             Do KK=1,Nc
                PP=Ind_CI_1b(KK,II)
                xbprefac=Prefactors_1b(KK,II)
                VOUT(KK)=VOUT(KK)+VIN(PP)* DSQRT(xbprefac)
             EndDO     
          END IF
!              write(6,*)"INd CI 1 1",Ind_CI_1b(:,1)
!              write(6,*)"VIN",VIN(:)
!              write(6,*)"VOUT",VOUT(:)
!              write(6,*)"Prefactors CI 1 1",Prefactors_1b(:,3)
!=====2-body===============================================
       ELSE
     
!         discriminate ijkl: ij,ji and kl,lk

          ihelp=i; jhelp=j;
          IF(i>j) THEN
             ihelp=j; jhelp=i;
          END IF
          khelp=k; lhelp=l;
          IF(k>l) THEN
             khelp=l; lhelp=k;
          END IF

!         discriminate ijkl: ij and kl

          suml=i+j; sumr=k+l;
          conj=0
          IF(ihelp>khelp ) THEN
             ihelphelp=ihelp; jhelphelp=jhelp
             ihelp=khelp; jhelp=lhelp; khelp=ihelphelp; lhelp=jhelphelp;
             conj=1
          END IF
          IF(ihelp==khelp .AND. jhelp>lhelp ) THEN
             ihelphelp=ihelp; jhelphelp=jhelp
             ihelp=khelp; jhelp=lhelp; khelp=ihelphelp; lhelp=jhelphelp;
             conj=1
          END IF

!         create index for mapping
          DO II=1,MaxTrm2b
             I_current_term= TERM_REQ_2B(II)
             PP=TERM_INDEX_2B(I_current_term)
             cL= INT(PP/1000000)
             cK= INT((PP-cL*1000000)/10000)
             cJ= INT((PP-cL*1000000-cK*10000)/100)
             cI= PP-cL*1000000-cK*10000-cJ*100
             IF(cI==ihelp.AND.cJ==jhelp.AND.cK==khelp.AND.cL==lhelp) THEN
              myII=II; MyPP=PP; mycI=cI; mycJ=cJ; mycK=cK; mycL=cL
             END IF
          END DO
        II=myII; PP=MyPP; cI=mycI; cJ=mycJ; cK=mycK; cL=mycL
!STR STR   WRITE(6,*) 'II',II,'PP',PP, 'cI', cI, 'cJ', cJ, 'cK', cK, 'cL', cL, 'conj', conj

!         loop over configurations

          IF(conj==0) THEN
             Do KK=1,Nc 
                PP=Ind_CI_2b(KK,II)
                xbprefac=Prefactors_2b(KK,II)
                VOUT(PP)=VOUT(PP)+VIN(KK)* DSQRT(xbprefac) !STR
             EndDO       
          ELSE 
             Do KK=1,Nc
                PP=Ind_CI_2b(KK,II)
                xbprefac=Prefactors_2b(KK,II)
                VOUT(KK)=VOUT(KK)+VIN(PP)* DSQRT(xbprefac) ! STR
             EndDO     
          END IF

       END IF

       end subroutine Produce_Cij

!==========================================================
!===CONSTRUCT CO- MATRIX==============================================================
!==========================================================

       SUBROUTINE construct_CO_matrix(PSI,VIN,COmat,COmat_bare,Nc)

          USE CI_ALL
          USE Parallel_CI 
       USE   W_INTERPARTICLE !ADD BY STR
          IMPLICIT NONE
          COMPLEX*16, INTENT(IN) :: VIN(Nc)
          COMPLEX*16, INTENT(OUT) :: COmat(2*Nc,2*ND*Morb)
          COMPLEX*16, INTENT(OUT) :: COmat_bare(2*Nc,2*ND*Morb)
          COMPLEX*16  :: PL(2*Nc,2*ND*Morb)
          COMPLEX*16, DIMENSION(ND,Morb), INTENT(IN) :: PSI
          COMPLEX*16 :: Vhelp1(Nc), Vhelp2(Nc)
          COMPLEX*16 :: Vhelp1a(Nc), Vhelp2a(Nc)
          COMPLEX*16 :: Vhelp1b(Nc), Vhelp2b(Nc)
          COMPLEX*16, DIMENSION(ND,Morb) :: h2_PSI
          COMPLEX*16 :: COmat_help(2*Nc,2*ND*Morb)
          COMPLEX*16, DIMENSION(Morb,Morb) :: RhoInvSQR
          INTEGER, DIMENSION(2) :: ind11, ind12, ind21, ind22
          INTEGER, DIMENSION(2) :: ind11a, ind12a, ind21a, ind22a
          COMPLEX*16, DIMENSION(dimL_orb,dimL_orb) :: Pb, Pbb, PLL
          COMPLEX*16, DIMENSION(dimL_orb/2,dimL_orb/2) :: P
          COMPLEX*16 :: VIN_help(Nc)
          COMPLEX*16 :: rhohelp
         COMPLEX*16, DIMENSION(NDX*NDY*NDZ)  :: WSL !ADDED BY STR
!==========================================================
          INTEGER :: i,j,k,s,l,q,Nc
!===========================================================
           xlambda0=xlambda_0

!         Initialize

          COmat=0.0d0; Vhelp1a=0.0d0; Vhelp2a=0.0d0; Vhelp1b=0.0d0; Vhelp2b=0.0d0
          CALL get_h2_PSI(PSI,h2_PSI)

!         upper left 
!                      OPEN(66,FILE='Vhelp1.out')
!STR                      OPEN(68,FILE='Vhelp2.out')
!STR                      OPEN(67,FILE='COmat.out')
         
          DO q=1,Morb
             DO k=1,Morb

!                      IF(k==q) THEN
!               1-body
                VIN_help=VIN; Vhelp1a=0.0d0
                CALL Produce_Cij(VIN_help,Vhelp1a,k,q,0,0,Nc) !Correct - proved by STR
                rhohelp=SUM(Conjg(VIN_help)*Vhelp1a)
!STR                WRITE(66,'(I4,2X,I4,2X,2E16.8,2E16.8)') k, q, rhohelp, Rho_ij(k,q)

                DO i=1,Nc
                   DO j=1,ND
                      COmat(i,(q-1)*ND+j)=COmat(i,(q-1)*ND+j)+Vhelp1a(i)*CONJG(h2_PSI(j,k))
                   END DO
                END DO
!                      END IF
!               2-body

                DO s=1,Morb
                   DO l=1,Morb
!                      IF(k==l.AND.s==k) THEN
!                      IF(q==s.AND.k==l.AND.q==k) GOTO 111
                      VIN_help=VIN; Vhelp2a=0.0d0
                      CALL Produce_Cij(VIN_help,Vhelp2a,k,s,q,l,Nc) !Correct - proved by STR
                rhohelp=SUM(Conjg(VIN)*Vhelp2a)
!STR                WRITE(68,'(I4,2X,I4,2X,I4,2X,I4,2X,2E16.8,2E16.8)') k,s,q,l, rhohelp,Rho_ijkl(k,s,q,l)
                 
!STR add          
            CALL  Get_WSL_omp(WSL,psi(:,s),psi(:,l))
                      DO i=1,Nc
                         DO j=1,ND
             
!                  COmat(i,(q-1)*ND+j)=COmat(i,(q-1)*ND+j)+xlambda_0*Vhelp2a(i)*CONJG(PSI(j,k))*CONJG(PSI(j,s))*PSI(j,l)
                            COmat(i,(q-1)*ND+j)=COmat(i,(q-1)*ND+j)+Vhelp2a(i)*CONJG(PSI(j,k))*WSL(j)
                             
                         END DO
                      END DO
!111           CONTINUE
!                     END IF
                   END DO
                END DO

             END DO

!               lower right

             DO i=1,Nc
                DO j=1,ND
                   COmat(Nc+i,Morb*ND+(q-1)*ND+j)=-CONJG(COmat(i,(q-1)*ND+j))
                END DO
             END DO

          END DO

!         upper right
          CALL get_h2_PSI(CONJG(PSI),h2_PSI)

          DO q=1,Morb
             DO k=1,Morb

!                      IF(k==q) THEN
!               1-body
                VIN_help=VIN; Vhelp1b=0.0d0
                CALL Produce_Cij(VIN_help,Vhelp1b,q,k,0,0,Nc) !Correct - proved by STR
                rhohelp=SUM(Conjg(VIN)*Vhelp1b)
   !STR             WRITE(66,'(I4,2X,I4,2X,2E16.8)') k, q, rhohelp 
!            write(6,*)"Produce Cij", k,q, rhohelp
!            write(6,*)"Vin", VIN
!            write(6,*)"Vout",Vhelp1b

!               pause

                DO i=1,Nc
                   DO j=1,ND
             
                      COmat(i,Morb*ND+(q-1)*ND+j)=COmat(i,Morb*ND+(q-1)*ND+j)+Vhelp1b(i)*CONJG(h2_PSI(j,k))
                   END DO
                END DO
!                       END IF
!               2-body

                DO s=1,Morb
                   DO l=1,Morb
!                      IF(k==l.AND.s==k) THEN
!                      IF(q==s.AND.k==l.AND.q==k) GOTO 112
                      VIN_help=VIN; Vhelp2b=0.0d0
                      CALL Produce_Cij(VIN_help,Vhelp2b,q,s,k,l,Nc)
                      rhohelp=SUM(Conjg(VIN)*Vhelp2b)
!STR                      WRITE(68,'(I4,2X,I4,2X,I4,2X,I4,2X,2E16.8)') q,s,k,l, rhohelp
            
!STR add          
            CALL  Get_WSL_omp(WSL,psi(:,s),psi(:,l))
                      DO i=1,Nc
                         DO j=1,ND
             
                        COmat(i,Morb*ND+(q-1)*ND+j)=COmat(i,Morb*ND+(q-1)*ND+j)+Vhelp2b(i)*PSI(j,k)*WSL(j)
!ORG                            COmat(i,Morb*ND+(q-1)*ND+j)=COmat(i,Morb*ND+(q-1)*ND+j)+xlambda_0*Vhelp2b(i)*PSI(j,k)*CONJG(PSI(j,s))*PSI(j,l)

                         END DO
                      END DO
!112       CONTINUE
!                      END IF
                   END DO
                END DO

             END DO

!               lower left

             DO i=1,Nc
                DO j=1,ND
                   COmat(Nc+i,(q-1)*ND+j)=-CONJG(COmat(i,Morb*ND+(q-1)*ND+j))
                END DO
             END DO

          END DO

          COmat_bare=COmat
          COmat=COmat

!           DO i=1,Nc
!              DO q=1,Morb
!                 DO j=1,ND
!
!                    WRITE(67,'(I4,2X,I4,2X,2E16.8)') i, (q-1)*ND+j, COmat(i,(q-1)*Morb+j)
!
!                 END DO
!              END DO
!          END DO

!     Get projector matrix

         CALL get_proj(PSI,Pb,P)

!     Multiply proj. matrix (double)

         CALL ZGEMM('N', 'N', dimL_orb, dimL_orb, dimL_orb, DCMPLX(1.0d0), Pb, &
     &             dimL_orb, Pb, dimL_orb, DCMPLX(0.0d0), Pbb, dimL_orb)

         CALL ZGEMM('N', 'N', 2*Nc, dimL_orb, dimL_orb, DCMPLX(1.0d0), COmat, &
     &             2*Nc, Pbb, dimL_orb, DCMPLX(0.0d0), PL, 2*Nc)

         COmat=PL

!      metric
!      get squareroot of inverse of rho

          CALL squarerootRInv(RhoInvSQR)

!      before
          COmat_help=0.0d0;
          DO s=1,Nc
             DO k=1,Morb
                DO i=1,ND
                   DO j=1,Morb

                  ! for the blocks
                  ind11=(/s,(j-1)*ND+i/)
                  ind12=(/s,Morb*ND+(j-1)*ND+i/)
                  ind21=(/Nc+s,(j-1)*ND+i/)
                  ind22=(/Nc+s,Morb*ND+(j-1)*ND+i/)
                  ind11a=(/s,(k-1)*ND+i/)
                  ind12a=(/s,Morb*ND+(k-1)*ND+i/)
                  ind21a=(/Nc+s,(k-1)*ND+i/)
                  ind22a=(/ Nc+s,Morb*ND+(k-1)*ND+i/)

                  COmat_help(ind11a(1),ind11a(2))=COmat_help(ind11a(1),ind11a(2))+COmat(ind11(1),ind11(2))*RhoInvSQR(j,k)

                  COmat_help(ind12a(1),ind12a(2))=COmat_help(ind12a(1),ind12a(2))+COmat(ind12(1),ind12(2))*CONJG(RhoInvSQR(j,k))

                  COmat_help(ind21a(1),ind21a(2))=COmat_help(ind21a(1),ind21a(2))+COmat(ind21(1),ind21(2))*RhoInvSQR(j,k)

                  COmat_help(ind22a(1),ind22a(2))=COmat_help(ind22a(1),ind22a(2))+COmat(ind22(1),ind22(2))*CONJG(RhoInvSQR(j,k))
                  END DO
               END DO
            END DO
         END DO

         COmat=COmat_help


!           DO i=1,Nc
!              DO q=1,Morb
!                 DO j=1,ND
!
!                    WRITE(67,'(I4,2X,I4,2X,2E16.8)') i, (q-1)*ND+j, COmat(i,(q-1)*Morb+j)
!
!                 END DO
!              END DO
!          END DO

!         CLOSE(66-68)

       END SUBROUTINE construct_CO_matrix

!==========================================================
!===CONSTRUCT OC- MATRIX==============================================================
!==========================================================

       SUBROUTINE construct_OC_matrix(PSI,VIN,OCmat,Nc)

          USE CI_ALL
          USE Parallel_CI 
          IMPLICIT NONE
          COMPLEX*16, INTENT(IN) :: VIN(Nc)
          COMPLEX*16, INTENT(OUT) :: OCmat(2*ND*Morb,2*Nc)
          COMPLEX*16, DIMENSION(ND,Morb), INTENT(IN) :: PSI

          COMPLEX*16 :: Vhelp1(Nc), Vhelp2(Nc)
          COMPLEX*16, DIMENSION(ND,Morb) :: h2_PSI
          COMPLEX*16 :: COmat(2*Nc,2*ND*Morb)
          COMPLEX*16 :: COmat_bare(2*Nc,2*ND*Morb)
          COMPLEX*16 :: PL(2*ND*Morb,2*Nc)
          COMPLEX*16, DIMENSION(Morb,Morb) :: RhoInvSQR
          COMPLEX*16 :: OCmat_help(2*ND*Morb,2*Nc)
          INTEGER, DIMENSION(2) :: ind11, ind12, ind21, ind22
          INTEGER, DIMENSION(2) :: ind11a, ind12a, ind21a, ind22a
          COMPLEX*16, DIMENSION(dimL_orb,dimL_orb) :: Pb, Pbb, PLL
          COMPLEX*16, DIMENSION(dimL_orb/2,dimL_orb/2) :: P
!==========================================================
          INTEGER :: i,j,k,s,l,q,Nc
!===========================================================

!         Initialize: COmat

          OCmat=0.0d0
          CALL construct_CO_matrix(PSI,VIN,COmat,COmat_bare,Nc)


          DO q=1,Morb
             DO i=1,Nc
                DO j=1,ND
             
!         upper left 
                OCmat((q-1)*ND+j,i)=CONJG(COmat_bare(i,(q-1)*ND+j))

!         lower right 
                OCmat(Morb*ND+(q-1)*ND+j,Nc+i)=CONJG(COmat_bare(Nc+i,Morb*ND+(q-1)*ND+j))

!         upper right 
                OCmat((q-1)*ND+j,Nc+i)=-CONJG(COmat_bare(Nc+i,(q-1)*ND+j))

!         lower left 
                OCmat(ND*Morb+(q-1)*ND+j,i)=-CONJG(COmat_bare(i,Morb*ND+(q-1)*ND+j))

                END DO
             END DO
          END DO

!     Get projector matrix

         CALL get_proj(PSI,Pb,P)

!     Multiply proj. matrix (double)

         CALL ZGEMM('N', 'N', dimL_orb, dimL_orb, dimL_orb, DCMPLX(1.0d0), Pb, &
     &             dimL_orb, Pb, dimL_orb, DCMPLX(0.0d0), Pbb, dimL_orb)

         CALL ZGEMM('N', 'N', dimL_orb, 2*Nc, dimL_orb, DCMPLX(1.0d0), Pbb, &
     &             dimL_orb, OCmat, dimL_orb, DCMPLX(0.0d0), PL, dimL_orb)

         OCmat=PL

!      metric
!      get squareroot of inverse of rho

          CALL squarerootRInv(RhoInvSQR)

!      before
          OCmat_help=0.0d0;
          DO k=1,Morb
             DO s=1,Nc
                DO i=1,ND
                   DO j=1,Morb

                  ! for the blocks
                  ind11=(/(j-1)*ND+i,s/)
                  ind12=(/(j-1)*ND+i,Nc+s/)
                  ind21=(/Morb*ND+(j-1)*ND+i,s/)
                  ind22=(/Morb*ND+(j-1)*ND+i,Nc+s/)
                  ind11a=(/(k-1)*ND+i,s/)
                  ind12a=(/(k-1)*ND+i,Nc+s/)
                  ind21a=(/Morb*ND+(k-1)*ND+i,s/)
                  ind22a=(/Morb*ND+(k-1)*ND+i,Nc+s/)

                  OCmat_help(ind11a(1),ind11a(2))=OCmat_help(ind11a(1),ind11a(2))+RhoInvSQR(k,j)*OCmat(ind11(1),ind11(2))

                  OCmat_help(ind12a(1),ind12a(2))=OCmat_help(ind12a(1),ind12a(2))+RhoInvSQR(k,j)*OCmat(ind12(1),ind12(2))

                  OCmat_help(ind21a(1),ind21a(2))=OCmat_help(ind21a(1),ind21a(2))+CONJG(RhoInvSQR(k,j))*OCmat(ind21(1),ind21(2))

                  OCmat_help(ind22a(1),ind22a(2))=OCmat_help(ind22a(1),ind22a(2))+CONJG(RhoInvSQR(k,j))*OCmat(ind22(1),ind22(2))

                  END DO
               END DO
            END DO
         END DO

         OCmat=OCmat_help
!        write(6,*)"OCmat_NEW", SUM(OCmat),SUM(OCmat(:,1)),SUM(OCmat(:,2)),SUM(OCmat(:,100))
!         pause
       END SUBROUTINE construct_OC_matrix

!==========action of h2 on PSI      
       SUBROUTINE get_h2_PSI(PSI,h2_PSI)

          USE CI_ALL
          USE Parallel_CI 
          USE DVR_ALL

          IMPLICIT NONE

          COMPLEX*16, DIMENSION(ND,Morb), INTENT(IN) :: PSI
          COMPLEX*16, DIMENSION(ND,Morb), INTENT(OUT) :: h2_PSI

          COMPLEX*16, DIMENSION(dimL_orb,dimL_orb) :: L_orb
          COMPLEX*16, DIMENSION(dimL_orb/2,dimL_orb/2) :: h2

          INTEGER :: i,j,q
!==========================================================================

          h2_PSI=0.0d0
          CALL h2_help(h2)

          DO q=1,Morb
             DO i=1,ND
                DO j=1,ND
      
                   h2_PSI(i,q)=h2_PSI(i,q)+h2(i,j)*PSI(j,q) 
           
                END DO
             END DO
          END DO

       END SUBROUTINE get_h2_PSI

!==========================================================================
!==========squareroot of inverse one-body reduced density matrix
!==========================================================================

       SUBROUTINE squarerootRInv(RhoInvSQR)

          USE CI_ALL
          USE Parallel_CI 
          USE DVR_ALL

          IMPLICIT NONE         

          COMPLEX*16, DIMENSION(Morb,Morb), INTENT(OUT) :: RhoInvSQR

          COMPLEX*16, DIMENSION(Morb,Morb) :: V, VI, Vhelp
          REAL*8, DIMENSION(Morb) :: w
         
          INTEGER :: i, j, k, info, infoLU, infoInv
          INTEGER, PARAMETER :: lwork=200, lrwork=500, liwork=100, lworkInv=20
          COMPLEX*16, DIMENSION(lwork) :: work
          REAL*8, DIMENSION(lrwork) :: rwork
          INTEGER, DIMENSION(liwork) :: iwork
          COMPLEX*16, DIMENSION(lworkInv) :: workInv
          INTEGER, DIMENSION(Morb) :: ipiv 

!     Eigensystem of InvZRIJ: w, V 

            
          V=InvZRIJ
          CALL ZHEEVD('V','U', Morb, V, Morb, w, work, lwork, rwork, lrwork, iwork, &
     &       liwork, info)

          if(info.ne.0) WRITE(6,*) 'info ZHEEVD in squarerootRInv:', info, 'work(1): ', & 
     &              work(1), 'lwork: ', lwork, 'rwork(1): ', rwork(1), &
     &             'lrwork: ', lrwork, 'iwork(1): ', iwork(1), &
     &             'liwork: ', liwork

!          WRITE(6,*)

!     Inverse of V
!     LU
          VI=V

          CALL ZGETRF( Morb, Morb, VI, Morb, ipiv, infoLU )

          if(info.ne.0) WRITE(6,*) 'info ZGETRF in squarerootRInv:', infoLU
           
          WRITE(6,*) 
!     INV
          CALL ZGETRI( Morb, VI, Morb, ipiv, workInv, lworkInv, infoInv )

           if(info.ne.0) WRITE(6,*) 'info ZGETRI in squarerootRInv:', infoInv, 'work(1):', workInv(1)

!          WRITE(6,*)

!     Calculate squareroot of inverse matrix

!STR2014          OPEN(91,file='rhoinvsqr.out')

          RhoInvSQR=0.0d0
          DO i=1,Morb
             DO j=1,Morb
                Vhelp(i,j)=V(i,j)*SQRT(w(j))           
             END DO
          END DO

         CALL ZGEMM('N', 'N', Morb, Morb, Morb, DCMPLX(1.0d0),Vhelp, &
     &                      Morb, VI, Morb, DCMPLX(0.0d0), RhoInvSQR, Morb)
 
!STR2014          DO i=1,Morb
!STR2014             DO j=1,Morb
!STR2014                WRITE(91,'(I4,I4,2E16.8)') i,j,RhoInvSQR(i,j)
!STR2014             END DO
!STR2014          END DO

!STR2014          CLOSE(91)

       END SUBROUTINE squarerootRInv

!==========================================================================
!==========squareroot of one-body reduced density matrix
!==========================================================================

       SUBROUTINE squarerootR(RhoSQR)

          USE CI_ALL
          USE Parallel_CI 
          USE DVR_ALL

          IMPLICIT NONE         

          COMPLEX*16, DIMENSION(Morb,Morb), INTENT(OUT) :: RhoSQR

          COMPLEX*16, DIMENSION(Morb,Morb) :: V, VI, Vhelp
          REAL*8, DIMENSION(Morb) :: w
         
          INTEGER :: i, j, k, info, infoLU, infoInv
          INTEGER, PARAMETER :: lwork=200, lrwork=500, liwork=100, lworkInv=20
          COMPLEX*16, DIMENSION(lwork) :: work
          REAL*8, DIMENSION(lrwork) :: rwork
          INTEGER, DIMENSION(liwork) :: iwork
          COMPLEX*16, DIMENSION(lworkInv) :: workInv
          INTEGER, DIMENSION(Morb) :: ipiv 

!     Eigensystem of InvZRIJ: w, V 

          V=AllZRIJ
          CALL ZHEEVD('V','U', Morb, V, Morb, w, work, lwork, rwork, lrwork, iwork, &
     &       liwork, info)

          if(info.ne.0) WRITE(6,*) 'info ZHEEVD in squarerootR:', info, 'work(1): ', & 
     &              work(1), 'lwork: ', lwork, 'rwork(1): ', rwork(1), &
     &             'lrwork: ', lrwork, 'iwork(1): ', iwork(1), &
     &             'liwork: ', liwork

          WRITE(6,*)"Allrij NO",w
           Do i=1,Morb
           if(w(i).lt.0d0) w(i)=abs(w(i))
           enddo
         IF(info.ne.0) write(6,*) "Problem in squarerootR ZHEEVD"

!     Inverse of V
!     LU
          VI=V

          CALL ZGETRF( Morb, Morb, VI, Morb, ipiv, infoLU )

!          WRITE(6,*) 'info ZGETRF in squarerootR:', infoLU
         IF(infoLU.ne.0) write(6,*) "Problem in squarerootR in ZGETRF", infoLU
           
!          WRITE(6,*) 
!     INV
          CALL ZGETRI( Morb, VI, Morb, ipiv, workInv, lworkInv, infoInv )

!          WRITE(6,*) 'info ZGETRI in squarerootR:', infoInv, 'work(1):', workInv(1)
         IF(infoInv.ne.0) write(6,*) "Problem in squarerootR in ZGETRI",infoInv

!          WRITE(6,*)

!     Calculate squareroot of inverse matrix

!STR2014          OPEN(92,file='rhosqr.out')

          RhoSQR=0.0d0
          DO i=1,Morb
             DO j=1,Morb
                Vhelp(i,j)=V(i,j)*SQRT(w(j))           
             END DO
          END DO

         CALL ZGEMM('N', 'N', Morb, Morb, Morb, DCMPLX(1.0d0),Vhelp, &
     &                      Morb, VI, Morb, DCMPLX(0.0d0), RhoSQR, Morb)
 
!STR2014          DO i=1,Morb
!STR2014             DO j=1,Morb
!STR2014                WRITE(92,'(I4,I4,2E16.8)') i,j,RhoSQR(i,j)
!STR2014             END DO
!STR2014          END DO
!STR2014          CLOSE(92)

       END SUBROUTINE squarerootR

!==========================================================================
!==========Construct the whole LR-MCTDHB matrix           
!==========================================================================

       SUBROUTINE construct_matrix(PSI,VIN,L,Nc,pathLR)

          USE CI_ALL
          USE Parallel_CI 
          USE DVR_ALL
          IMPLICIT NONE
          COMPLEX*16, INTENT(IN) :: VIN(Nc)
          COMPLEX*16, DIMENSION(ND,Morb), INTENT(IN) :: PSI
          COMPLEX*16, INTENT(OUT) :: L(dimL,dimL)

          COMPLEX*16 :: L_orb(dimL_orb,dimL_orb)
          COMPLEX*16 :: OCmat(dimL_orb,2*Nc)
          COMPLEX*16 :: COmat(2*Nc,dimL_orb)
          COMPLEX*16 :: COmat_bare(2*Nc,dimL_orb)
          COMPLEX*16 :: CImat(Nc,Nc)
!          COMPLEX*16, DIMENSION(dimL,dimL) :: L1, L3

          COMPLEX*16 :: Vhelp1(Nc), Vhelp2(Nc)
          COMPLEX*16, DIMENSION(ND,Morb) :: h2_PSI
          COMPLEX*16, DIMENSION(dimL_orb/2,dimL_orb/2) :: h2
          INTEGER, DIMENSION(2) :: ind11, ind12, ind21, ind22
          INTEGER, DIMENSION(2) :: ind11a, ind12a, ind21a, ind22a
          REAL*8    ::   RCIeigenvals(Nc),RCIeigenvecs(Nc,Nc),exp_Val
          COMPLEX*16 :: Tkin(ND,ND)
!==========================================================
          INTEGER :: i,j,k,s,q,P,Nc
          character*10 pathLR
!===========================================================

          L=0.0d0

!         Get submatrices
          CALL construct_orb_matrix(PSI,L_orb,h2)
          IF(Morb.gt.1) THEN
             CALL diag_CI_ex(RCIeigenvecs,RCIeigenvals,Nc,pathLR)
             write(6,*)"Pure MCHB CI eigenproblem is solved see it in Pure_MCHB_CI_spectrum.out"
             CALL construct_CI_matrix(CImat,Nc)
             CALL construct_OC_matrix(PSI,VIN,OCmat,Nc)
             CALL construct_CO_matrix(PSI,VIN,COmat,COmat_bare,Nc)
          END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           Exp Values Of Operators x x*x etc
           Tkin=0d0
           Do I=1,ND
           Tkin(i,i)=Ort_x(i)
           ENDDO
!              Tkin=Op_X !STR October 2012 modifications
!              DO j=1,ND
!              Tkin(j,j)=Tkin(j,j)+REAL(VTRAP_EXT(j))
!              END DO
         CALL Get_Exp_Val_1bOp(PSI,Tkin,Rho_ij,exp_Val)
           write(6,*)" <PSI_GS|x|PSI_GS>=",exp_Val
           Tkin=0d0
           Do I=1,ND
           Tkin(i,i)=Ort_x(i)**2
           ENDDO
         CALL Get_Exp_Val_1bOp(PSI,Tkin,Rho_ij,exp_Val)
           write(6,*)" <PSI_GS|x^2|PSI_GS>=",exp_Val
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         orb part
          DO i=1,dimL_orb
             DO j=1,dimL_orb
                L(i,j)=L_orb(i,j)
             END DO
          END DO
         write(6,*)"Is OrbMTRX COMPLX or NOT?",SUM(real(L_orb)),SUM(abs(imag(L_orb)))
!             pause
 
     IF(Morb.gt.1) THEN

!         CI part
          DO i=1,Nc
             DO j=1,Nc
                L(dimL_orb+i,dimL_orb+j)=CImat(i,j)
                L(dimL_orb+Nc+i,dimL_orb+Nc+j)=-CONJG(CImat(i,j))
             END DO
          END DO

!         OC part
          DO i=1,dimL_orb
             DO j=1,2*Nc
                L(i,dimL_orb+j)=OCmat(i,j)
             END DO
          END DO

!         CO part
          DO i=1,2*Nc
             DO j=1,dimL_orb
!ORG                L(dimL_orb+i,j)=COmat(i,j)*dx*0
                L(dimL_orb+i,j)=COmat(i,j)
             END DO
          END DO

      END IF

!     Analyze symmetries
!     Calculate transformed matrices L1 and L3
             write(6,*)" Nc before ..",Nc
!          CALL analyse_matrix(L, L1, L3,Nc)

       END SUBROUTINE construct_matrix

!=================================================
!======Diagonalize full matrix====================
!=================================================

      SUBROUTINE diag_full_matrix(PSI,L,wo,u,v,Cu,Cv,Nc,&
                     norm_mat_orb,norm_mat_CI,pathLR)

          USE CI_ALL
          USE Parallel_CI 
          USE DVR_ALL

          IMPLICIT NONE         

          INTEGER :: i, j, k, q, s, ll, mm, j1
  
          COMPLEX*16, DIMENSION(ND,Morb), INTENT(IN) :: PSI
          COMPLEX*16, DIMENSION(dimL,dimL), INTENT(IN) :: L
          INTEGER, INTENT(IN) :: Nc
         
!          COMPLEX*16, DIMENSION(dimL), INTENT(OUT) :: wo_org 
          COMPLEX*16, DIMENSION(numeig), INTENT(OUT) :: wo 
          COMPLEX*16, DIMENSION(ND,numeig,Morb), INTENT(OUT) :: u, v
          COMPLEX*16, DIMENSION(Nc,numeig), INTENT(OUT) :: Cu, Cv

          COMPLEX*16, DIMENSION(dimL,dimL) :: Lh
          COMPLEX*16, DIMENSION(dimL) :: w
          REAL*8, DIMENSION(dimL) :: wh
          INTEGER, PARAMETER :: lwork=24576, lwork2=3201, lwork4=9600
 !        INTEGER, PARAMETER :: lwork=245760, lwork2=32010, lwork4=9600
!         INTEGER, PARAMETER :: lwork=2045760, lwork2=302010, lwork4=90600
         INTEGER :: info
         COMPLEX*16 :: z
         COMPLEX*16, DIMENSION(dimL-1) :: tau
         INTEGER, DIMENSION(dimL) :: ind
         COMPLEX*16, DIMENSION(lwork) :: work
         COMPLEX*16, DIMENSION(lwork2) :: work2
         COMPLEX*16, DIMENSION(lwork4) :: work4
         REAL*8, DIMENSION(dimL) :: rwork3
         COMPLEX*16, DIMENSION(dimL*dimL) :: work3
 
         LOGICAL, DIMENSION(dimL) :: select
         INTEGER, DIMENSION(numeig) :: ifaill, ifailr
         COMPLEX*16 :: vl
         COMPLEX*16, DIMENSION(dimL,numeig) :: vr
         INTEGER :: mn,problematicStates(numeig)

         COMPLEX*16, DIMENSION(numeig,numeig) :: norm_mat, orth_mat
         COMPLEX*16, DIMENSION(numeig) :: norm_mat_orb, norm_mat_CI
         INTEGER, DIMENSION(numeig) :: indvec, indorder, indorder2
         REAL*8, DIMENSION(numeig) :: indvecr
         INTEGER :: nz,isft
         character*10 pathLR
        
         
         OPEN (20,FILE=trim(pathLR)//path_sep//'MC_w.out')
         OPEN (21,FILE=trim(pathLR)//path_sep//'MC_wtest.out')
         OPEN (19,FILE=trim(pathLR)//path_sep//'MC_wstr.out')
         OPEN (22,FILE=trim(pathLR)//path_sep//'MC_ev_u.out')
         OPEN (23,FILE=trim(pathLR)//path_sep//'MC_ev_v.out')
         OPEN (24,FILE=trim(pathLR)//path_sep//'MC_ev_Cu.out')
         OPEN (25,FILE=trim(pathLR)//path_sep//'MC_ev_Cv.out')
         OPEN (28,FILE=trim(pathLR)//path_sep//'MC_norm.out')
         OPEN (29,FILE=trim(pathLR)//path_sep//'MC_norm_diff.out')
         print *, "Diagonalization of LR matrix - can take some time..."
!     Diagonalize matrix (eigenvalues)
!              zgehrd (N, ILO, IHI, A, LDA, TAU, WORK, LWORK, INFO)
!         CALL ZGEHRD(dimL, 1, dimL, L, dimL, tau, work, lwork, info)
!              zgehd2 (N, ILO, IHI, A, LDA, TAU, WORK, INFO) 
         CALL ZGEHD2(dimL, 1, dimL, L, dimL, tau, work, info)

         Lh=L  ! save projected matrix for later (eigenvector calculation) (upper Hessenberg form)


         if(info.ne.0) WRITE(6,*) 'Problem in lr_ZGEHD2: ', info
!         WRITE(6,*) 'info_lr_diag: ', info, 'work(1): ', work(1), & 
!     &             'lwork: ', lwork
!         WRITE(6,*) 'dimL', dimL
       
 
         CALL ZHSEQR('E', 'N', dimL, 1, dimL, L, dimL, w, z, 1, &
     &                work2, lwork2, info)
         
!         WRITE(6,*)
         if(info.ne.0) WRITE(6,*) 'Problem in lr_ZHSEQR:', info
!         WRITE(6,*) 'info_lr_diag_2": ', info, 'work(1): ', work2(1), &
!     &             'lwork: ', lwork2
         WRITE(6,*) "LR-MCTDHB eigenvalues are found see it all roots in MC_wtest.out"
!c         stop "done eigenvalues"
!     Sort eigenvalues

!         DO i=1,dimL
!            WRITE(21,'(I8,1X,2E16.8)') i, w(i)
!         END DO

         wh=REAL(w)

         DO i=1,dimL
            ind(i)=i
         END DO

         CALL SHELL(dimL,wh,ind) ! "ind" gives location of i-th value

!  Writing UNsorted eigenvalues to  (21,FILE='MC_wtest.out')

         DO i=1,dimL
         WRITE(21,'(I8,1X,F26.16,2E16.8,I8,1X,F26.16,2E16.8)')i,Real(w(i))+Energy,w(i),ind(i),Real(w(ind(i)))+Energy,w(ind(i))
         END DO

! Writing sorted eigenvalues to (20,FILE='MC_w.out')
       
         DO i=dimL/2+1,dimL
!         WRITE(20,'(I8,1X,I8,1X,2E16.8)') i-dimL/2, ind(i), w(ind(i))
         WRITE(20,'(I8,1X,I8,1X,2E16.8,1x,F26.16)') i-dimL/2, ind(i), w(ind(i)),Real(w(ind(i))+Energy)
!            WRITE(6,'(I8,1X,I8,1X,2E16.8)') i-dimL/2, ind(i), w(ind(i))
         END DO

         wo=0.0d0
         isft=200 !additional negative roots are analized
         if(Morb.eq.1) isft=0 ! for GP less roots
         DO i=dimL/2+1,dimL/2+numeig
         wo(i-dimL/2)=w(ind(i-isft))
         END DO
         indvec=ind(dimL/2+1-isft:dimL/2+numeig-isft) ! isft=0 for positive ones
!JL ORG         indvec=ind(dimL/2+1:dimL/2+numeig) ! positive ones
        write(6,*)"Anls of the roots from till:",dimL/2+1-isft,dimL/2+numeig-isft
          IF(isft.ne.0) write(6,*)isft," roots from -w manifold are considerd see line 2494"
         WRITE(6,*) "LR-MCTDHB eigenvalues are found see sorted roots in MC_w.out"
!     Compute selected eigenvectors
!STR            j=0
!STR            indvec=ind(dimL)
!STR            wo=w(ind(dimL))
!STR         DO i=1,dimL
!STR            IF(Real(w(ind(i))).GE.-1.000000000d0)  THEN
!STR            IF(Abs(w(ind(i))).GE.0.0000000001d0)  THEN
!STR            j=j+1 
!STR            IF(j.LE.numeig)  THEN
!STR            wo(j)=w(ind(i))
!STR            indvec(j)=ind(i)
!STR            ENDIF
!STR            ENDIF
!STR            ENDIF
!STR         END DO

!         DO i=1,numeig
!            WRITE(19,'(I8,1X,I8,F20.16,1X,F20.16)') i, indvec(i),Abs(wo(i)+Energy),Real(wo(i))
!         END DO
!     Define subset of eigenvectors to be computed
!          stop
!JL ORG         indvec=ind(dimL/2+1:dimL/2+numeig) ! positive ones
!         indvec=ind(1:mv) ! first ones
        
!     Compute order of subset: 

         DO i=1,numeig
            indorder(i)=i
            indorder2(i)=i
         END DO

         indvecr=REAL(indvec) ! subset of ind 731 732 730
         CALL SHELL(numeig,indvecr,indorder) ! gives the places of each rank 3 1 2  
         indvecr=REAL(indorder)  
         CALL SHELL(numeig,indvecr,indorder2) ! gives the rank of each place 2 3 1        
!     Calculate number of zero vectors
 
!         nz=0;
!         DO i=1,numeig
!            IF(ABS(wo(i))<1.0d-1) THEN
!               nz=nz+1;
!            END IF
!!               WRITE(6,*) 'nz', nz, ABS(wo(i))
!         END DO

!     Eigenvector calculation
  
         vl=0.0d0; vr=0.0d0
         select=.FALSE.
!           print *, indvec
         select(indvec)=.TRUE.
         mn=1
         work3=0.0d0; rwork3=0.0d0

         CALL ZHSEIN('R', 'Q', 'N', select, dimL, Lh, dimL, w, vl, & 
     &                  1, vr, dimL, numeig, mn, work3, rwork3, ifaill, &
     &                  ifailr, info)

  !       WRITE(6,*)
         if(info.ne.0) WRITE(6,*) 'Problem in info_lr_ev: ', info
 !        WRITE(6,*)

         CALL ZUNMHR('L', 'N', dimL, numeig, 1, dimL, Lh, dimL, tau, vr, &
     &                  dimL, work4, lwork4, info)

!        WRITE(6,*)
        if(info.ne.0) WRITE(6,*) 'Problem in info_lr_ev_2: ', info
!, 'work(1): ', work4(1), &
!    &             'lwork: ', lwork4
!        WRITE(6,*)

!     Check norm and orthogonalization

         CALL norm_full(PSI,vr,norm_mat,orth_mat,indorder2,u,v,Cu,Cv,Nc,norm_mat_orb,norm_mat_CI)

!     Write to (22,FILE='MC_ev_u.out') (23,FILE='MC_ev_v.out')

         DO j=1,numeig
            j1=ind(j+dimL/2-isft)
            j1=j
            DO q=1,Morb
               DO i=1,ND
                  WRITE(22,'(E16.8,1X,I4,1X,I4,1X,2E16.8,1X,I4,1X,F16.8)') &
     &                     ort_x(i), j, q, u(i,j,q),indvec(j), Real(wo(j))
                  WRITE(23,'(E16.8,1X,I4,1X,I4,1X,2E16.8,1X,I4,1X,F16.8)') &
     &                     ort_x(i), j, q, v(i,j,q),indvec(j), Real(wo(j))
               END DO
            END DO
         END DO
        
         DO j=1,numeig
            j1=ind(j+dimL/2-isft)
            DO q=1,Nc
                  WRITE(24,'(I4,1X,I4,1X,2F16.8)') &
     &                     q, j, Cu(q,j)
                  WRITE(25,'(I4,1X,I4,1X,2F16.8)') &
     &                     q, j, Cv(q,j)
            END DO
         END DO

! Writing          OPEN (24,FILE='MC_ev_Cu.out')
!         OPEN (25,FILE='MC_ev_Cv.out')
         DO i=1,numeig
            DO j=1,numeig
               WRITE(28,'(I8,1X,I8,1X,2E16.8,1X,2E16.8,1X)') i, j, &
     &      norm_mat(i,j), orth_mat(i,j) 
            END DO
         END DO

         DO i=1,numeig
        WRITE(29,'(I8,1X,F20.16,1X,2E16.8,1X,2E16.8,1X,F20.16)') i,abs(w(ind(i+dimL/2-isft))+Energy), &
     &                          norm_mat_orb(i), norm_mat_CI(i), &
     &                          abs(norm_mat_orb(i))-abs(norm_mat_CI(i))
         END DO

       CLOSE(20); 
       CLOSE(21); 
       CLOSE(22); 
       CLOSE(23); 
       CLOSE(24); 
       CLOSE(25); 
       CLOSE(28);
       CLOSE(29);
         DO i=1,numeig
       WRITE(19,'(I8,1X,I8,3F20.16)') i, indvec(i),Real(norm_mat(i,i)),&
         Real(w(ind(i+dimL/2-isft))+Energy),abs(w(ind(i+dimL/2-isft))+Energy)
!            WRITE(19,'(I8,1X,I8,3F20.16)') i, indvec(i),Real(norm_mat(i,i)),Abs(wo(i)+Energy),Real(wo(i))
         END DO
       CLOSE(19) !STR 
!STR   PRINT OUT problematic eigenalues ABS(e) around 0 and norm is negative
       i=0
       DO k=1,numeig
       IF(Real(norm_mat_orb(k)+norm_mat_CI(k)).le.0.0001d0) THEN
       i=i+1
       problematicStates(i)=k
!       norm_mat(:,k)=0d0
!       norm_mat(k,:)=0d0
!       orth_mat(k,:)=0d0
!       orth_mat(:,k)=0d0
       END IF
       END DO
       write(6,*)"Total number of the problematic States (probalbly from -w manifold) is:",i
       write(6,*)"Total number of the problematic States (probalbly from -w manifold) is:",i
       write(6,*)"They are",(problematicStates(k),k=1,i)
       write(6,'(80a)')"For analysis see file MC_anlsplot(ALL).out were for every required root you find"
       WRITE(6,'(50a,100a)') "#  i,  Re(w_i), Re(w_i+Eref), Re(Norm_Orb),",&
      "Re(Norm_CI),Abs(<PSI_ref|OP1(x)|LR_PSI_i>),Abs( <PSI_ref|OP2(x)|LR_PSI_i>)" 
!       write(6,*)"For other roots orthogonality  is=",Numeig-i ,SUM(ABS(norm_mat))
!       write(6,*)"For other roots orthonotmality is=",0, SUM(ABS(orth_mat))


      END SUBROUTINE diag_full_matrix

!=================================================
!======Evaluate====================
!=================================================

      SUBROUTINE evaluate(PSI,VIN,wo,u,v,Cu,Cv,Nc, &
                     norm_mat_orb,norm_mat_CI,pathLR)

          USE CI_ALL
          USE Parallel_CI 
          USE DVR_ALL
          USE rR_hW

          IMPLICIT NONE         

          INTEGER :: i, j, k, q, s, ll, mm
  
          COMPLEX*16, DIMENSION(ND,Morb), INTENT(IN) :: PSI
          COMPLEX*16, INTENT(IN) :: VIN(Nc)
          INTEGER, INTENT(IN) :: Nc
          COMPLEX*16, DIMENSION(numeig), INTENT(IN) :: wo 
          COMPLEX*16, DIMENSION(ND,numeig,Morb), INTENT(IN) :: u, v
          COMPLEX*16, DIMENSION(Nc,numeig), INTENT(IN) :: Cu, Cv
         COMPLEX*16, DIMENSION(numeig) :: norm_mat_orb, norm_mat_CI

         COMPLEX*16, DIMENSION(ND,numeig) :: dens_orb_x, dens_CI_x
         COMPLEX*16, DIMENSION(ND) :: dens_orb_ORG
         COMPLEX*16, DIMENSION(ND,numeig,Morb) :: ut, vt
         COMPLEX*16, DIMENSION(Morb,Morb) :: RhoInvSQR, RhoSQR
         COMPLEX*16, DIMENSION(ND,2) :: tr 
         COMPLEX*16, DIMENSION(Morb,Morb,numeig) :: CCu, CCu_r, CCv, CCv_r
         COMPLEX*16, DIMENSION(numeig) :: CCu_norm
         COMPLEX*16, DIMENSION(2) :: qqq
         COMPLEX*16 :: VIN_help(Nc), Vhelp(Nc)
         COMPLEX*16, DIMENSION(numeig,2) :: weight_orb, weight_CI
         COMPLEX*16, DIMENSION(Morb,Morb,2) :: orb_int
         character*10 pathLR

         OPEN (198,FILE=trim(pathLR)//path_sep//'MC_we.out')
         OPEN (12,FILE=trim(pathLR)//path_sep//'MC_anlsplot.out')
         OPEN (11,FILE=trim(pathLR)//path_sep//'MC_anlsplotALL.out')
         OPEN (199,FILE=trim(pathLR)//path_sep//'MC_dens_x.out')
         OPEN (197,FILE=trim(pathLR)//path_sep//'MC_rho.out')
         OPEN (196,FILE=trim(pathLR)//path_sep//'MC_we_help.out')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        density in position space (x)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        ORBITAL!!!!!!!!!!!!!!!!!!!!!!

!        first get metric

         CALL squarerootRInv(RhoInvSQR)

!        multiply metric with amplitudes

         ut=0.0d0; vt=0.0d0

         DO k=1,numeig
            DO j=1,Morb
               DO i=1,Morb

                  ut(:,k,j)=ut(:,k,j)+RhoInvSQR(j,i)*u(:,k,i)
                  vt(:,k,j)=vt(:,k,j)+RhoInvSQR(j,i)*v(:,k,i)

               END DO
            END DO
         END DO

!        now density

         dens_orb_x=0.0d0  
         DO k=1,numeig
         dens_orb_ORG=0.0d0  
            DO i=1,Morb
               DO j=1,Morb

                dens_orb_x(:,k)=dens_orb_x(:,k)+  &
     &            AllZRIJ(i,j)*PSI(:,i)*(ut(:,k,j)+vt(:,k,j))
               dens_orb_ORG(:)=dens_orb_ORG(:)+ AllZRIJ(i,j)*PSI(:,i)*Conjg(PSI(:,j))
               END DO
            END DO
         END DO

!      CI!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      Produce Cij

       DO k=1,numeig
          DO i=1,Morb
             DO j=1,Morb

                VIN_help=VIN; Vhelp=0.0d0; 
                IF(Morb.gt.1) THEN
                   CALL Produce_Cij(VIN_help,Vhelp,i,j,0,0,Nc)
                END IF
                CCu(j,i,k)=SUM(Cu(:,k)*CONJG(Vhelp(:)))
                CCv(i,j,k)=SUM(Cv(:,k)*Vhelp(:))

             END DO
          END DO
       END DO

!      Multiply with orbitals

       dens_CI_x=0.0d0  
       DO k=1,numeig
          DO i=1,Morb
             DO j=1,Morb

                dens_CI_x(:,k)=dens_CI_x(:,k)+  &
     &            CONJG(PSI(:,i))*PSI(:,j)*(CCu(i,j,k)+CCv(i,j,k))

             END DO
          END DO
       END DO

!      write all to (199,FILE='MC_dens_x.out')

       DO k=1,numeig
          DO i=1,ND

      WRITE(199,'(E16.8,3X,I8,3X,2E16.8,3X,2E16.8,3X,2E16.8,1X,F26.16,1X,F16.10)') ort_x(i), &
     & k,dens_orb_x(i,k),dens_CI_x(i,k),dens_orb_ORG(i),Real(wo(k)),weight(i)

          END DO
!STR 21.02.2013       write(6,*)Real(wo(k)),"DNS CI",SUM(ABS(dens_CI_x(:,k))),"Orb", SUM(ABS(dens_orb_x(:,k))),"ORG", ABS(SUM(dens_orb_ORG(:)))
       END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! WEIGHTS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!        ORBITAL!!!!!!!!!!!!!!!!!!!!!!

!        first get square roots of densities
         IF(Morb.gt.1) THEN
            CALL squarerootR(RhoSQR)
         ELSE 
            RhoSQR(1,1)=SQRT(REAL(AllZRIJ(1,1)))
!JL ORG            RhoSQR(1,1)=REAL(AllZRIJ(1,1))
         END IF

!        multiply square roots with amplitudes

         ut=0.0d0; vt=0.0d0

         DO k=1,numeig
            DO j=1,Morb
               DO i=1,Morb

                  ut(:,k,j)=ut(:,k,j)+RhoSQR(j,i)*u(:,k,i)
                  vt(:,k,j)=vt(:,k,j)+RhoSQR(j,i)*v(:,k,i)
                
               END DO
            END DO
         END DO

!        weights for different external driving: f=x,x^2
             tr=(0d0,0d0)
            Do i=1,ND
            tr(i,1)=ort_x(i)
            tr(i,2)=ort_x(i)**2
            ENDDO

         weight_orb=0.0d0  
         DO k=1,numeig
            DO i=1,Morb

               weight_orb(k,1)=weight_orb(k,1) +  &
     &   SUM(tr(:,1) *(PSI(:,i)*CONJG(ut(:,k,i))+CONJG(PSI(:,i))*CONJG(vt(:,k,i))))
!     &            SUM((PSI(:,i)*CONJG(ut(:,k,i))+CONJG(PSI(:,i))*CONJG(vt(:,k,i))))
!OLD     &            SUM(ort_x(:)*(PSI(:,i)*CONJG(ut(:,k,i))+CONJG(PSI(:,i))*CONJG(vt(:,k,i))))*dx
              weight_orb(k,2)=weight_orb(k,2) +  &
     &   SUM(tr(:,2)*(PSI(:,i)*CONJG(ut(:,k,i))+CONJG(PSI(:,i))*CONJG(vt(:,k,i))))
!OLD     &            SUM(ort_x(:)*ort_x(:)*(PSI(:,i)*CONJG(ut(:,k,i))+CONJG(PSI(:,i))*CONJG(vt(:,k,i))))*dx
            END DO

!          if(k==107) then
!          write(6,*)" Orb 107 ", k, sum(abs(imag(ut(:,k,1)))),sum(abs(imag(ut(:,k,2))))
!          qqq(1)=SUM(tr(:,1) *(PSI(:,1)*CONJG(ut(:,k,1))+CONJG(PSI(:,1))*CONJG(vt(:,k,1)))) 
!          qqq(2)=SUM(tr(:,1) *(PSI(:,2)*CONJG(ut(:,k,2))+CONJG(PSI(:,2))*CONJG(vt(:,k,2)))) 
!          write(6,'(a10,F10.5,a10,F10.6)')"u1|f|psi1=",abs(qqq(1))," u2|f|psi2=",abs(qqq(2))
!           pause
!          endif

         END DO

!      CI!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      Produce Cij

       DO k=1,numeig
          DO i=1,Morb
             DO j=1,Morb

                VIN_help=VIN; Vhelp=0.0d0; 
                IF(Morb.gt.1) THEN
                    CALL Produce_Cij(VIN_help,Vhelp,i,j,0,0,Nc)
                END IF
                CCu_r(i,j,k)=SUM(CONJG(Cu(:,k))*Vhelp(:))
                CCv_r(j,i,k)=SUM(CONJG(Cv(:,k)*Vhelp(:)))
!            if(k==265) then
!           write(6,*)" Prod_CI", i,j
!            do ll=1,Nc
!          write(6,'(i2,a4,F10.5,a6,F10.6,a6,F10.6)')ll,"Cu=",real(Cu(ll,k))," Vhelp=",real(Vhelp(ll)), "VIn=",real(VIN(ll))
!            enddo
!            endif
             END DO
          END DO
       END DO

!      Orbital integrals for different external driving: f=1,x,x^2

       DO i=1,Morb
          DO j=1,Morb

!OLD             orb_int(i,j,1)=SUM(CONJG(PSI(:,i))*PSI(:,j)*ort_x(:)*dx)
!OLD             orb_int(i,j,2)=SUM(CONJG(PSI(:,i))*PSI(:,j)*ort_x(:)*ort_x(:)*dx)
!             orb_int(i,j,1)=SUM(CONJG(PSI(:,i))*PSI(:,j))
             orb_int(i,j,1)=SUM(CONJG(PSI(:,i))*PSI(:,j)*tr(:,1))
             orb_int(i,j,2)=SUM(CONJG(PSI(:,i))*PSI(:,j)*tr(:,2))

          END DO
       END DO

!      Weights for different external driving: f=1,x,x^2

         weight_CI=0.0d0  
         DO k=1,numeig
            DO i=1,Morb
               DO j=1,Morb

                  weight_CI(k,1)=weight_CI(k,1) +  &
     &               orb_int(i,j,1)*(CCu_r(i,j,k)+CCv_r(i,j,k))
                  weight_CI(k,2)=weight_CI(k,2) +  &
     &               orb_int(i,j,2)*(CCu_r(i,j,k)+CCv_r(i,j,k))
!            IF(k==6) WRITE(6,*) i,j,orb_int(i,j,1),orb_int(i,j,2)
            IF(k==265) WRITE(196,*) ABS(wo(k)), i,j,CCu_r(i,j,k),CCv_r(i,j,k)
               END DO
            END DO
         END DO
     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! Density matrix of excitations (C_u only)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
!      Produce Cij

       DO k=1,numeig
          DO i=1,Morb
             DO j=1,Morb

                VIN_help=Cu(:,k); Vhelp=0.0d0; 
                IF(Morb.gt.1) THEN
                   CALL Produce_Cij(VIN_help,Vhelp,i,j,0,0,Nc)
                END IF
                CCu(j,i,k)=SUM(Cu(:,k)*CONJG(Vhelp(:)))

             END DO
          END DO
                CCu_norm(k)=SUM(Cu(:,k)*CONJG(Cu(:,k)))
       END DO

!   WRITE spectrum and weights
!         DO i=1,dimL/2
       WRITE(11,'(50a,50a,80a)') &
        "#      i,                 Re(w_i),              ",&
        "Re(w_i+Eref),         Re(Norm_Orb),     Re(Norm_CI)",&
        ",   Abs(<PSI_ref|OP1(x)|LR_PSI_i>),  Abs( <PSI_ref|OP2(x)|LR_PSI_i>)" 
               WRITE(12,'(50a,54a,80a)')& 
        "#      i,                 Re(w_i),              ",&
        "Re(w_i+Eref),         Re(Norm_Orb),     Re(Norm_CI),  ",&
        " Abs(<PSI_ref|OP1(x)|LR_PSI_i>),  Abs( <PSI_ref|OP2(x)|LR_PSI_i>)" 
         DO i=1,numeig

       WRITE(198,'(I8,E26.16,1X,F26.16,1X,F26.16,1X,2E16.8,2E16.8,2E16.8,2E16.8,2E16.8)') &
     & i, abs(wo(i))+Energy,abs(weight_orb(i,1)+weight_CI(i,1)),abs(weight_orb(i,2)+weight_CI(i,2)) &
     & ,weight_orb(i,1),weight_orb(i,2),weight_CI(i,1),weight_CI(i,2),wo(i)                                    ! 198=MC_we.out
!     & i, wo(i),abs(weight_orb(i,1)+weight_CI(i,1)),abs(weight_orb(i,2)+weight_CI(i,2)) &
!     & ,weight_orb(i,1),weight_orb(i,2),weight_CI(i,1),weight_CI(i,2)
!JL ORG            WRITE(198,'(I8,1X,2E16.8,1X,2E16.8,1X,2E16.8,1X,2E16.8,1X,2E16.8)') &
!JL ORG     & i, wo(i), weight_orb(i,1),weight_orb(i,2),weight_CI(i,1),weight_CI(i,2)
       WRITE(11,'(I8,2(F26.16,1X),6(F16.8,1X))') &
     & i,Real(wo(i)),Real(wo(i))+Energy,Real(norm_mat_orb(i)),Real(norm_mat_CI(i)),&
     & abs(weight_orb(i,1)), abs(weight_CI(i,1)),abs(weight_orb(i,2)), abs(weight_CI(i,2))  !11,FILE='MC_anlsplotALL.out'
        IF(real(norm_mat_orb(i)+norm_mat_CI(i)).ge.0.9d0) &
     & WRITE(12,'(I8,F26.16,1X,E26.16,1X,F20.16,1X,F20.16,1X,F20.16,1X,F20.16)') &
     & i,Real(wo(i)),Real(wo(i))+Energy,Real(norm_mat_orb(i)),Real(norm_mat_CI(i)),&
     & abs(weight_orb(i,1)+weight_CI(i,1)),abs(weight_orb(i,2)+weight_CI(i,2)) 
         END DO

!   WRITE Density matrices to (197,FILE='MC_rho.out')

       DO k=1,numeig
          DO i=1,Morb
             DO j=1,Morb

                WRITE(197,'(I8,1X,I8,1X,I8,1X,2E16.8,1X,2E16.8,1X,F26.16)') k,i,j,CCu(j,i,k),CCu_norm(k), Real(wo(k)) 

             END DO
          END DO
       END DO

!       CLOSE(196-199);
       CLOSE(196);
       CLOSE(197);
       CLOSE(198);
       CLOSE(199);
       CLOSE(11);
       CLOSE(12);

      END SUBROUTINE evaluate

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Check norm and orthogonalization of u and v
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE norm_full(PSI,vr,norm_mat,orth_mat,indorder2,u,v,Cu,Cv,Nc,norm_mat_orb,norm_mat_CI)

          USE CI_ALL
          USE Parallel_CI 
          USE DVR_ALL

         IMPLICIT NONE

         COMPLEX*16, DIMENSION(ND,Morb), INTENT(IN) :: PSI
         COMPLEX*16, DIMENSION(dimL,numeig), INTENT(IN) :: vr
         INTEGER, DIMENSION(numeig), INTENT(IN) :: indorder2
          INTEGER, INTENT(IN) :: Nc

         COMPLEX*16, DIMENSION(numeig,numeig), INTENT(OUT) :: norm_mat
         COMPLEX*16, DIMENSION(numeig,numeig), INTENT(OUT) :: orth_mat
         COMPLEX*16, DIMENSION(ND,numeig,Morb), INTENT(OUT) :: u, v
         COMPLEX*16, DIMENSION(Nc,numeig), INTENT(OUT) :: Cu, Cv
         COMPLEX*16, DIMENSION(numeig), INTENT(OUT) :: norm_mat_orb
         COMPLEX*16, DIMENSION(numeig), INTENT(OUT) :: norm_mat_CI

         COMPLEX*16, DIMENSION(dimL_orb,dimL_orb) :: Pb
         COMPLEX*16, DIMENSION(dimL_orb/2,dimL_orb/2) :: P
         COMPLEX*16, DIMENSION(ND) :: uhelp, vhelp

         INTEGER :: i,j,q,k,help
         COMPLEX*16, DIMENSION(numeig,Morb) :: ol

!     create u,v
         DO j=1,numeig
            DO q=1,Morb
               DO i=1,ND
                  u(i,j,q)=vr(i+(q-1)*ND,indorder2(j))
                  v(i,j,q)=vr(i+(Morb+q-1)*ND,indorder2(j))
               END DO
            END DO
         END DO

!!     project 
!
!         CALL get_proj(PSI,Pb,P)
!
!         DO j=1,numeig
!            DO q=1,Morb
!            
!               CALL ZGEMV('N', ND, ND, DCMPLX(1.0d0), P, &
!     &             ND, u(:,j,q), 1, DCMPLX(0.0d0), uhelp, 1)
!
!               u(:,j,q)=uhelp
!
!               CALL ZGEMV('N', ND, ND, DCMPLX(1.0d0), CONJG(P), &
!     &             ND, v(:,j,q), 1, DCMPLX(0.0d0), vhelp, 1)
!
!               v(:,j,q)=vhelp
!
!            END DO
!         END DO

!     create Cu,Cv
         DO j=1,numeig
             DO i=1,Nc
                  Cu(i,j)=vr(i+2*ND*Morb,indorder2(j))
                  Cv(i,j)=vr(i+Nc+2*ND*Morb,indorder2(j))
             END DO
         END DO

!     norm: u,v
         norm_mat=0.0d0; orth_mat=0.0d0

         DO q=1,Morb
            DO i=1,numeig
               DO j=1,numeig

                  norm_mat(i,j)= norm_mat(i,j)+ &
     &                SUM(CONJG(u(:,i,q))*u(:,j,q) &
     &               -CONJG(v(:,i,q))*v(:,j,q))

                  orth_mat(i,j)= orth_mat(i,j)+ &
     &                SUM(CONJG(v(:,i,q))*CONJG(u(:,j,q)) &
     &               -CONJG(u(:,i,q))*CONJG(v(:,j,q)))

               END DO
            END DO
         END DO

!         norm_mat=norm_mat*dx; orth_mat=orth_mat*dx
         norm_mat=norm_mat; orth_mat=orth_mat

!     norm: Cu,Cv
            DO i=1,numeig
               DO j=1,numeig

                  norm_mat(i,j)= norm_mat(i,j)+ &
     &                SUM(CONJG(Cu(:,i))*Cu(:,j) &
     &               -CONJG(Cv(:,i))*Cv(:,j))

                  orth_mat(i,j)= orth_mat(i,j)+ &
     &                SUM(CONJG(Cv(:,i))*CONJG(Cu(:,j)) &
     &               -CONJG(Cu(:,i))*CONJG(Cv(:,j)))

               END DO
            END DO
!        write(6,*)("NormOrth(8,",i,")=",norm_mat(8,i),i=1,numeig)
!        write(6,*)("Orth(",i,")=",orth_mat(i,i),i=1,numeig)
!        write(6,*)("Normmat(",i,")=",norm_mat(i,i),i=1,numeig)
            DO i=1,numeig
            IF(ABS(norm_mat(i,i))<1.0d-1) write(6,*)"o! Root N.",i," has ZERO norm",norm_mat(i,i)
            END DO


!      normalize states: u,v
         DO q=1,Morb
            DO i=1,numeig

               IF(ABS(norm_mat(i,i))>1.0d-8) THEN
                  u(:,i,q)=u(:,i,q)/SQRT(norm_mat(i,i))
                  v(:,i,q)=v(:,i,q)/SQRT(norm_mat(i,i))
               END IF

            END DO
         END DO


!      normalize states: Cu,Cv
         DO q=1,Nc
            DO i=1,numeig

               IF(ABS(norm_mat(i,i))>1.0d-8) THEN
                  Cu(q,i)=Cu(q,i)/SQRT(norm_mat(i,i))
                  Cv(q,i)=Cv(q,i)/SQRT(norm_mat(i,i))
               END IF

            END DO
         END DO

!     norm again: u, v
         norm_mat=0.0d0; orth_mat=0.0d0
         norm_mat_orb=0.0d0; norm_mat_CI=0.0d0

         DO q=1,Morb
            DO i=1,numeig
                  norm_mat_orb(i)= norm_mat_orb(i)+ &
     &                SUM(CONJG(u(:,i,q))*u(:,i,q) &
     &               -CONJG(v(:,i,q))*v(:,i,q))
               DO j=1,numeig

                  norm_mat(i,j)= norm_mat(i,j)+ &
     &                SUM(CONJG(u(:,i,q))*u(:,j,q) &
     &               -CONJG(v(:,i,q))*v(:,j,q))

                  orth_mat(i,j)= orth_mat(i,j)+ &
     &                SUM(CONJG(v(:,i,q))*CONJG(u(:,j,q)) &
     &               -CONJG(u(:,i,q))*CONJG(v(:,j,q)))

               END DO
            END DO
         END DO

         norm_mat=norm_mat; orth_mat=orth_mat
         norm_mat_orb=norm_mat_orb;

!     norm again: Cu,Cv
            DO i=1,numeig
                  norm_mat_CI(i)= norm_mat_CI(i)+ &
     &                SUM(CONJG(Cu(:,i))*Cu(:,i) &
     &               -CONJG(Cv(:,i))*Cv(:,i))
               DO j=1,numeig

                  norm_mat(i,j)= norm_mat(i,j)+ &
     &                SUM(CONJG(Cu(:,i))*Cu(:,j) &
     &               -CONJG(Cv(:,i))*Cv(:,j))

                  orth_mat(i,j)= orth_mat(i,j)+ &
     &                SUM(CONJG(Cv(:,i))*CONJG(Cu(:,j)) &
     &               -CONJG(Cu(:,i))*CONJG(Cv(:,j)))

               END DO
            END DO

!        write(6,*)("CIOrb(",i,")=",norm_mat_CI(i)+norm_mat_orb(i),i=1,numeig)
!        write(6,*)("UTorth(1,",i,")=",orth_mat(1,i),i=1,numeig)


      END SUBROUTINE norm_full

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Check norm and orthogonalization of u and v: orbital part
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE norm_orb(PSI,vr,norm_mat,orth_mat,indorder2,u,v)

          USE CI_ALL
          USE Parallel_CI 
          USE DVR_ALL

         IMPLICIT NONE

         COMPLEX*16, DIMENSION(ND,Morb), INTENT(IN) :: PSI
         COMPLEX*16, DIMENSION(dimL_orb,numeig), INTENT(IN) :: vr
         INTEGER, DIMENSION(numeig), INTENT(IN) :: indorder2

         COMPLEX*16, DIMENSION(numeig,numeig), INTENT(OUT) :: norm_mat
         COMPLEX*16, DIMENSION(numeig,numeig), INTENT(OUT) :: orth_mat
         COMPLEX*16, DIMENSION(ND,numeig,Morb), INTENT(OUT) :: u, v

         COMPLEX*16, DIMENSION(dimL_orb,dimL_orb) :: Pb
         COMPLEX*16, DIMENSION(dimL_orb/2,dimL_orb/2) :: P
         COMPLEX*16, DIMENSION(ND) :: uhelp, vhelp

         INTEGER :: i,j,q,k,help
         COMPLEX*16, DIMENSION(numeig,Morb) :: ol

!     create u,v
         DO j=1,numeig
            DO q=1,Morb
               DO i=1,ND
                  u(i,j,q)=vr(i+(q-1)*ND,indorder2(j))
                  v(i,j,q)=vr(i+(Morb+q-1)*ND,indorder2(j))
               END DO
            END DO
         END DO

!!     project 
!
!         CALL get_proj(PSI,Pb,P)
!
!         DO j=1,numeig
!            DO q=1,Morb
!            
!               CALL ZGEMV('N', ND, ND, DCMPLX(1.0d0), P, &
!     &             ND, u(:,j,q), 1, DCMPLX(0.0d0), uhelp, 1)
!
!               u(:,j,q)=uhelp
!
!               CALL ZGEMV('N', ND, ND, DCMPLX(1.0d0), CONJG(P), &
!     &             ND, v(:,j,q), 1, DCMPLX(0.0d0), vhelp, 1)
!
!               v(:,j,q)=vhelp
!
!            END DO
!         END DO

!     norm: u,v
         norm_mat=0.0d0; orth_mat=0.0d0

         DO q=1,Morb
            DO i=1,numeig
               DO j=1,numeig

                  norm_mat(i,j)= norm_mat(i,j)+ &
     &                SUM(CONJG(u(:,i,q))*u(:,j,q) &
     &               -CONJG(v(:,i,q))*v(:,j,q))

                  orth_mat(i,j)= orth_mat(i,j)+ &
     &                SUM(CONJG(v(:,i,q))*CONJG(u(:,j,q)) &
     &               -CONJG(u(:,i,q))*CONJG(v(:,j,q)))

               END DO
            END DO
         END DO

!OLD         norm_mat=norm_mat*dx; orth_mat=orth_mat*dx
         norm_mat=norm_mat; orth_mat=orth_mat

!      normalize states: u,v
         DO q=1,Morb
            DO i=1,numeig

!               IF(ABS(norm_mat(i,i))>1.0d-1) THEN
!                  u(:,i,q)=u(:,i,q)/SQRT(norm_mat(i,i))
                  v(:,i,q)=v(:,i,q)/SQRT(norm_mat(i,i))
!               END IF

            END DO
         END DO


!     norm again: u, v
         norm_mat=0.0d0; orth_mat=0.0d0

         DO q=1,Morb
            DO i=1,numeig
               DO j=1,numeig

                  norm_mat(i,j)= norm_mat(i,j)+ &
     &                SUM(CONJG(u(:,i,q))*u(:,j,q) &
     &               -CONJG(v(:,i,q))*v(:,j,q))

                  orth_mat(i,j)= orth_mat(i,j)+ &
     &                SUM(CONJG(v(:,i,q))*CONJG(u(:,j,q)) &
     &               -CONJG(u(:,i,q))*CONJG(v(:,j,q)))

               END DO
            END DO
         END DO

!OLD         norm_mat=norm_mat*dx; orth_mat=orth_mat*dx
         norm_mat=norm_mat; orth_mat=orth_mat

      END SUBROUTINE norm_orb

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     analyse symmetries of the matrix
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE analyse_matrix(L,L1,L3,Nc)

          USE CI_ALL
          USE Parallel_CI 
          USE DVR_ALL
!         USE mkl95_lapack
!         USE mkl95_precision

         IMPLICIT NONE

         COMPLEX*16, DIMENSION(dimL,dimL), INTENT(IN) :: L
          INTEGER, INTENT(IN) :: Nc
         COMPLEX*16, DIMENSION(dimL,dimL), INTENT(OUT) :: L1, L3

         COMPLEX*16, DIMENSION(dimL,dimL) :: sigma1, sigma3, C1, C2
         INTEGER :: i,j

             write(6,*)" Nc in analyze ..",dimL,Nc,dimL_orb
                  pause

         sigma1=0.0d0; sigma3=0.0d0; L1=0.0d0; L3=0.0d0

         DO i=1,dimL_orb/2

            sigma1(i,dimL_orb/2+i)=1.0d0
            sigma1(dimL_orb/2+i,i)=1.0d0

            sigma3(i,i)=1.0d0
            sigma3(dimL_orb/2+i,dimL_orb/2+i)=-1.0d0

         END DO

         DO i=1,Nc

            sigma1(dimL_orb+i,dimL_orb+Nc+i)=1.0d0
            sigma1(dimL_orb+Nc+i,dimL_orb+i)=1.0d0

            sigma3(dimL_orb+i,dimL_orb+i)=1.0d0
            sigma3(dimL_orb+Nc+i,dimL_orb+Nc+i)=-1.0d0

         END DO

         CALL ZGEMM('N', 'N', dimL, dimL, dimL, DCMPLX(1.0d0),sigma1, &
     &                      dimL, L, dimL, DCMPLX(0.0d0), C1, dimL)

         CALL ZGEMM('N', 'N', dimL, dimL, dimL, DCMPLX(1.0d0),C1, &
     &                      dimL, sigma1, dimL, DCMPLX(0.0d0), C2, dimL)

         L1=C2
         CALL ZGEMM('N', 'N', dimL, dimL, dimL, DCMPLX(1.0d0),sigma3, &
     &                      dimL, L, dimL, DCMPLX(0.0d0), C1, dimL)

         CALL ZGEMM('N', 'N', dimL, dimL, dimL, DCMPLX(1.0d0),C1, &
     &                      dimL, sigma3, dimL, DCMPLX(0.0d0), C2, dimL)

         L3=C2

!STR        DO i=1,dimL
!STR            DO j=1,dimL
!STR               WRITE(30,'(I8,1X,I8,1X,2E16.8)') i, j, &
!STR     &                          L(i,j)!L1(i,j)+CONJG(L(i,j))
!STR               WRITE(31,'(I8,1X,I8,1X,2E16.8)') i, j, &
!STR     &                          L3(i,j)-CONJG(L(j,i))
!STR               WRITE(32,'(I8,1X,I8,1X,2E16.8)') i, j, &
!STR     &                          sigma1(i,j)
!STR            END DO
!STR         END DO
!STR         CLOSE(30); CLOSE(31); CLOSE(32)

      END SUBROUTINE analyse_matrix

      SUBROUTINE Get_Exp_Val_1bOp(PSI,Op,rhoij,exp_Val)

          USE CI_ALL
          USE Parallel_CI 
          USE W_INTERPARTICLE !ADD BY STR
          USE DVR_ALL
          USE rR_hW

          IMPLICIT NONE         

          INTEGER :: i, j, k

          COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb), INTENT(IN) :: PSI
          COMPLEX*16, DIMENSION(NDX*NDY*NDZ,NDX*NDY*NDZ), INTENT(IN) :: Op
          COMPLEX*16, DIMENSION(Morb,Morb) :: rhoij
          REAL*8, INTENT(OUT) :: exp_Val
!           write(6,*)"<fx>ij=",Matmul(Transpose(Conjg(PSI)),Matmul(Op,PSI)) 
          exp_val=Abs(SUM(Matmul(Transpose(Conjg(PSI)),Matmul(Op,PSI))*rhoij))

      END SUBROUTINE Get_Exp_Val_1bOp

END MODULE LINEAR_RESPONSE

