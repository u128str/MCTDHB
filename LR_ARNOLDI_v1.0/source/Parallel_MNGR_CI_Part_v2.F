! Subroutine to provides jobs to different process 
! Operates in the assumptio that number of processes (processors) in not
! large than Morb+Morb*(Morb+1)/2;
! in principle there sould be several possibilities:
!                          a) delta-function - only linear part is
!                          needed to be parallelized;
!                          b)  True separable potential 

!V1 - in redistribution of the JOBS it is now taken onto account that Master is doing whole 1B part plus few terms from 2B part
! in the limiting case when NPROC is large Master is doing ALL 1B terms plus ONE (the first) two-body term!!!!!
!
!V2 - parallilization of the 1body block - one hase to introduce work_1b and work_2b integer arrays
! Allocation of arrays for CI_Prod module - fills ones these arrays and
! then reuse them each time when needed for CI part.
          SUBROUTINE  MNGR_CI_Part(NprocORG)
          USE SHARED_DIMS
          USE Parallel_CI
          USE rR_hW
          USE CI_prod
          USE CI_ALL
          IMPLICIT NONE
          INCLUDE 'mpif.h'
c==============================================================================
          INTEGER :: MYID
          INTEGER :: Nterms,FromN,TillN,Master_From, Master_Till,MaxTrm
          INTEGER :: I,J,K,icntr,Nproc,M,NprocORG,ia,ib
          INTEGER ::  ierr,IPRC,TRM
c==============================================================================
c          INTEGER, ALLOCATABLE  :: Ind_CI_1b(:,:)
c          INTEGER, ALLOCATABLE  :: Ind_CI_2b(:,:)
c          Real*8,  ALLOCATABLE  :: Prefactors_1b(:,:)
c          Real*8,  ALLOCATABLE  :: Prefactors_2b(:,:)
c==============================================================================
!       I USE "Parallel" arrays: Proc_Job(Nproc) 
!                                Proc_Iorb_Starts -> summation over Iorb in Func_Str 
!                                Proc_Iorb_Finish
!                                Proc_Jdim_Starts -> summation over Dimensionalities in Func_St
!                                Proc_Jdim_Finish
!                                Proc_From
!                                Proc_Till -> part which is computed in   Get_Op_PSI(Vec,WOp_X,NDX,1,From,Till)
c==============================================================================
c============================== If (.FALSE.) W(x-x')=delta(x-x')*xlambda0 is used
c============================== or if (.TRUE.) - true interparticle potential 
c============================== W(r-r')*xlambda0 is used
c============================== defined by Wxx(I,J)=W(r,r') in DVR grid
c============================== or by vector W2xx(J)=W(ort(1)-ort(J))
c============================== or for separable case: WOP_X*WOP_Y*WOP_Z
       WRITE(6,*)"CI Parallel Menager CI Part ",NprocORG," processors"
       WRITE(6,*)"CI MaxTrm1b",MaxTrm1b," MaxTrm2b",MaxTrm2b

        CI_Proc_From=0
        CI_Proc_Till=0


        IF(NprocORG.eq.1) THEN
        Nproc=NprocORG
        CI_Proc_From(1)=1
        CI_Proc_Till(1)=MaxTrm2b
        write(6,*)"CI Master does all 2b terms from",CI_Proc_From(1),
     . " till",CI_Proc_Till(1),"!"
c        return
        ENDIF

        Nproc=NprocORG

        MaxTrm=MaxTrm2b+MaxTrm1b
        Master_From=1
        Master_Till=1
         
       WRITE(6,*)"CI MaxTrm",MaxTrm

!       IF(NPROC.GT.MaxTrm) stop "NPROC>MaxTrm: reduce nproc!!!"

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
      write(6,*)"CI oRG Proc ",MYID+1," works from",FromN,
     ." till",TillN,"!"
       ENDDO


c       EXCEPTIONAL_ID=-1 ! negative means all jobs are nicely redistributed between PROCS
        MYID_TRM=0
       DO I=1,Nproc
      IF(CI_Proc_From(I).le.MaxTrm1b.and.CI_Proc_Till(I).gt.MaxTrm1b)
     . EXCEPTIONAL_ID=I-1
      IF(CI_Proc_From(I).le.MaxTrm1b.and.CI_Proc_Till(I).gt.MaxTrm1b)
     . MYID_TRM(I)=1

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
C Proc            1  works from           1  till          42 ! tRM:           3
C Proc            2  works from          43  till         105 ! tRM:           2
C Proc            3  works from         106  till         168 ! tRM:           2
C Proc            4  works from         169  till         231 ! tRM:           2


c          CI_Proc_Till(1)=
c          CI_Proc_Till(1)=18
c          CI_Proc_From(2)=19
c          CI_Proc_Till(2)=80
c          CI_Proc_From(3)=81
c          CI_Proc_Till(3)=146
c          CI_Proc_From(4)=147
c          CI_Proc_Till(1)=CI_Proc_Till(1)-10
c          CI_Proc_From(1)=CI_Proc_From(1)-10
c          CI_Proc_Till(1)=CI_Proc_Till(1)-10
c          CI_Proc_From(2)=CI_Proc_From(2)-10
c

      write(6,*)"CI Exceptional Proc:",EXCEPTIONAL_ID+1,"works 1b [",
     .CI_Proc_From(EXCEPTIONAL_ID+1),":",MaxTrm1b,"] 2b [1",
     .":",CI_Proc_Till(EXCEPTIONAL_ID+1),"]" 
       DO I=1,NprocORG
       write(6,*)"CI Proc",I," works from",CI_Proc_From(I),
     ." till",CI_Proc_Till(I),"!"," tRM:",MYID_TRM(I)
       ENDDO
c        stop
c==================== Allocation of memory for Prefactors_1b,Prefactors_2b,Ind_CI_1b,Ind_CI_2b arrays =====================
c==================== These arrays, different on each node will be used after CI_Production EQ TRUE =======================
            IF(Morb.eq.1) CI_PRD=.TRUE. !For GP there is no need for parallelization
            CI_Prd=.TRUE.                                 

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
C         nYdim= ! number of the ci vectors needed work on each node - depends on MYID


           TRM=MYID_TRM(MYID+1)
           TRM_choice: SELECT CASE (TRM)
                    CASE (1) ! Worker on 1b
            nYdim=TillN-FromN+1
        ALLOCATE( Prefactors_1b(nXdim,nYdim) )
        IF(ALLOCATED(Prefactors_1b).eqv..FALSE.)
     . write(6,*)" MEM for Prefactors_1b is NOT ok!",SIZE(Prefactors_1b)
        ALLOCATE( Ind_CI_1b(nXdim,nYdim) )
        IF(ALLOCATED(Ind_CI_1b).eqv..FALSE.)
     . write(6,*)" MEM for Ind_CI_1b is NOT ok!",SIZE(Ind_CI_1b)
        IF(ALLOCATED(Prefactors_1b).eqv..FALSE.) stop "Prefactors_1b"
        IF(ALLOCATED(Ind_CI_1b).eqv..FALSE.)     stop "Ind_CI_1b"
            Ind_CI_1b=1
           Prefactors_1b=0d0
                    CASE (2) ! Worker on 2b
            nYdim=TillN-FromN+1
        ALLOCATE( Prefactors_2b(nXdim,nYdim) )
        IF(ALLOCATED(Prefactors_2b).eqv..FALSE.)
     . write(6,*)" MEM for Prefactors_2b is NOT ok!",SIZE(Prefactors_2b)
        ALLOCATE( Ind_CI_2b(nXdim,nYdim) )
        IF(ALLOCATED(Ind_CI_2b).eqv..FALSE.)
     . write(6,*)" MEM for Ind_CI_2b is NOT ok!",SIZE(Ind_CI_2b)
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
        IF(ALLOCATED(Prefactors_1b).eqv..FALSE.)
     . write(6,*)" MEM for Prefactors_1b is NOT ok!",SIZE(Prefactors_1b)
        ALLOCATE( Ind_CI_1b(nXdim,nYdim) )
        IF(ALLOCATED(Ind_CI_1b).eqv..FALSE.)
     . write(6,*)" MEM for Ind_CI_1b is NOT ok!",SIZE(Ind_CI_1b)
        IF(ALLOCATED(Prefactors_1b).eqv..FALSE.) stop "Prefactors_1b"
        IF(ALLOCATED(Ind_CI_1b).eqv..FALSE.)     stop "Ind_CI_1b"
            Ind_CI_1b=1
           Prefactors_1b=0d0
c       write(6,*)"MYID:",MYID,"Ind_CI_1b",nYdim,"fr:",FromN," tl:",TillN
c            pause
        FromN=1
        TillN=CI_Proc_Till(EXCEPTIONAL_ID+1)
            nYdim=TillN-FromN+1
        ALLOCATE( Prefactors_2b(nXdim,nYdim) )
        IF(ALLOCATED(Prefactors_2b).eqv..FALSE.)
     . write(6,*)" MEM for Prefactors_2b is NOT ok!",SIZE(Prefactors_2b)
        ALLOCATE( Ind_CI_2b(nXdim,nYdim) )
        IF(ALLOCATED(Ind_CI_2b).eqv..FALSE.)
     . write(6,*)" MEM for Ind_CI_2b is NOT ok!",SIZE(Ind_CI_2b)
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

       END  SUBROUTINE MNGR_CI_Part
