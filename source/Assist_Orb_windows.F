C=============== Assistent MPI subroutive -assists in RHS|psi>  evaluation 
       subroutine Assist_Orb
       USE SHARED_DIMS
       USE DVR_ALL
       USE Parallel_Orb
       USE rR_hW
       USE W_INTERPARTICLE
       implicit NONE
c=================== MPI 
       INCLUDE 'mpif.h'
       INTEGER ::  ierr,MYID,numprocs
c===================================================================
         REAL*8 :: local_time
c===================================================================
       INTEGER ::  n,Iorb,IPRC,m,nn,icntr,FromN,TillN,i
c       COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb) :: PSI
c       COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb) :: PSI1,PSI2,temp
c       COMPLEX*16, DIMENSION(NDX*NDY*NDZ) ::tempX,tempY,tempZ

c       COMPLEX*16,  DIMENSION(NDX*NDY*NDZ,Morb) ::  PSI, OPsi
c       COMPLEX*16, ALLOCATABLE ::  PSI(:,:), OPsi(:,:)
       COMPLEX*16, ALLOCATABLE ::  PSI(:,:),temp(:,:) !,temp1(:,:)
c       COMPLEX*16, ALLOCATABLE :: tempX(:), tempY(:), tempZ(:)

       COMPLEX*16 :: RhoAll, ZDOTC,Z
       EXTERNAL Get_Op_PSI,ZSCAL,VTRAP_EXT_TD
c===================================================================
         allocate(PSI(NDX*NDY*NDZ,Morb),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in FUNCSTR Asistent"
         allocate(temp(NDX*NDY*NDZ,Morb),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in FUNCSTR Asistent"
!         allocate(temp1(NDX*NDY*NDZ,Morb),stat=ierr)
!         if(ierr /= 0)write(*,*)"allocation error in FUNCSTR Asistent"
c========================= MPI =====================================
         call MPI_COMM_RANK(MPI_COMM_WORLD,MYID,ierr)
         call MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,ierr)
         CALL Share_r_R
c===================================================================
               n=NDX*NDY*NDZ*MORB
               m=NDX*NDY*NDZ
c===================================================================
1       continue 
      call MPI_BCAST(JOB_TYPE,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c====================  JOB_TYPE = 1 !  RHS ordinary   ============== 
c====================  JOB_TYPE = 2 !  RHS ENERGY EVALUATION ======= 
c====================  JOB_TYPE = 999 !  Terminates Jobs =========== 
c===================================================================
          IF(JOB_TYPE.eq.1) ENRG_EVAL=.FALSE.
          IF(JOB_TYPE.eq.2) ENRG_EVAL=.TRUE.
c===================================================================
          JOBTYPE: SELECT CASE (JOB_TYPE)
                   CASE (1,2)
        IPRC=MYID+1
c        CALL MPI_ScatterV(PSI,Orb_Block,Orb_Displ,MPI_DOUBLE_COMPLEX,
c     .  temp,Orb_Block(IPRC),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(PSI,n,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
c===================== GET EXTERNAL TRAP ===========================
      call MPI_BCAST(local_time,1,MPI_DOUBLE_PRECISION,0,
     .               MPI_COMM_WORLD,ierr)
c================ WSL Part =========================================
        FromN=Proc_WSL_Starts(IPRC)
        TillN=Proc_WSL_Finish(IPRC)
          IF(TillN.ne.-1) THEN ! Case Wsl - PMI paralellized
c        CALL Get_OPSI_WSL_delta_OMP (PSI,temp,FromN,TillN)
c        CALL Get_OPSI_WSL_based_OMP (PSI,temp,FromN,TillN)
       IF(DIM_MCTDHB.eq.1.and.Wxx_TYPE.eq.0.and.NDX.le.1024) THEN
       CALL Get_OPSI_WSL_disbalanced_OMP (PSI,temp,FromN,TillN,
     &                                    local_time) !AIS 15JUL2014
       ELSE
       CALL Get_OPSI_WSL_balanced_OMP (PSI,temp,FromN,TillN,local_time) !AIS 15JUL2014
       ENDIF
c=====================================================================
!               temp1=temp
      CALL MPI_REDUCE(temp,temp,n,MPI_DOUBLE_COMPLEX,MPI_SUM,0,
     .                                          MPI_COMM_WORLD,ierr)
        IF(ENRG_EVAL.eqv..TRUE.) THEN
!                 WIJKL1=WIJKL
      CALL MPI_REDUCE(WIJKL,WIJKL,Rdim1
     . ,MPI_DOUBLE_COMPLEX,MPI_SUM,0, MPI_COMM_WORLD,ierr)
!      CALL MPI_REDUCE(WIJKL1,WIJKL,Rdim1,
!     . MPI_DOUBLE_COMPLEX,MPI_SUM,0, MPI_COMM_WORLD,ierr)
        ENDIF
          ENDIF
c==================================================================
          temp=PSI
       CALL OPERS(temp,local_time)
c================ Gathering togeather MPI_REDUCE ===================
       CALL  MPI_GatherV(temp,Orb_Block(IPRC),MPI_DOUBLE_COMPLEX,PSI,
     . Orb_Block,Orb_Displ,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)

                   goto 1 
                   CASE (999)
c       write(6,*)"Orbital  Assistent Finishes 
c     .with P_R_O_P_E_R  CASE 999!!!"
c       write(6,*)"Orbital  Assistent ",MYID," is DONE"
                   GOTO 9999
                   END SELECT JOBTYPE

c================= END MPI
9999   continue
        deallocate(PSI,temp)
!        deallocate(temp1)
c       write(6,*)"Orbital  Assistent ",MYID," is DONE"
       return
       end subroutine Assist_Orb
