C=============== Short Iterative Lanzcos Vin(o) -> Vout(time)
        subroutine SIL_PROPG(time,VIN,VOUT,Error_SIL,MAXIT)
!        USE PASS_ARG
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
c=================== MPI ==================================
       INCLUDE 'mpif.h'
       INTEGER ::  ierr,MYID,numprocs
c==========================================================
       integer :: i,j,k,l,P,MAXIT
       INTEGER ::  n,FromN,TillN
c====================== For SIL ===========================
        integer ::  maxsil,NL,Iter_Total,maxsil1,minsil
        COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: VSIL(:,:)
        REAL*8   , DIMENSION(125) :: SILALPHA,SILBETA,SILPROJ
        REAL*8 :: time,Error_SIL, E_state,xerr,t,t1,Error_SIL_WORK
        COMPLEX*16 :: Z,Z1
        COMPLEX*16, DIMENSION(125) :: SCVECP,SCVECM
c==============================================

      real*4 start,finish,exec_time ,finish1,iter_time       
      real*8   time_start,time_finish
      real*8  xnrm

      character*20 message
      character*16 lname 
      character*26 lnamenew,lnamerest  
      LOGICAL  CNV,SIL

      real*8 ::   FKTR,CNK,DZASUM,DZNRM2,DNRM2,DDOT
      COMPLEX*16 :: ZDOTC

      external FKTR,CNK
      external DZASUM,DZNRM2,DNRM2,ZDOTC,DDOT,ZSCAL
      external Get_Eigen_SIL

c===========================================================
         call MPI_COMM_RANK(MPI_COMM_WORLD,MYID,ierr)
         call MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,ierr)
c         write(6,*)"Master  Process ",MYID," of ", numprocs
c===========================================================
         Error_SIL_WORK=Error_SIL/1.0d0
       xerr=1.0d0
       t=time
c============================================
       maxsil=40
       maxsil=Time_maxsil
       minsil=4
       minsil=Time_minsil
c       maxsil=40
c       minsil=40
       MaxSil1=MaxSil
       ALLOCATE(VSIL(Nconf,maxsil))
       IF(ALLOCATED(VSIL).eqv..FALSE.) 
     . write(6,*) " MEM for SIL", maxsil," is NOT ok!", SIZE(VSIL)
c==========================================================================
c        write(6,*)"H  =",(ZHJK(TRM_REQ1B(i)),i=1,MaxTrm1b)
c        write(6,*)"W  =",(ZWIJKL(TRM_REQ2B(i)),i=1,MaxTrm2b)
c         pause
c         write(6,*)"SIL Vin 0",VIN
c===================== S I L  =============================================
              SIL=.TRUE.
              Iter_Total=0
c==========================================================================
99           continue 
              VSIL=ZERO
c              Error_SIL=1.0d-8
c              Error_SIL=1.0d-5
              CNV=.FALSE.
c================== On-line SIL integrator/diagonalizer
              VIN=VIN/DZNRM2(Nconf,VIN,1)
c==================== First Iteration out of DO-LOOP
c          call cpu_time(start)
             time_start=MPI_WTIME(ierr)

             VSIL(:,1)=VIN

         CALL HPSI(VIN,VOUT)

            VSIL(:,2)=VOUT

       IF(Morb.eq.1) THEN
       SIL=.FALSE.
       call MPI_BCAST(SIL,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
       VIN(1) =Zoner
       VOUT(1)=Zoner
c         write(6,*)"SIL Vin 2",VIN
       goto 102 
       ENDIF

          SILALPHA(1)=DREAL(ZDOTC(Nconf,VSIL(:,1),1,VSIL(:,2),1))
            VSIL(:,2)=VSIL(:,2)-VSIL(:,1)*SILALPHA(1)
           SILBETA(1)=DZNRM2(Nconf,VSIL(:,2),1)
            VSIL(:,2)=VSIL(:,2)/SILBETA(1)
              E_state=SILALPHA(1)
            Iter_Total=Iter_Total+1


        xerr=SILBETA(1)*t

           Do J=2,maxsil-1
c==========================  VSIL(:,I-1)=|PSI>
                  VIN=VSIL(:,J)
                  VOUT=ZERO
             CALL HPSI(VIN,VOUT)
                VSIL(:,J+1)=VOUT
c==========================  VSIL(:,I)=H|PSI> r=AV_j
               VSIL(:,J+1)=VSIL(:,J+1)-VSIL(:,J-1)*SILBETA(J-1) 
c==========================  r=r-V_{J-1} Beta_{j-1}
             SILALPHA(J)=DREAL(ZDOTC(Nconf,VSIL(:,J),1,VSIL(:,J+1),1))
               VSIL(:,J+1)=VSIL(:,J+1)-VSIL(:,J)*SILALPHA(J) 
c==========================  r=r-V_{J} Alpha_{j}
                  SILBETA(J)=DZNRM2(Nconf,VSIL(:,J+1),1)
                  VSIL(:,J+1)=VSIL(:,J+1)/SILBETA(J)
c================ Error evaluation =================================
              xerr=xerr*SILBETA(J)*t/J
              t1=t*(Error_SIL_WORK/xerr)**(1.0d0/J)

       Iter_Total=Iter_Total+1
       NL=J
       SILPROJ=SILALPHA

!     . "SIL subspace covers in [%]",
       IF(ABS(t1).ge.ABS(t).and.J.ge.minsil) THEN
       CNV=.TRUE.
         IF(print_level.ge.3) then
       write(6,'(I3,a20,F30.6,a49,F15.10)')J,
     . "SIL vectors cover ", 100*(Error_SIL_WORK/xerr)**(1.0d0/J),
     . "[%] of Krylov subspace needed for tau/2=",time
            ENDIF !PRINT_LEVEL
       goto 101
       ENDIF
       EndDO
c======================= END SIL DO -LOOP till maxsil iterations
c================ CASE  when MaxSil is achieved
       CNV=.TRUE.
       write(6,'(a45,I16,f20.16,f20.16)')
     $"NOT CONVERGED -- MAXSIL is achived !",NL,t,t1
       if (abs(t1).lt.abs(t)) then 
          time=t1
          write(6,'(a34,f20.16,a23,f20.16)')
     $"instead desired time interval tau=",t," SIL uses smaller tau=",t1
          t=t1
       end if
c================== CASE  when not enough SIL iterations 
101         CONTINUE
           if(100*(Error_SIL_WORK/xerr)**(1.0d0/J).lt.100) then
         IF(print_level.ge.3) then
       write(6,'(a17,I5,a20,F36.16,a20,F15.10)')
     . "SIL subspace NL=",J," covers in [%]",
     . 100*(Error_SIL_WORK/xerr)**(1.0d0/J)," of desired tau=",time
         ENDIF !PRINT_LEVEL
           endif

       SIL=.FALSE.
       call MPI_BCAST(SIL,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
c================= Solve Eigenproblem for current Krylov subspace DIM==J
c           write(6,*)"SIL BETA",(SILBETA(i),i=1,NL-1)
         CALL
     .  Get_Eigen_SIL(NL,SILPROJ,SILBETA,E_state,Error_SIL_WORK,
     .  CNV,maxsil,t,SCVECP,SCVECM)

c       call MPI_BCAST(SIL,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
           VIN=ZERO
           VOUT=ZERO
           DO J=1,NL
c          VIN=VIN+SCVECP(J)*VSIL(:,J)
c           VOUT=VOUT+SCVECM(J)*VSIL(:,J)
           VIN=VIN+SCVECM(J)*VSIL(:,J)
           VOUT=VOUT+SCVECP(J)*VSIL(:,J)
           ENDDO

c         VIN=ZERO
c         Do I=1,Nconf
c         Do K=1,NL
c         VIN(I)=VIN(I)+VSIL(I,K)*SILPROJ(K)
c         EndDo
c         EndDo

c           VOUT=VIN*EXP( PreFac*time) ! Clean and normally integrated Vout(time)
c           VIN =VIN*EXP(-PreFac*time) ! Clean and normally integrated Vout(-time) needed for error evaluation
102          continue
c       write(6,*)"Relax norm ",Real(Prefac),epsilon(0d0)
              IF (Abs(Real(Prefac)).ge.epsilon(0d0)) THEN
              xnrm=DZNRM2(Nconf,VIN,1)
              VIN=VIN/DZNRM2(Nconf,VIN,1)
              xnrm=DZNRM2(Nconf,VOUT,1)
              VOUT=VOUT/DZNRM2(Nconf,VOUT,1)
         IF(print_level.ge.3) then
        write(6,*)"Sil Relax CI vectors are renormalized",xnrm
           ENDIF !PRINT_LEVEL
              endif
c        write(6,*)Iter_Total,"In sil CMPL NORM :",xnrm
c        write(6,*)" Sil vec",SILPROJ
c             pause


c==============================================
        DEALLOCATE(VSIL)
c        VOUT=VIN
c         write(6,*)"SIL OUT"
        end subroutine SIL_PROPG
