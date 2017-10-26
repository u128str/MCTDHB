C================= Diagonalization of SIL subspace Tridiagonal matrix is C================= assumed to be constructed out off DIAG0,DIAG1
C       call  Get_Eigen_SIL(N,SILALPHA,SILBETA,E_state,Error_SIL,CNV)
       subroutine Get_Eigen_SIL(NL,D0,D1,E_state,Error_SIL,CNV,MaxSil,t,
     .                       SCVECP,SCVECM)
        USE SHARED_DIMS
        implicit NONE
        REAL*8 :: Error_SIL, E_state,t,ABSTOL,VL,VU
        REAL*8   , DIMENSION(125) :: D0,D1 
        COMPLEX*16, DIMENSION(125) :: SCVECP,SCVECM
        LOGICAL CNV
        integer NL,INFO,MaxSIl,I,J,K,IL,IU,M
        REAL*8   , DIMENSION(NL)     :: EVALS
        REAL*8   , DIMENSION(NL-1)   :: OffDg2
c        REAL*8   , DIMENSION(2*NL-2) :: WORK     
        REAL*8   , DIMENSION(5*NL) :: WORK     
        INTEGER   , DIMENSION(5*NL) :: IWORK     
        INTEGER   , DIMENSION(NL) :: IFAIL     
        REAL*8   , DIMENSION(NL,NL)  :: EV

c        COMPLEX*16, PARAMETER :: ZERO =(0d0,0d0)
        external DSTEQR
c========================== GP case 
        IF(Morb.eq.1) THEN
           SCVECP=ZERO
           SCVECM=ZERO
           SCVECP(1)=ZONER
           SCVECM(1)=ZONER
         return 
        ENDIF
c========================== GP case 

         EV=0.0d0
         Do i=1,NL
         EVALS(I)=D0(I)
         EV(I,I)=1.0d0
         Enddo
         Do i=1,NL-1
         OffDg2(I)=D1(I)
         Enddo
c         EV=0.0d0
c          write(6,*)"DIAG",EVALS
c          write(6,*)"OFFDIAG",OffDg2
c         call DSTEQR('I',NL,EVALS,OffDg2,EV,NL,WORK,INFO)
c        SUBROUTINE DSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )

c         call DSTEVX( JOBZ, RANGE,N, D, E, VL, VU, IL, IU, ABSTOL,
c     $                   M, W, Z, LDZ, WORK, IWORK, IFAIL, INFO )
           ABSTOL=1.0d-14
         call DSTEVX('V','A',NL, D0,OffDg2,VL,VU, IL, IU, ABSTOL,
     $                   M, EVALS, EV, NL, WORK, IWORK, IFAIL, INFO )

       IF(M.ne.NL) write(6,*)"Problems in SIL diag ierr",INFO

       IF(INFO.ne.0) write(6,*)"Problems in SIL diag ierr",INFO
c          write(6,*)"EVALS",EVALS
c          Nroot=NL
c          if(NL.GE.20) Nroot=20
c          Do K=1,Nroot
c          E_dinam(K)=EVALS(K)
c          ENDDO


c===============================================================
            IF (Abs(Imag(Prefac)).ge.0.1d-10) THEN 
c             CASE ((0.0d0,-1.0d0),(0.0d0,+1.0d0))
c========================= Propagation if CNV=.TRUE.
        IF(CNV.eqv..TRUE.) THEN
         Do K=1,NL
         SCVECP(k)=ZERO
         SCVECM(K)=ZERO
         DO J=1,NL
         SCVECP(K)=SCVECP(K)+ EV(K,J)*EV(1,J)* EXP(PreFac*EVALS(J)*t)
         SCVECM(K)=SCVECM(K)+ EV(K,J)*EV(1,J)* EXP(-PreFac*EVALS(J)*t)
c         D0(I)=D0(I)+ EV(I,J)*EV(1,J)* EXP(PreFac*EVALS(J)*t)
         ENDDO
         EndDo

        return
        ENDIF
            ENDIF
c===============================================================
            IF (Abs(Real(Prefac)).ge.0.1d-10) THEN 
c             CASE ((0.0d0,-1.0d0),(0.0d0,+1.0d0))
c================= Relaxation ==================================
c         write(6,'(a20,E16.10,a2,E16.10,a2,E16.10)')" Relaxation",
c     .EVALS(1)," ",EVALS(2),"  ",E_state
         IF(print_level.ge.3) then
      write(6,'(F10.5,a20,5F10.5)')E_state," Relaxation",
     .(EVALS(i),i=1,MIN(5,NL))
         ENDIF !PRINT_LEVEL

          IF(DABS(E_state-EVALS(STATE)).le.Error_SIL) CNV=.TRUE. 

c          IF(NL.eq.MaxSil-1) THEN
           DO i=1,NL
c           D0(i)=EV(i,1)
           SCVECP(i)=EV(i,STATE)
           SCVECM(i)=EV(i,STATE)
           ENDDO
c         END IF
c           IF(NL.eq.MaxSil-1) write(6,*)EV(:,1)         

        return
            ENDIF
            
        end subroutine Get_Eigen_SIL
