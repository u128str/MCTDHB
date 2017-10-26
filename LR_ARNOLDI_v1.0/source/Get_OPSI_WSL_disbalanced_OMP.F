c============== This subroutine computes 2body part of the right hand side part of the equations for DPsi/Dt -
c============== usig WSL- based algorithm - one needs to get an access
c============== to all elements of the 2-body density density matrix:
c============== RhoAll(I.J.K.L) I,J,K,L=1..Morb
Corg       SUBROUTINE  Get_OPSI_WSL_based (PSI,OPSI,ZRIJKL,
Corg     .             NDX,NDY,NDZ,Morb,
Corg     .             delta,xlambda0,
Corg     .             Fij,F,weight,ort,NtVec,
Corg     .             FromN,TillN)

      SUBROUTINE  Get_OPSI_WSL_disbalanced_OMP (PSI,OPSI,FromN,TillN,
     &                                                         WSLTime)
      USE SHARED_DIMS
      USE W_INTERPARTICLE
      USE rR_hW
      USE CI_All
      USE omp_lib
      USE DVR_ALL

      IMPLICIT NONE 
      real*8 :: WSLTime
      COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb)  :: PSI
      COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb)  :: OPSI
c      COMPLEX*16, DIMENSION(NDX*NDY*NDZ) :: T1,T2,T3

      COMPLEX*16,  DIMENSION(NDX*NDY*NDZ)  :: t1,t2,t3
      INTEGER :: TID,this_thread,wdim,ierr
      COMPLEX*16, DIMENSION(100) :: sumvout1=Zoner
       
c============================== Timings  =============================
      real*4   start,finish, exec_time , cnk_time       
c=====================================================================
      INTEGER :: K,S,L,Q,I,n1,Nadr,ii,n11,SL,P
c========================= ============================================
      INTEGER  ::  FromN,TillN!,rdim,rdim1
c============================ Data with 2-body density elements =============================
      COMPLEX*16 :: RhoAll,Z,ZDOTC,DZNRM2

      EXTERNAL RhoAll,Nadr,ZDOTC,DZNRM2
      
      COMPLEX*16 WIJKLTEMP(RDIM1)

c============================ RhoAll(K,L,Q,S,ZRIJKL) complex function
c============================ providing required element of density
c============================ matrix, assuming that all unique elements of R_ijkl are packed as (i<=j)<=(k<=l)
c=====================================================================
      OPSI=ZERO
      
      IF(ENRG_EVAL.eqv..TRUE.) WIJKLTEMP=ZERO

c!$OMP PARALLEL DO 
c!$OMP& PRIVATE(P,L,S,Z,n1,n11,K,Q,SL,T1,T2,T3)
c!$OMP& FIRSTPRIVATE(PSI,Morb,ENRG_EVAL,xlambda0,FromN,TillN)
c!$OMP& REDUCTION (+:WIJKLTEMP) REDUCTION(+:OPSI)
c      SLloop: DO SL=1,Morb*(Morb+1)/2
      SLloop: DO SL=FromN,TillN
        P=TERM_INDEX_1B(SL)
        L= INT(P/100)
        S= P-L*100
        CALL  Get_WSL(T3,psi(:,S),psi(:,l),WSLTime)
        Qloop: DO Q=1,Morb
                T1=T3*PSI(:,Q)
                T2=Conjg(T3)*PSI(:,Q)
        Kloop:  DO K=1,Morb
                IF(ENRG_EVAL.eqv..TRUE.) THEN
                   n11= Nadr(K,S,Q,L,Morb)
                   n1= ABS(n11)
                   IF(n11.ge.0)  THEN 
                      WIJKLTEMP(n1)=WIJKLTEMP(n1)+
     .                ZDOTC(NDX*NDY*NDZ,PSI(:,K),1,T1,1)
                   ENDIF
                ENDIF
                Z=RhoAll(K,S,Q,L)
                OPSI(:,K)=OPSI(:,K)+Z*T1 
                IF(L.eq.S) cycle Kloop
c============================ CYCLE is to awoid double counting in the
c============================ case of W11
                IF(ENRG_EVAL.eqv..TRUE.) THEN
                    n11= Nadr(K,L,Q,S,Morb)
                    n1=  ABS(n11)
                    IF(n11.ge.0) THEN
                       WIJKLTEMP(n1)=WIJKLTEMP(n1)+
     .                 ZDOTC(NDX*NDY*NDZ,PSI(:,K),1,T2,1)
                    ENDIF
                ENDIF
                Z=RhoAll(K,L,Q,S)
                OPSI(:,K)=OPSI(:,K)+Z*T2
             END DO Kloop
           END DO Qloop
         END DO SLloop
c!$OMP END PARALLEL DO
      IF(ENRG_EVAL.eqv..TRUE.) WIJKL=WIJKLTEMP
      end subroutine  Get_OPSI_WSL_disbalanced_OMP

