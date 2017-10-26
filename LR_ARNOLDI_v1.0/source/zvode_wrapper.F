      subroutine ZVODE_wrapper(Func,PsiDim,Psi,AbsTime,IntPeriod,
     +           TolError,ErrorCode,ZWORK,Arr)

      INTEGER PsiDim
      INTEGER ML
      INTEGER MU,Arr

      EXTERNAL JEX, JEXB,Func

      DOUBLE COMPLEX Psi(PsiDim), RPAR, WTRU
      DOUBLE PRECISION RWORK(20+PsiDim)
      INTEGER IWORK(30+PsiDim)
      DOUBLE COMPLEX ZWORK(10*PsiDim + (5)*PsiDim)
      DOUBLE PRECISION ABERR, ATOL, AbsTime, IntPeriod,TOUT
      INTEGER ErrorCode
C    OLD WHEELS 
      INTEGER ITOL,ITASK,IOPT
      DOUBLE PRECISION RTOL,TolError,TolErrorArray(PsiDim)
 
      RTOL = 0.d0 
C ITOL   = 1 or 2 according as ATOL (below) is a scalar or array.
      IF (Arr.eq.1) then
         ITOL = 1
      ELSE IF (Arr.eq.2) then
         ITOL = 2
         TolErrorArray=Tolerror
      ELSE
         write(6,*) 'Integration Order must be set to 1 or 2 for ZVODE!'
      ENDIF

      ITASK = 1
      IOPT=0
      ISTATE=1
      IWORK(1) = 1
      IWORK(2) = 1
      RPAR = DCMPLX(0.0D0,1.0D0)
      NERR = 0

      LZW=10*PsiDim + (5)*PsiDim
      LRW=20+PsiDim
      LIW=30+PsiDim

      TOUT=AbsTime+IntPeriod
      IF (ITOL.eq.1) then
         CALL ZVODE(Func,PsiDim,Psi,AbsTime,TOUT,
     +           ITOL,RTOL,TolError,ITASK,ISTATE,IOPT,
     +           ZWORK,LZW,RWORK,LRW,IWORK,LIW,JEX,25)
      ELSE IF (ITOL.eq.2) then
          CALL ZVODE(Func,PsiDim,Psi,AbsTime,TOUT,
     +           ITOL,RTOL,TolErrorArray,ITASK,ISTATE,IOPT,
     +           ZWORK,LZW,RWORK,LRW,IWORK,LIW,JEX,25)
      ENDIF
  
CCCCCCCC TAILS ARE REMOVED
C        ,RPAR,IPAR)
        ErrorCode=ISTATE

      RETURN
      END

      SUBROUTINE JEX (NEQ, T, W, ML, MU, PD, NRPD, RPAR, IPAR)
C     Dense Jacobian routine.
      DOUBLE COMPLEX W(NEQ), PD(NRPD,NEQ), RPAR
      DOUBLE PRECISION T

      write(6,*) 'I am computing a non-existant Jacobian !'
      stop

      RETURN
      END
