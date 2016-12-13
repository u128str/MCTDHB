************************************************************************
*
* SUBROUTINE ODP5CR
*
* Core of the Runge-Kutta integrator of order 5
* (Dormand-Prince formula)
*
* REFERENCE : Hairer/Norsett/Wanner
*             Solving Ordinary Differential Equations I, 2.ed.
*             Springer Verlag, 1993
*
************************************************************************
*
* INPUT:
* ------
* ODEMOD  external subroutine f that calculates the RHS of the
*         differential equation dy/dt = f(y,t).
*         Must have this interface:
*
*         subroutine ODEMOD(t,y,f,mc,mr,mi,ml)
*         real*8     t       in: time of RHS-evaluation
*         complex*16 y(*)    in: current values of the dynamic variables
*         complex*16 f(*)    out: result of RHS-evaluation
*         complex*16 mc(*)   .
*         real*8     mr(*)   . the well-known
*         integer    mi(*)   . MCTDH-arrays
*         logical    ml(*)   .
*
* ndgl    INTEGER     number of dynamic variables
* y       COMPLEX*16  initial values of the dynamic variables
* t       REAL*8      start time of integration
* tend    REAL*8      end time of integration
* hmax    REAL*8      maximum allowed stepsize
* h       REAL*8      initial stepsize (for h=0, we will guess one)
* tol     REAL*8      accuracy parameter
*
* uround  REAL*8      machine precision
* safe    REAL*8      safety factor for stepsize prediction
* beta    REAL*8      parameter for Lund-stabilisation
* facdec  REAL*8      factor for maximal decrease of stepsize
* facinc  REAL*8      factor for maximal increase of stepsize
*
* y1      COMPLEX*16  storage for intermediate result
* k1..k6, COMPLEX*16  storage for Runge-Kutta stages
* ys
* wt      REAL*8      storage for scaling vector
*
* nfcn    INTEGER     number of evaluations of righthand side
* nstep   INTEGER     number of steps
* naccpt  INTEGER     number of accepted steps
* nrejct  INTEGER     number of rejected steps
*
* intgfl  INTEGER     control/error flag
*                     0: this is the first call
*                     1: this is a follow-up call
*
* lpost   LOGICAL     whether to call a post-processing function after
*                     each accepted step
* POSTFN  EXTERNAL    the post-processing function. Interface:
*                     POSTFN(y,mc,mr,mi,ml)
*
*
* OUTPUT:
* -------
* y       COMPLEX*16  final values of dynamic variables
* t       REAL*8      actual time reached in integration
* h       REAL*8      last stepsize
* hopt    REAL*8      last optimal stepsize (not trimmed)
*
* nfcn    INTEGER     number of evaluations of righthand side
* nstep   INTEGER     number of steps
* naccpt  INTEGER     number of accepted steps
* nrejct  INTEGER     number of rejected steps
*
* intgfl  INTEGER     control/error flag
*                        1: OK
*                     else: error!
* ierr    INTEGER     error code
* rerr    REAL*8      error parameter
*
************************************************************************

      subroutine ODP5CR(
     $     ODEMOD, ndgl, y, t, tend, hmax, h, hopt, tol,
     $     uround, safe, beta, facdec, facinc,
     $     y1, k1, k2, k3, k4, k5, k6, ys, wt,
     $     nfcn, nstep, naccpt, nrejct,
     $     intgfl, ierr, rerr,
     $     lpost, POSTFN)

      implicit none

      external   ODEMOD, POSTFN
      integer    ndgl, nfcn, nstep, naccpt, nrejct, intgfl, ierr
      real*8     t, tend, hmax, h, hopt, tol, uround, wt(ndgl),
     $           safe, beta, facdec, facinc, rerr
      complex*16 y(ndgl), y1(ndgl), ys(ndgl), k1(ndgl), k2(ndgl),
     $           k3(ndgl), k4(ndgl), k5(ndgl), k6(ndgl)
      logical    lpost

      integer    stepnr, iasti, iord
      real*8     err, expo1, fac, fac1, facold, facdc1, facin1,
     $           posneg, hnew, told, wgtnrm
      complex*16 ztmp
      logical    lastfl, lrjcfl, lcntfl

      REAL*8
     @ a21,a31,a32,a41,a42,a43,a51,a52,a53,a54,a61,
     @ a62,a63,a64,a65,b1,b3,b4,b5,b6,c2,c3,c4,c5,e1,e3,e4,e5,e6,e7,
     $ ONE,ZERO
  
      PARAMETER( ONE = 1.d0, ZERO = 0.d0 )

C  ----------------------------------
C  --- Coefficients of the method ---
C  ----------------------------------
      PARAMETER (
     @ c2  = 0.2d0,
     @ c3  = 0.3d0,
     @ c4  = 0.8d0,
     @ c5  =  8.d0/9.d0 )

      PARAMETER (
     @ a21 =   0.2d0,
     @ a31 =    3.d0/40.d0,
     @ a32 =    9.d0/40.d0,
     @ a41 =   44.d0/45.d0,
     @ a42 = - 56.d0/15.d0,
     @ a43 =   32.d0/9.d0 )

      PARAMETER (
     @ a51 =  19372.d0/6561.d0,
     @ a52 = -25360.d0/2187.d0,
     @ a53 =  64448.d0/6561.d0,
     @ a54 = -  212.d0/729.d0,
     @ a61 =   9017.d0/3168.d0,
     @ a62 = -  355.d0/33.d0,
     @ a63 =  46732.d0/5247.d0,
     @ a64 =     49.d0/176.d0,
     @ a65 = - 5103.d0/18656.d0 )

      PARAMETER (
     @ b1  =    35.d0/384.d0,
     @ b3  =   500.d0/1113.d0,
     @ b4  =   125.d0/192.d0,
     @ b5  = -2187.d0/6784.d0,
     @ b6  =    11.d0/84.d0 )

      PARAMETER (
     @ e1  =     71.d0/57600.d0,
     @ e3  = -   71.d0/16695.d0,
     @ e4  =     71.d0/1920.d0,
     @ e5  = -17253.d0/339200.d0,
     @ e6  =     22.d0/525.d0,
     @ e7  = -    1.d0/40.d0 )

      SAVE

C  =======================
C  === INITIALIZATIONS ===
C  =======================
      stepnr = 1
      IF (intgfl .EQ. 1) THEN
         err    = .1d0        ! zur Sicherheit
         GOTO 320             ! Zwischenlabel, (kein Sprung in IF-Block)
      ENDIF
C  -----------------------------------------
C  --- Some factors for stepsize-control ---
C  -----------------------------------------
      facold  = 1.d-4
      expo1   = 0.2d0-beta*0.75d0
      facdc1  = ONE/facdec
      facin1  = ONE/facinc
C  ------------------------------
C  --- Some more preparations ---
C  ------------------------------
      posneg = SIGN(ONE,tend-t)
      iasti  = 0
      call zcopy(ndgl,y,1,y1,1)
C  ------------------------
C  --- Initialize flags ---
C  ------------------------
      lastfl = .FALSE.
      lrjcfl = .FALSE.
C  --------------------------------
C  --- Initial model evaluation ---
C  --------------------------------
*  can be omitted since k2=dtpsi
*      CALL ODEMOD(t,y1,k2,mc,mr,mi,ml)
*      nfcn = nfcn+1

C  ===============
C  === RESTART ===
C  ===============
      lastfl = .FALSE.
C  ------------------------------------
C  --- Preparations for ODEHIN-CALL ---
C  ------------------------------------
      hnew   = h
      hmax   = ABS(hmax)
      iord   = 5

      IF (hnew .EQ. ZERO) THEN
         call ODEHIN(ODEMOD, ndgl, t, y1, tend, posneg, k2, k3, k4,
     $        iord, hmax, uround, wt, tol, nfcn, hnew)
      ENDIF
      hnew = ABS(hnew)*posneg

C  ========================================
C  === TOP OF THE MAIN INTEGRATION LOOP ===
C  ========================================
   1  CONTINUE

         lcntfl = .FALSE.
         h      =  hnew
         hopt   =  hnew
         IF (lrjcfl) THEN
            CALL zcopy(ndgl,y,1,y1,1)
         ELSE
C        -----------------------------------
C        ---  Copy result to y, k2 to k1 ---
C        -----------------------------------
            CALL zcopy(ndgl,y1,1, y,1)
            CALL zcopy(ndgl,k2,1,k1,1)
         ENDIF
C     ========================
C     === LEAVE INTEGRATOR ===
C     ========================
         IF ((intgfl .NE. 0) .OR. lastfl ) THEN
C        ---------------------------------
C        --- Leave the core-integrator ---
C        ---------------------------------
            if (lastfl) intgfl = 1
            RETURN
         ENDIF
C     -----------------------------------
C     --- Is the stepsize too small ? ---
C     -----------------------------------
         IF (0.1d0*ABS(h) .LE. ABS(t)*uround) THEN
            intgfl = -51
            rerr   =  t
            ierr   =  LOG10(h)
            GOTO 1
         ENDIF
C     -----------------------------------------------
C     --- Adjust h, if integrator is near the end ---
C     -----------------------------------------------
         IF ((t+1.01d0*h-tend)*posneg .GT. ZERO) THEN
            h      =  tend - t
            lastfl = .TRUE.
         ENDIF
C     ===============================================
C     === CALCULATE THE 6 STAGES FOR THE NEW STEP ===
C     ===============================================
         nstep = nstep + 1
C     ---------------
C     --- Stage 1 ---
C     ---------------
         ztmp = DCMPLX(a21*h)
         CALL zaxpy(ndgl,ztmp,k1,1,y1,1)

         CALL ODEMOD(t+c2*h,y1,k2)
C     ---------------
C     --- Stage 2 ---
C     ---------------
         CALL zcopy(ndgl,     k1,1,y1,1)
         ztmp = DCMPLX(a31*h)
         CALL zscal(ndgl,ztmp,     y1,1)
         ztmp = DCMPLX(a32*h)
         CALL zaxpy(ndgl,ztmp,k2,1,y1,1)
         ztmp = DCMPLX(ONE  )
         CALL zaxpy(ndgl,ztmp, y,1,y1,1)

         CALL ODEMOD(t+c3*h,y1,k3)
C     ---------------
C     --- Stage 3 ---
C     ---------------
         CALL zcopy(ndgl,     k1,1,y1,1)
         ztmp = DCMPLX(a41*h)
         CALL zscal(ndgl,ztmp,     y1,1)
         ztmp = DCMPLX(a42*h)
         CALL zaxpy(ndgl,ztmp,k2,1,y1,1)
         ztmp = DCMPLX(a43*h)
         CALL zaxpy(ndgl,ztmp,k3,1,y1,1)
         ztmp = DCMPLX(ONE  )
         CALL zaxpy(ndgl,ztmp, y,1,y1,1)

         CALL ODEMOD(t+c4*h,y1,k4)
C     ---------------
C     --- Stage 4 ---
C     ---------------
         CALL zcopy(ndgl,     k1,1,y1,1)
         ztmp = DCMPLX(a51*h)
         CALL zscal(ndgl,ztmp,     y1,1)
         ztmp = DCMPLX(a52*h)
         CALL zaxpy(ndgl,ztmp,k2,1,y1,1)
         ztmp = DCMPLX(a53*h)
         CALL zaxpy(ndgl,ztmp,k3,1,y1,1)
         ztmp = DCMPLX(a54*h)
         CALL zaxpy(ndgl,ztmp,k4,1,y1,1)
         ztmp = DCMPLX(ONE  )
         CALL zaxpy(ndgl,ztmp, y,1,y1,1)

         CALL ODEMOD(t+c5*h,y1,k5)
C     ---------------
C     --- Stage 5 ---
C     ---------------
         CALL zcopy(ndgl,     k1,1,ys,1)
         ztmp = DCMPLX(a61*h)
         CALL zscal(ndgl,ztmp,     ys,1)
         ztmp = DCMPLX(a62*h)
         CALL zaxpy(ndgl,ztmp,k2,1,ys,1)
         ztmp = DCMPLX(a63*h)
         CALL zaxpy(ndgl,ztmp,k3,1,ys,1)
         ztmp = DCMPLX(a64*h)
         CALL zaxpy(ndgl,ztmp,k4,1,ys,1)
         ztmp = DCMPLX(a65*h)
         CALL zaxpy(ndgl,ztmp,k5,1,ys,1)
         ztmp = DCMPLX(ONE  )
         CALL zaxpy(ndgl,ztmp, y,1,ys,1)

         CALL ODEMOD(t+h,ys,k6)
C     ---------------
C     --- Stage 6 ---
C     ---------------
         CALL zcopy(ndgl,     k1,1,y1,1)
         ztmp = DCMPLX(b1*h)
         CALL zscal(ndgl,ztmp,     y1,1)
         ztmp = DCMPLX(b3*h)
         CALL zaxpy(ndgl,ztmp,k3,1,y1,1)
         ztmp = DCMPLX(b4*h)
         CALL zaxpy(ndgl,ztmp,k4,1,y1,1)
         ztmp = DCMPLX(b5*h)
         CALL zaxpy(ndgl,ztmp,k5,1,y1,1)
         ztmp = DCMPLX(b6*h)
         CALL zaxpy(ndgl,ztmp,k6,1,y1,1)
         ztmp = DCMPLX(ONE )
         CALL zaxpy(ndgl,ztmp, y,1,y1,1)

         CALL ODEMOD(t+h,y1,k2)
         nfcn = nfcn + 6

C     -----------------------
C     --- The error-stage ---
C     -----------------------
         ztmp = DCMPLX(e4*h)
         CALL zscal(ndgl,ztmp,     k4,1)
         ztmp = DCMPLX(e1*h)
         CALL zaxpy(ndgl,ztmp,k1,1,k4,1)
         ztmp = DCMPLX(e3*h)
         CALL zaxpy(ndgl,ztmp,k3,1,k4,1)
         ztmp = DCMPLX(e5*h)
         CALL zaxpy(ndgl,ztmp,k5,1,k4,1)
         ztmp = DCMPLX(e6*h)
         CALL zaxpy(ndgl,ztmp,k6,1,k4,1)
         ztmp = DCMPLX(e7*h)
         CALL zaxpy(ndgl,ztmp,k2,1,k4,1)
C     --------------------------
C     --- Estimate the error ---
C     --------------------------
         CALL SETEWV(ndgl,y,y1,tol,uround,wt)
         err = wgtnrm(ndgl,k4,wt)
C     --------------------
C     --- Compute hnew ---
C     --------------------
         fac1  = err**expo1
C     --------------------------
C     --- Lund-stabilization ---
C     --------------------------
         fac   = fac1/facold**beta
C     ----------------------------------------------
C     --- We require  facdec <= hnew/h <= facinc ---
C     ----------------------------------------------
         fac   = MAX(facin1,MIN(facdc1,fac/safe))
         hnew  = h/fac
C     ========================
C     === STEP IS ACCEPTED ===
C     ========================
 320     CONTINUE
         IF (err .LE. ONE) THEN
            IF (intgfl .EQ. 1) THEN
               intgfl = 0
               GOTO 32
            ENDIF
            lrjcfl = .FALSE.
            facold =  MAX(err,1.0d-4)
            naccpt =  naccpt + 1
C        -------------------
C        --- Adjust time ---
C        -------------------
            told = t
            t    = t + h

   32       CONTINUE
*           write into MCTDH's STEPS-file
!            CALL WriteStep(stepnr,1,h,err,told)
*           call post-processing function if required
            if (lpost) call POSTFN(y1)
            stepnr = stepnr+1

C        ==========================
C        === PREPARE NEXT STEP  ===
C        ==========================
C        ---------------------------
C        --- Calculation of hnew ---
C        ---------------------------
            IF (ABS(hnew) .GT. hmax) hnew = posneg*hmax
            IF (lrjcfl) hnew = posneg*MIN(ABS(hnew),ABS(h))
C     ========================
C     === STEP IS REJECTED ===
C     ========================
         ELSE
*           write into MCTDH's STEPS-file
!            CALL WriteStep(stepnr,0,h,err,t)
            stepnr = stepnr+1
C        ---------------------------
C        --- Calculation of hnew ---
C        ---------------------------
            hnew   = h/MIN(facdc1,fac1/safe)
C        ---------------------------------------
C        --- Adjust the flags and statistics ---
C        ---------------------------------------
            lrjcfl = .TRUE.
            lastfl = .FALSE.
*           IF (naccpt .GE. 1) nrejct = nrejct + 1
            nrejct = nrejct+1
         ENDIF
C  ====================================================
C  === HERE IS THE END OF THE MAIN INTEGRATION-LOOP ===
C  ====================================================
      GOTO 1

      END                                                    ! of ODP5CR
