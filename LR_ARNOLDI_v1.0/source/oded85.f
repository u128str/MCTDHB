************************************************************************
*
* SUBROUTINE ODP5CR
*
* Core of the Runge-Kutta integrator DoPri853.
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
* k1..k10 COMPLEX*16  storage for Runge-Kutta stages
* kerr,
* kerr2   COMPLEX*16  storage for the error stages
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

      subroutine OD85CR(
     $     ODEMOD, ndgl, y, t, tend, hmax, h, hopt, tol,
     $     uround, safe, beta, facdec, facinc,
     $     y1, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10,
     $     kerr, kerr2, wt, nfcn, nstep, naccpt, nrejct,
     $     intgfl, ierr, rerr,
     $     lpost, POSTFN)

      implicit none

      external   ODEMOD, POSTFN
      integer    ndgl, nfcn, nstep, naccpt, nrejct, intgfl, ierr
      real*8     t, tend, hmax, h, hopt, tol, uround, wt(ndgl),
     $           safe, beta, facdec, facinc, rerr
      complex*16 y(ndgl), y1(ndgl), k1(ndgl), k2(ndgl),
     $           k3(ndgl), k4(ndgl), k5(ndgl), k6(ndgl), k7(ndgl),
     $           k8(ndgl), k9(ndgl), k10(ndgl),
     $           kerr(ndgl), kerr2(ndgl)
      logical    lpost

      integer    stepnr, iasti, iord
      real*8     err, expo1, fac, fac1, facold, facdc1, facin1,
     $           posneg, hnew, told, wgtnrm, err2, deno
      complex*16 ztmp
      logical    lastfl, lrjcfl, lcntfl

      REAL*8
     @ a21,a31,a32,a41,a43,a51,a53,a54,a61,a64,a65,a71,a74,a75,a76,
     @ a81,a84,a85,a86,a87,a91,a94,a95,a96,a97,a98,a101,a104,a105,a106,
     @ a107,a108,a109,a111,a114,a115,a116,a117,a118,a119,a1110,a121,
     @ a124,a125,a126,a127,a128,a129,a1210,a1211,b1,b6,b7,
     @ b8,b9,b10,b11,b12,bh1,bh2,bh3,c2,c3,c4,c5,c6,c7,c8,c9,c10,
     @ c11,e1,e6,e7,e8,e9,e10,e11,e12,
     $ ONE,ZERO

      PARAMETER( ONE = 1.d0, ZERO = 0.d0 )

C  ----------------------------------
C  --- Coefficients of the method ---
C  ----------------------------------
      PARAMETER (
     @ c2  = 0.526001519587677318785587544488D-01,
     @ c3  = 0.789002279381515978178381316732D-01,
     @ c4  = 0.118350341907227396726757197510D+00,
     @ c5  = 0.281649658092772603273242802490D+00,
     @ c6  = 0.333333333333333333333333333333D+00,
     @ c7  = 0.25D+00,
     @ c8  = 0.307692307692307692307692307692D+00,
     @ c9  = 0.651282051282051282051282051282D+00,
     @ c10 = 0.6D+00,
     @ c11 = 0.857142857142857142857142857142D+00 )

      PARAMETER (
     @ b1 =   5.42937341165687622380535766363D-02,
     @ b6 =   4.45031289275240888144113950566D+00,
     @ b7 =   1.89151789931450038304281599044D+00,
     @ b8 =  -5.8012039600105847814672114227D+00,
     @ b9 =   3.1116436695781989440891606237D-01,
     @ b10 = -1.52160949662516078556178806805D-01,
     @ b11 =  2.01365400804030348374776537501D-01,
     @ b12 =  4.47106157277725905176885569043D-02 )

      PARAMETER (
     @ bh1 = 0.244094488188976377952755905512D+00,
     @ bh2 = 0.733846688281611857341361741547D+00,
     @ bh3 = 0.220588235294117647058823529412D-01 )

      PARAMETER (
     @ e1 =  0.1312004499419488073250102996D-01,
     @ e6 = -0.1225156446376204440720569753D+01,
     @ e7 = -0.4957589496572501915214079952D+00,
     @ e8 =  0.1664377182454986536961530415D+01,
     @ e9 = -0.3503288487499736816886487290D+00,
     @ e10 =  0.3341791187130174790297318841D+00,
     @ e11 =  0.8192320648511571246570742613D-01,
     @ e12 = -0.2235530786388629525884427845D-01 )

      PARAMETER (
     @ a21 =    5.26001519587677318785587544488D-02,
     @ a31 =    1.97250569845378994544595329183D-02,
     @ a32 =    5.91751709536136983633785987549D-02,
     @ a41 =    2.95875854768068491816892993775D-02,
     @ a43 =    8.87627564304205475450678981324D-02,
     @ a51 =    2.41365134159266685502369798665D-01,
     @ a53 =   -8.84549479328286085344864962717D-01,
     @ a54 =    9.24834003261792003115737966543D-01,
     @ a61 =    3.7037037037037037037037037037D-02,
     @ a64 =    1.70828608729473871279604482173D-01,
     @ a65 =    1.25467687566822425016691814123D-01,
     @ a71 =    3.7109375D-02,
     @ a74 =    1.70252211019544039314978060272D-01,
     @ a75 =    6.02165389804559606850219397283D-02,
     @ a76 =   -1.7578125D-02 )

      PARAMETER (
     @ a81 =    3.70920001185047927108779319836D-02,
     @ a84 =    1.70383925712239993810214054705D-01,
     @ a85 =    1.07262030446373284651809199168D-01,
     @ a86 =   -1.53194377486244017527936158236D-02,
     @ a87 =    8.27378916381402288758473766002D-03,
     @ a91 =    6.24110958716075717114429577812D-01,
     @ a94 =   -3.36089262944694129406857109825D+00,
     @ a95 =   -8.68219346841726006818189891453D-01,
     @ a96 =    2.75920996994467083049415600797D+01,
     @ a97 =    2.01540675504778934086186788979D+01,
     @ a98 =   -4.34898841810699588477366255144D+01,
     @ a101 =   4.77662536438264365890433908527D-01,
     @ a104 =  -2.48811461997166764192642586468D+00,
     @ a105 =  -5.90290826836842996371446475743D-01,
     @ a106 =   2.12300514481811942347288949897D+01,
     @ a107 =   1.52792336328824235832596922938D+01,
     @ a108 =  -3.32882109689848629194453265587D+01,
     @ a109 =  -2.03312017085086261358222928593D-02 )

      PARAMETER (
     @ a111 =  -9.3714243008598732571704021658D-01,
     @ a114 =   5.18637242884406370830023853209D+00,
     @ a115 =   1.09143734899672957818500254654D+00,
     @ a116 =  -8.14978701074692612513997267357D+00,
     @ a117 =  -1.85200656599969598641566180701D+01,
     @ a118 =   2.27394870993505042818970056734D+01,
     @ a119 =   2.49360555267965238987089396762D+00,
     @ a1110 = -3.0467644718982195003823669022D+00,
     @ a121 =   2.27331014751653820792359768449D+00,
     @ a124 =  -1.05344954667372501984066689879D+01,
     @ a125 =  -2.00087205822486249909675718444D+00,
     @ a126 =  -1.79589318631187989172765950534D+01,
     @ a127 =   2.79488845294199600508499808837D+01,
     @ a128 =  -2.85899827713502369474065508674D+00,
     @ a129 =  -8.87285693353062954433549289258D+00,
     @ a1210 =  1.23605671757943030647266201528D+01,
     @ a1211 =  6.43392746015763530355970484046D-01 )

      SAVE
!       write(6,*) 'IN RK8 CORE BEFORE',nfcn, nstep, naccpt, nrejct,
!     $     intgfl, ierr, rerr,
!     $     lpost



C  =======================
C  === INITIALIZATIONS ===
C  =======================
      stepnr = 1
      IF (intgfl .EQ. 1) THEN
         err    = .1d0
         GOTO 320
      ENDIF
C  -----------------------------------------
C  --- Some factors for stepsize-control ---
C  -----------------------------------------
      facold  = 1.d-4  
      expo1   = 0.125d0-beta*0.75d0
      facdc1  = ONE/facdec
      facin1  = ONE/facinc
C  ------------------------------
C  --- Some more preparations ---
C  ------------------------------
      posneg = SIGN(ONE,tend-t)
      iasti  = 0
      CALL zcopy(ndgl,y,1,k5,1)
C  ------------------------
C  --- Initialize flags ---
C  ------------------------
      lastfl = .FALSE. 
      lrjcfl = .FALSE.    
C  --------------------------------
C  --- Initial model evaluation ---
C  --------------------------------
*  can be omitted since k4 already contains dtpsi
*      CALL ODEMOD(t,k5,k4,mc,mr,mi,ml)
*      nfcn = nfcn+1       

C  ===============
C  === RESTART ===
C  ===============
      lastfl = .FALSE.
C  ------------------------------------
C  --- Preparations for ODEHIN-call ---  
C  ------------------------------------
      hnew   = h
      hmax   = ABS(hmax)     
      iord   = 8

      IF (hnew .EQ. ZERO) THEN     
         call ODEHIN(ODEMOD, ndgl, t, k5, tend, posneg, k4, k2, k3,
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
         IF ( .NOT. lrjcfl ) THEN
C        -----------------------------------
C        ---  Copy result to y, k4 to k1 ---
C        -----------------------------------
            CALL zcopy(ndgl,k5,1, y,1)
            CALL zcopy(ndgl,k4,1,k1,1)
         ENDIF
C     ========================
C     === LEAVE INTEGRATOR ===
C     ========================
         IF ((intgfl .NE. 0) .OR. lastfl ) THEN
C        ---------------------------------
C        --- Leave the core-integrator ---
C        ---------------------------------   
!       write(6,*) 'IN RK8 CORE AFTER',nfcn, nstep, naccpt, nrejct,
!     $     intgfl, ierr, rerr,
!     $     lpost


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
C     ================================================
C     === CALCULATE THE 12 STAGES FOR THE NEW STEP === 
C     ================================================
         nstep = nstep + 1   
C     ---------------
C     --- Stage 1 ---
C     --------------- 
         CALL zcopy(ndgl,      y,1,y1,1)
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
         ztmp = DCMPLX(ONE)
         CALL zaxpy(ndgl,ztmp, y,1,y1,1)

         CALL ODEMOD(t+c3*h,y1,k3)
C     ---------------
C     --- Stage 3 ---
C     ---------------     
         CALL zcopy(ndgl,     k1,1,y1,1) 
         ztmp = DCMPLX(a41*h)
         CALL zscal(ndgl,ztmp,     y1,1)
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
         CALL zcopy(ndgl,     k1,1,y1,1) 
         ztmp = DCMPLX(a61*h)
         CALL zscal(ndgl,ztmp,     y1,1)
         ztmp = DCMPLX(a64*h)
         CALL zaxpy(ndgl,ztmp,k4,1,y1,1)  
         ztmp = DCMPLX(a65*h)
         CALL zaxpy(ndgl,ztmp,k5,1,y1,1)
         ztmp = DCMPLX(ONE  )
         CALL zaxpy(ndgl,ztmp, y,1,y1,1)

         CALL ODEMOD(t+c6*h,y1,k6)
C     ---------------
C     --- Stage 6 ---
C     ---------------  
         CALL zcopy(ndgl,     k1,1,y1,1) 
         ztmp = DCMPLX(a71*h)
         CALL zscal(ndgl,ztmp,     y1,1)
         ztmp = DCMPLX(a74*h)
         CALL zaxpy(ndgl,ztmp,k4,1,y1,1) 
         ztmp = DCMPLX(a75*h)
         CALL zaxpy(ndgl,ztmp,k5,1,y1,1)  
         ztmp = DCMPLX(a76*h)
         CALL zaxpy(ndgl,ztmp,k6,1,y1,1)
         ztmp = DCMPLX(ONE  )
         CALL zaxpy(ndgl,ztmp, y,1,y1,1)

         CALL ODEMOD(t+c7*h,y1,k7)
C     ---------------
C     --- Stage 7 ---
C     ---------------  
         CALL zcopy(ndgl,     k1,1,y1,1) 
         ztmp = DCMPLX(a81*h)
         CALL zscal(ndgl,ztmp,     y1,1)
         ztmp = DCMPLX(a84*h)
         CALL zaxpy(ndgl,ztmp,k4,1,y1,1) 
         ztmp = DCMPLX(a85*h)
         CALL zaxpy(ndgl,ztmp,k5,1,y1,1)  
         ztmp = DCMPLX(a86*h)
         CALL zaxpy(ndgl,ztmp,k6,1,y1,1)
         ztmp = DCMPLX(a87*h)
         CALL zaxpy(ndgl,ztmp,k7,1,y1,1) 
         ztmp = DCMPLX(ONE  )
         CALL zaxpy(ndgl,ztmp, y,1,y1,1)

         CALL ODEMOD(t+c8*h,y1,k8)
C     ---------------
C     --- Stage 8 ---
C     ---------------  
         CALL zcopy(ndgl,     k1,1,y1,1) 
         ztmp = DCMPLX(a91*h)
         CALL zscal(ndgl,ztmp,     y1,1)
         ztmp = DCMPLX(a94*h)
         CALL zaxpy(ndgl,ztmp,k4,1,y1,1) 
         ztmp = DCMPLX(a95*h)
         CALL zaxpy(ndgl,ztmp,k5,1,y1,1)  
         ztmp = DCMPLX(a96*h)
         CALL zaxpy(ndgl,ztmp,k6,1,y1,1)
         ztmp = DCMPLX(a97*h)
         CALL zaxpy(ndgl,ztmp,k7,1,y1,1) 
         ztmp = DCMPLX(a98*h)
         CALL zaxpy(ndgl,ztmp,k8,1,y1,1) 
         ztmp = DCMPLX(ONE  )
         CALL zaxpy(ndgl,ztmp, y,1,y1,1)

         CALL ODEMOD(t+c9*h,y1,k9)
C     ---------------
C     --- Stage 9 ---
C     --------------- 
         CALL zcopy(ndgl,     k1,1,y1,1) 
         ztmp = DCMPLX(a101*h)
         CALL zscal(ndgl,ztmp,     y1,1)
         ztmp = DCMPLX(a104*h)
         CALL zaxpy(ndgl,ztmp,k4,1,y1,1) 
         ztmp = DCMPLX(a105*h)
         CALL zaxpy(ndgl,ztmp,k5,1,y1,1)  
         ztmp = DCMPLX(a106*h)
         CALL zaxpy(ndgl,ztmp,k6,1,y1,1)
         ztmp = DCMPLX(a107*h)
         CALL zaxpy(ndgl,ztmp,k7,1,y1,1) 
         ztmp = DCMPLX(a108*h)
         CALL zaxpy(ndgl,ztmp,k8,1,y1,1) 
         ztmp = DCMPLX(a109*h)
         CALL zaxpy(ndgl,ztmp,k9,1,y1,1) 
         ztmp = DCMPLX(ONE   )
         CALL zaxpy(ndgl,ztmp, y,1,y1,1)

         CALL ODEMOD(t+c10*h,y1,k10)
C     ----------------
C     --- Stage 10 ---
C     ---------------- 
         CALL zcopy(ndgl,      k1,1,y1,1) 
         ztmp = DCMPLX(a111*h )
         CALL zscal(ndgl,ztmp,      y1,1)
         ztmp = DCMPLX(a114*h )
         CALL zaxpy(ndgl,ztmp, k4,1,y1,1) 
         ztmp = DCMPLX(a115*h )
         CALL zaxpy(ndgl,ztmp, k5,1,y1,1)  
         ztmp = DCMPLX(a116*h )
         CALL zaxpy(ndgl,ztmp, k6,1,y1,1)
         ztmp = DCMPLX(a117*h )
         CALL zaxpy(ndgl,ztmp, k7,1,y1,1) 
         ztmp = DCMPLX(a118*h )
         CALL zaxpy(ndgl,ztmp, k8,1,y1,1) 
         ztmp = DCMPLX(a119*h )
         CALL zaxpy(ndgl,ztmp, k9,1,y1,1) 
         ztmp = DCMPLX(a1110*h)
         CALL zaxpy(ndgl,ztmp,k10,1,y1,1) 
         ztmp = DCMPLX(ONE)
         CALL zaxpy(ndgl,ztmp,  y,1,y1,1) 

         CALL ODEMOD(t+c11*h,y1,k2)
C     ----------------
C     --- Stage 11 ---
C     ----------------  
         CALL zcopy(ndgl,      k1,1,y1,1) 
         ztmp = DCMPLX(a121*h )
         CALL zscal(ndgl,ztmp,      y1,1)
         ztmp = DCMPLX(a124*h )
         CALL zaxpy(ndgl,ztmp, k4,1,y1,1) 
         ztmp = DCMPLX(a125*h )
         CALL zaxpy(ndgl,ztmp, k5,1,y1,1)  
         ztmp = DCMPLX(a126*h )
         CALL zaxpy(ndgl,ztmp, k6,1,y1,1)
         ztmp = DCMPLX(a127*h )
         CALL zaxpy(ndgl,ztmp, k7,1,y1,1) 
         ztmp = DCMPLX(a128*h )
         CALL zaxpy(ndgl,ztmp, k8,1,y1,1) 
         ztmp = DCMPLX(a129*h )
         CALL zaxpy(ndgl,ztmp, k9,1,y1,1) 
         ztmp = DCMPLX(a1210*h)
         CALL zaxpy(ndgl,ztmp,k10,1,y1,1) 
         ztmp = DCMPLX(a1211*h)
         CALL zaxpy(ndgl,ztmp, k2,1,y1,1) 
         ztmp = DCMPLX(ONE)
         CALL zaxpy(ndgl,ztmp,  y,1,y1,1)

         CALL ODEMOD(t+h,y1,k3)
         nfcn = nfcn + 11
C     ----------------
C     --- Stage 12 ---
C     ----------------  
         CALL zcopy(ndgl,      k1,1,k4,1) 
         ztmp = DCMPLX(b1 )
         CALL zscal(ndgl,ztmp,      k4,1)
         ztmp = DCMPLX(b6 )
         CALL zaxpy(ndgl,ztmp, k6,1,k4,1)
         ztmp = DCMPLX(b7 )
         CALL zaxpy(ndgl,ztmp, k7,1,k4,1) 
         ztmp = DCMPLX(b8 )
         CALL zaxpy(ndgl,ztmp, k8,1,k4,1) 
         ztmp = DCMPLX(b9 )
         CALL zaxpy(ndgl,ztmp, k9,1,k4,1) 
         ztmp = DCMPLX(b10)
         CALL zaxpy(ndgl,ztmp,k10,1,k4,1) 
         ztmp = DCMPLX(b11)
         CALL zaxpy(ndgl,ztmp, k2,1,k4,1) 
         ztmp = DCMPLX(b12)
         CALL zaxpy(ndgl,ztmp, k3,1,k4,1) 

         CALL zcopy(ndgl,      y,1,k5,1)
         ztmp = DCMPLX(h)
         CALL zaxpy(ndgl,ztmp,k4,1,k5,1)

C     -----------------------
C     --- The error-stage ---  
C     -----------------------
         CALL zcopy(ndgl,      k1,1,kerr,1) 
         ztmp = DCMPLX(e1 )
         CALL zscal(ndgl,ztmp,      kerr,1)
         ztmp = DCMPLX(e6 )
         CALL zaxpy(ndgl,ztmp, k6,1,kerr,1)
         ztmp = DCMPLX(e7 )
         CALL zaxpy(ndgl,ztmp, k7,1,kerr,1) 
         ztmp = DCMPLX(e8 )
         CALL zaxpy(ndgl,ztmp, k8,1,kerr,1) 
         ztmp = DCMPLX(e9 )
         CALL zaxpy(ndgl,ztmp, k9,1,kerr,1) 
         ztmp = DCMPLX(e10)
         CALL zaxpy(ndgl,ztmp,k10,1,kerr,1) 
         ztmp = DCMPLX(e11)
         CALL zaxpy(ndgl,ztmp, k2,1,kerr,1) 
         ztmp = DCMPLX(e12)
         CALL zaxpy(ndgl,ztmp, k3,1,kerr,1) 

C     +-------------------------+
C     +--   Result     -> k5  --+
C     +--   Last stage -> y1  --+
C     +-------------------------+
C     +--   Extra stage -> k4 --+
C     +-------------------------+
C     +--   Later on:         --+
C     +--   RHS(k5)    -> k4  --+
C     +--   RHS(y1)    -> k3  --+ 
C     +-------------------------+
C     ------------------------
C     --- Error estimation ---
C     ------------------------
         CALL SETEWV(ndgl,y,k5,tol,uround,wt)
         err   = wgtnrm(ndgl,kerr,wt) 

         CALL zcopy(ndgl,     k4,1,kerr2,1)
         ztmp = DCMPLX(-bh1)
         CALL zaxpy(ndgl,ztmp,k1,1,kerr2,1)
         ztmp = DCMPLX(-bh2)
         CALL zaxpy(ndgl,ztmp,k9,1,kerr2,1)
         ztmp = DCMPLX(-bh3)
         CALL zaxpy(ndgl,ztmp,k3,1,kerr2,1) 

         err2  = wgtnrm(ndgl,kerr2,wt)  
         err   = err**2
         deno  = err + 0.01d0*err2**2
         IF (deno .LE. ZERO) deno = ONE 
         err   = h*(err/SQRT(deno))

C     ---------------------------
C     --- Computation of hnew --- 
C     ---------------------------
         fac1  = err**expo1
C     --------------------------
C     --- LUND-stabilization ---
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
C        -------------------
C        --- Adjust time ---
C        -------------------
            told   = t
            t      = t + h 
C        -----------------------------
C        --- New model evaluation  ---
C        -----------------------------
            CALL ODEMOD(t,k5,k4)
            nfcn   = nfcn + 1    
            facold = MAX(err,1.0d-4)
            naccpt = naccpt + 1

   32       CONTINUE
*           write into MCTDH's STEPS-file
!            CALL WriteStep(stepnr,1,h,err,told)
*           call post-processor if required
            if (lpost) call POSTFN(k5)
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

      END                                                    ! of OD85CR


      subroutine dummypost(psi)

      implicit none

      complex*16 psi(*)

*     do nothing
      return

      end
