************************************************************************
*
* SubRoutine ODEHIN
*
* computes an initial stepsize guess for integration
*
* REFERENCE : Hairer/Norsett/Wanner
*             Solving Ordinary Differential Equations I, 2.ed.
*             Springer Verlag, 1993
*
************************************************************************
*
* INPUT:
* ------
* ODEMOD  EXTERNAL    function to evaluate the RHS of the ODE
* ndgl    INTEGER     dimension of differential equation system
* t       REAL*8      initial time
* y       COMPLEX*16  initial values of dynamic variables
* tend    REAL*8      final time of integration 
* posneg  REAL*8      direction of integration
* dy0     COMPLEX*16  initial derivatives of dynamics variables
* dy1,y1  COMPLEX*16  additional storage for EULER step
* iord    INTEGER     order of integration method to be used
* hmax    REAL*8      maximum permitted stepsize
* small   REAL*8      underflow indicator
* wt      REAL*8      storage for weight vector
* tol     REAL*8      accuracy parameter
* nfcn    INTEGER     statistics variable
*
* OUPUT:
* ------
* nfcn    INTEGER     statistics variable
* hini    REAL*8      first stepsize to be tried in integration
*
************************************************************************

      subroutine ODEHIN(ODEMOD, ndgl, t, y, tend, posneg, dy0, dy1,
     $     y1, iord, hmax, small, wt, tol, nfcn, hini)

      implicit none

      external   ODEMOD
      integer    ndgl, iord, nfcn
      real*8     t, tend, posneg, hmax, small, wt(ndgl), tol, hini
      complex*16 y(ndgl), dy0(ndgl), dy1(ndgl), y1(ndgl)
      
      real*8     wgtnrm, wgtdst
      real*8     ynrm, dynrm, h, hmx, ddynrm, dnrmmx, htry, ONE
      complex*16 ztmp

      PARAMETER( ONE = 1.d0 )

C --------------------
C --- PREPARATIONS ---
C --------------------
      call SETEWV(ndgl,y,y,tol,small,wt)
      ynrm  = wgtnrm(ndgl, y, wt)
      dynrm = wgtnrm(ndgl, dy0, wt)
C -------------------------------------------------------
C --- COMPUTE STEPSIZE FOR EXPLICIT EULER SUCH THAT   ---
C --- THE INCREMENT IS SMALL COMPARED TO THE SOLUTION ---
C -------------------------------------------------------
      IF ((dynrm .LE. small) .OR. (ynrm .LE. small)) THEN
         h = 1.0d-6
      ELSE
         h = MAX(hmax*small,ABS((ynrm/dynrm)*0.01d0))
      ENDIF
      hmx = MIN(ABS(tend-t),hmax)
      h   = MIN(h,hmx)
      h   = h*posneg
C  --------------------------------------
C  --- PERFORM AN EXPLICIT EULER STEP ---
C  --------------------------------------
      CALL zcopy(ndgl,y,1,y1,1)
      ztmp = DCMPLX(h)
      CALL zaxpy(ndgl,ztmp,dy0,1,y1,1)
C  ----------------
C  --- EVALUATE ---
C  ----------------
      CALL ODEMOD(t+h,y1,dy1)
      nfcn=nfcn+1
C  ------------------------------------------------------
C  --- ESTIMATE THE SECOND DERIVATIVE OF THE SOLUTION ---
C  ------------------------------------------------------
      CALL SETEWV(ndgl,y,y1,tol,small,wt)
      ddynrm = wgtdst(ndgl,dy0,dy1,wt)
      ddynrm = ddynrm/h
C ---------------------------------
C --- COMPUTE INITIAL STEP SIZE ---
C ---------------------------------
      dnrmmx = MAX(ABS(ddynrm),dynrm)
      IF (dnrmmx .LE. 1.d-15) THEN
         htry = MAX(1.0d-6,ABS(h)*1.0d-3)
      ELSE
         htry = (0.01d0/dnrmmx)**(ONE/iord)
      ENDIF
      hini = MIN(1.0d2*h,htry,hmx)
      hini = hini*posneg

      RETURN
      END                                                   !  of ODEHIN




************************************************************************
*
* SubRoutine SETEWV
*
* set the error weight vector according to the given tolerance
*
************************************************************************
*
* INPUT:
* ------
* ndim    INTEGER     dimension of the vectors
* y0,y1   COMPLEX*16  scaling vectors
* tol     REAL*8      tolerance parameter
* small   REAL*8      underflow indicator
*
* OUTPUT:
* -------
* wt      REAL*8      error weight vector
*
************************************************************************

      subroutine SETEWV(ndim, y0, y1, tol, small, wt)

      implicit none

      integer    ndim, i
      real*8     tol, small, wt(ndim)
      complex*16 y0(ndim), y1(ndim)

C -----------------------------------
C --- BOTH TOLERANCES ARE SCALARS ---
C -----------------------------------
      do i = 1,ndim
         wt(i) = tol + tol*MAX(ABS(y0(i)), ABS(y1(i)), small)
      enddo

      return
      end                                                    ! of SETEWV




************************************************************************
*
* Function WGTNRM
*
* computes the weighted norm of a vector
*
************************************************************************
*
* INPUT:
* ------
* ndim    INTEGER     dimension of the vectors
* vec     COMPLEX*16  vector, whose norm is to be computed
* wt      REAL*8      vector of weights
*
* OUTPUT :
* --------
* wgtnrm  REAL*8     function value
*
************************************************************************

      REAL*8 FUNCTION wgtnrm(ndim, vec, wt)

      IMPLICIT NONE

      INTEGER    i,ndim
      REAL*8     tmp,wt(ndim)
      COMPLEX*16 vec(ndim)

      tmp = 0.d0
      DO i = 1,ndim
*        tmp = tmp + (dimag(vec(i))**2 + dble(vec(i))**2) / (wt(i)**2)
         tmp = tmp + (abs(vec(i))/wt(i))**2
      ENDDO
      wgtnrm = SQRT(tmp/ndim)

      RETURN
      END                                                    ! of WGTNRM





************************************************************************
*
* Function WGTDST
*
* computes the weighted distance of two vectors
*
************************************************************************
*
* INPUT:
* ------
* ndim       INTEGER     dimension of the vectors
* vec1,vec2  COMPLEX*16  vectors, whose distance is to be computed
* wt         REAL*8      vector of weights
*
* OUTPUT :
* --------
* wgtdst     REAL*8      function value
*
************************************************************************

      REAL*8 FUNCTION wgtdst(ndim, vec1, vec2, wt)

      IMPLICIT NONE

      INTEGER    i,ndim
      REAL*8     tmp,wt(ndim)
      COMPLEX*16 vec1(ndim),vec2(ndim)

      tmp = 0.d0
      DO i = 1,ndim
         tmp = tmp + (abs(vec1(i)-vec2(i))/wt(i))**2
      ENDDO
      wgtdst = SQRT(tmp/ndim)

      RETURN
      END                                                    ! of WGTDST
