c================ Factorial
       FUNCTION FKTR(N)
       real*8 FKTR
       integer N,k
        FKTR=1.0d0
        Do k=1,N
        FKTR=FKTR*k
        EndDo
       end
c================ LoGFactorial
       FUNCTION LogFKTR(N)
       real*8 logFKTR
       integer N,k
        logFKTR=0.0d0
        Do k=2,N
        LogFKTR=LogFKTR+ Log(1.0d0*k)
        EndDo
       end
c================ Binomial C(N,K) old
C       FUNCTION CNK(N,K)
C       real*8 CNK,LogFKTR
C       integer N,K
C       external LogFKTR
C       if(N.lt.K) then
C       CNK=0.0d0
C       return 
C       endif
C       CNK=1.0d0
C       CNK= Exp(LogFKTR(N)-LogFKTR(K)-LogFKTR(N-K))
C       write(6,*)"CNK",CNK
Cc       end
c================ Binomial C(N,K)
       FUNCTION CNK(N,K)
       real*8 :: CNK,dGAMMA,beta
       integer :: N,K
c       external ddGAMMA
c       write(6,*)N,K," in CNK"
       if(N.lt.K) then
       CNK=0.0d0
       return 
       endif
cc       CNK=1.0d0
c       write(6,*)" Gamma",DBLE(N+1),Gamma(DBLE(N+1))
c       write(6,*)" Gamma",DBLE(K+1),Gamma(DBLE(K+1))
Cc       write(6,*)" Gamma",DBLE(N-K+1),Gamma(DBLE(N-K+1))
c       CNK= dGamma(DBLE(N+1))/dGamma(DBLE(K+1))/dGamma(DBLE(N-K+1))
c       write(6,*)N,K,"OLD CNK",CNK
c      write(6,*)"Bt",DBLE(N-k+1),DBLE(k+1),1/beta(DBLE(N-K+1),DBLE(K+1))
       CNK= 1/beta(DBLE(N-K+1),DBLE(K+1))/(N+1)
c       write(6,*)N,K,"CNK",CNK
c        pause
       end

c============================ DGAMMA
! Gamma function in double precision
!
      function ddgamma(x)
      implicit real*8 (a - h, o - z)
      real*8 p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13
      real*8  ddgamma,w,y,x
      integer n,k
      parameter (
     &    p0 = 0.999999999999999990d+00, 
     &    p1 = -0.422784335098466784d+00, 
     &    p2 = -0.233093736421782878d+00, 
     &    p3 = 0.191091101387638410d+00, 
     &    p4 = -0.024552490005641278d+00, 
     &    p5 = -0.017645244547851414d+00, 
     &    p6 = 0.008023273027855346d+00)
      parameter (
     &    p7 = -0.000804329819255744d+00, 
     &    p8 = -0.000360837876648255d+00, 
     &    p9 = 0.000145596568617526d+00, 
     &    p10 = -0.000017545539395205d+00, 
     &    p11 = -0.000002591225267689d+00, 
     &    p12 = 0.000001337767384067d+00, 
     &    p13 = -0.000000199542863674d+00)
      n = nint(x - 2)
      w = x - (n + 2)
      y = ((((((((((((p13 * w + p12) * w + p11) * w + p10) * 
     &    w + p9) * w + p8) * w + p7) * w + p6) * w + p5) * 
     &    w + p4) * w + p3) * w + p2) * w + p1) * w + p0
      if (n .gt. 0) then
          w = x - 1
          do k = 2, n
              w = w * (x - k)
          end do
      else
          w = 1
          do k = 0, -n - 1
              y = y * (x + k)
          end do
      end if
      ddgamma = w / y
c       write(6,*)" Gamma",x,ddgamma
      end

      function dlgama ( x )
cc DLGAMA evaluates log ( Gamma ( X ) ) for a real argument.
c
c  Discussion:
c
c    This routine calculates the LOG(GAMMA) function for a positive real
c    argument X.  Computation is based on an algorithm outlined in
c    references 1 and 2.  The program uses rational functions that
c    theoretically approximate LOG(GAMMA) to at least 18 significant
c    decimal digits.  The approximation for X > 12 is from reference
c    3, while approximations for X < 12.0 are similar to those in
c    reference 1, but are unpublished.
c
c  Modified:
c
c    03 April 2007
c
c  Author:
c
c    William Cody, Laura Stoltz
c
c  Reference:
c
c    William Cody, Kenneth Hillstrom,
c    Chebyshev Approximations for the Natural Logarithm of the 
c    Gamma Function,
c    Mathematics of Computation,
c    Volume 21, Number 98, April 1967, pages 198-203.
c
c    Kenneth Hillstrom,
c    ANL/AMD Program ANLC366S, DGAMMA/DLGAMA,
c    May 1969.
c
c    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, 
c    Charles Mesztenyi, John Rice, Henry Thatcher, 
c    Christoph Witzgall,
c    Computer Approximations,
c    Wiley, 1968,
c    LC: QA297.C64.
c
c  Parameters:
c
c    Input, double precision X, the argument of the function.
c
c    Output, double precision DLGAMA, the value of the function.
c
      implicit none

      double precision c(7)
      double precision corr
      double precision d1
      double precision d2
      double precision d4
      double precision dlgama
      double precision eps
      double precision frtbig
      double precision four
      double precision half
      integer i
      double precision one
      double precision pnt68
      double precision p1(8)
      double precision p2(8)
      double precision p4(8)
      double precision q1(8)
      double precision q2(8)
      double precision q4(8)
      double precision res
      double precision sqrtpi
      double precision thrhal
      double precision twelve
      double precision two
      double precision x
      double precision xbig
      double precision xden
      double precision xinf
      double precision xm1
      double precision xm2
      double precision xm4
      double precision xnum
      double precision y
      double precision ysq
      double precision zero
c
c  Mathematical constants
c
      data one /1.0D+00/
      data half /0.5D+00/
      data twelve /12.0D+00/
      data zero /0.0D+00/
      data four /4.0D+00/
      data thrhal /1.5D+00/
      data two /2.0D+00/
      data pnt68 /0.6796875D+00/
      data sqrtpi /0.9189385332046727417803297D+00/
c
c  Machine dependent parameters
c
      data xbig /2.55D+305/
      data xinf /1.79D+308/
      data eps /2.22D-16/
      data frtbig /2.25D+76/
c
c  Numerator and denominator coefficients for rational minimax
c  approximation over (0.5,1.5).
c
      data d1/-5.772156649015328605195174D-01/
      data p1/
     &   4.945235359296727046734888D+00,
     &   2.018112620856775083915565D+02,
     &   2.290838373831346393026739D+03,
     &   1.131967205903380828685045D+04,
     &   2.855724635671635335736389D+04,
     &   3.848496228443793359990269D+04,
     &   2.637748787624195437963534D+04,
     &   7.225813979700288197698961D+03/
      data q1/
     &   6.748212550303777196073036D+01,
     &   1.113332393857199323513008D+03,
     &   7.738757056935398733233834D+03,
     &   2.763987074403340708898585D+04,
     &   5.499310206226157329794414D+04,
     &   6.161122180066002127833352D+04,
     &   3.635127591501940507276287D+04,
     &   8.785536302431013170870835D+03/
c
c  Numerator and denominator coefficients for rational minimax
c  Approximation over (1.5,4.0).
c
      data d2/4.227843350984671393993777D-01/
      data p2/
     &   4.974607845568932035012064D+00,
     &   5.424138599891070494101986D+02,
     &   1.550693864978364947665077D+04,
     &   1.847932904445632425417223D+05,
     &   1.088204769468828767498470D+06,
     &   3.338152967987029735917223D+06,
     &   5.106661678927352456275255D+06,
     &   3.074109054850539556250927D+06/
      data q2/
     &   1.830328399370592604055942D+02,
     &   7.765049321445005871323047D+03,
     &   1.331903827966074194402448D+05,
     &   1.136705821321969608938755D+06,
     &   5.267964117437946917577538D+06,
     &   1.346701454311101692290052D+07,
     &   1.782736530353274213975932D+07,
     &   9.533095591844353613395747D+06/
c
c  Numerator and denominator coefficients for rational minimax
c  Approximation over (4.0,12.0).
c
      data d4/1.791759469228055000094023D+00/
      data p4/
     &   1.474502166059939948905062D+04,
     &   2.426813369486704502836312D+06,
     &   1.214755574045093227939592D+08,
     &   2.663432449630976949898078D+09,
     &   2.940378956634553899906876D+10,
     &   1.702665737765398868392998D+11,
     &   4.926125793377430887588120D+11,
     &   5.606251856223951465078242D+11/
      data q4/
     &   2.690530175870899333379843D+03,
     &   6.393885654300092398984238D+05,
     &   4.135599930241388052042842D+07,
     &   1.120872109616147941376570D+09,
     &   1.488613728678813811542398D+10,
     &   1.016803586272438228077304D+11,
     &   3.417476345507377132798597D+11,
     &   4.463158187419713286462081D+11/
c
c  Coefficients for minimax approximation over (12, INF).
c
      data c/
     &  -1.910444077728D-03,
     &   8.4171387781295D-04,
     &  -5.952379913043012D-04,
     &   7.93650793500350248D-04,
     &  -2.777777777777681622553D-03,
     &   8.333333333333333331554247D-02,
     &   5.7083835261D-03/

      y = x

      if ( zero .lt. y .and. y .le. xbig ) then

        if ( y .le. eps ) then

          res = - dlog ( y )
c
c  EPS < X <= 1.5.
c
        else if ( y .le. thrhal ) then

          if ( y .lt. pnt68 ) then
            corr = - dlog ( y )
            xm1 = y
          else
            corr = zero
            xm1 = ( y - half ) - half
          end if

          if ( y .le. half .or. pnt68 .le. y ) then

            xden = one
            xnum = zero
            do i = 1, 8
              xnum = xnum * xm1 + p1(i)
              xden = xden * xm1 + q1(i)
            end do

            res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) )

          else

            xm2 = ( y - half ) - half
            xden = one
            xnum = zero
            do i = 1, 8
              xnum = xnum * xm2 + p2(i)
              xden = xden * xm2 + q2(i)
            end do

            res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) )

          end if
c
c  1.5 < X <= 4.0.
c
        else if ( y .le. four ) then

          xm2 = y - two
          xden = one
          xnum = zero
          do i = 1, 8
            xnum = xnum * xm2 + p2(i)
            xden = xden * xm2 + q2(i)
          end do

          res = xm2 * ( d2 + xm2 * ( xnum / xden ) )
c
c  4.0 < X <= 12.0.
c
        else if ( y .le. twelve ) then

          xm4 = y - four
          xden = -one
          xnum = zero
          do i = 1, 8
            xnum = xnum * xm4 + p4(i)
            xden = xden * xm4 + q4(i)
          end do

          res = d4 + xm4 * ( xnum / xden )
c
c  Evaluate for 12 <= argument.
c
        else

          res = zero

          if ( y .le. frtbig ) then

            res = c(7)
            ysq = y * y

            do i = 1, 6
              res = res / ysq + c(i)
            end do

          end if

          res = res / y
          corr = dlog ( y )
          res = res + sqrtpi - half * corr
          res = res + y * ( corr - one )

        end if
c
c  Return for bad arguments.
c
      else

        res = xinf

      end if
c
c  Final adjustments and return.
c
      dlgama = res

      return
      end

       function beta ( a, b )
c*********************************************************************72
c
cc BETA returns the value of the Beta function.
c
c  Formula:
c
c    BETA(A,B) = ( GAMMA ( A ) * GAMMA ( B ) ) / GAMMA ( A + B )
c              = Integral ( 0 <= T <= 1 ) T**(A-1) (1-T)**(B-1) dT.
c
c  Modified:
c
c    10 July 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision A, B, the parameters of the function.
c    0.0 < A,
c    0.0 < B.
c
c    Output, double precision BETA, the value of the function.
c
      implicit none

      double precision a
      double precision b
      double precision beta
      double precision dlgama
c
c  Check.
c
      if ( a <= 0.0D+00 .or. b <= 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BETA - Fatal error!'
        write ( 6, * ) 'A,B',A,B
        write ( *, '(a)' ) '  Both A and B must be greater than 0.'
        stop
      end if

      beta = exp ( dlgama ( a ) 
     &           + dlgama ( b ) 
     &           - dlgama ( a + b ) )

      return
      end
       
