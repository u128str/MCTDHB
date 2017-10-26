!C --- linux_suppl.inc          
!C
!C     g77 does not have intrinsic sind, cosd....
!C     They have to be declared external
!C
!      real sind, cosd, tand
!      real asind, acosd, atand, atan2d
!      external sind, cosd, tand
!      external asind, acosd, atand, atan2d
!C
!+++++
!C
!C     linux supplemental fuctions
!C
!C sind
      real function sind (x)
      implicit none
      real x
      real degrad
      parameter (degrad = 3.141592654 / 180.)

      sind = sin (degrad * x)
      return
      end
!C
!C cosd
      real function cosd (x)
      implicit none
      real x
      real degrad
      parameter (degrad = 3.141592654 / 180.)

      cosd = cos (degrad * x)
      return
      end
!C
!C tand
      real function tand (x)
      implicit none
      real x
      real degrad
      parameter (degrad = 3.141592654 / 180.)

      tand = tan (degrad * x)
      return
      end
!C
!C asind
      real function asind (x)
      implicit none
      real x
      real raddeg
      parameter (raddeg =  180. / 3.141592654 )

      asind = raddeg * asin (x)
      return
      end
!C
!C acosd
      real function acosd (x)
      implicit none
      real x
      real raddeg
      parameter (raddeg =  180. / 3.141592654 )

      acosd = raddeg * acos (x)
      return
      end
!C
!C atand
      real function atand (x)
      implicit none
      real x
      real raddeg
      parameter (raddeg =  180. / 3.141592654 )

      atand = raddeg * atan (x)
      return
      end
!C
!C atan2d
      real function atan2d (x, y)
      implicit none
      real x, y
      real raddeg
      parameter (raddeg =  180. / 3.141592654 )

      atan2d = raddeg * atan2 (x, y)
      return
      end
!C Heaviside
!C x<a H(a)=0 
!C x>=a H(a)=1  
      double precision function hvsd(x)
      implicit none
      double precision x
      hvsd = 1.d0
      IF (x<0.d0) hvsd=0.d0
      return
      end
