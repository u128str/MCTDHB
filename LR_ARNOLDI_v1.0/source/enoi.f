C#######################################################################
C
C  Routines related to 1D ENO (Essentially Non Oscillatory)
C  Interpolation.
C  
C  ENO is basically a polynomial interpolation where the range
C  in the spatial coordinate used to perform the interpolation
C  is dynamically selected to avoid discontinuities or huge
C  oscillations. This, hopefully, gives a non oscillatory interpolant.
C
C  OV 02/06
C#######################################################################


C#########################################################################
      subroutine enoi(xgrid,n,diff,xint,rmax,rmin,dfac,yint)

      implicit none
      
c     Arguments
c        xgrid:  input, vector of dimension n storing the coordinate grid
c                       to be interpolated
c        vgrid:  input, vector of dimension n storing the value of the
c                       function, so that vgrid(i)=f(xgrid(i))
c            n:  input, dimension of vectors xgrid and vgrid
c         diff:  input, matrix containing the finite differences
c         xint:  input, x coordinate to be interpolated
c         rmax:  input, max degree in eno interpolation
c         rmin:  input, min degree in eno interpolation
c         dfac:  input, the larger this factor, the more probable
c                       is a symmetric interpolation. Def is 50
c         yint: output, interpolation at point xint, f_int(xint)

      integer  n
      real*8   xgrid(n),diff(n,n),xint,yint,dfac
      
      integer  i,k,kmin,kmax,rmin,rmax
      real*8   poly
      
c Find k, the lowest index of the segment in xgrid where the point 
c to be interpolated, xint, lies
      k=1
      do i=1,n-1
         if( xint .gt. xgrid(i+1) ) k=k+1
      enddo

c Find kmin and kmax, the lowest and highest indexes of the points 
c used to compute the interpolating polynomial. We must supply rmin, 
c and rmax, the minimum and maximum allowed interpolating orders
      call rstencil(diff,n,k,rmin,rmax,dfac,kmin,kmax)


c Interpolate
      yint= poly(xgrid,diff,n,kmin,kmax,xint)


c yint contains now the desired result, exit subroutine
      return
      end
C#########################################################################




      
C#########################################################################
      subroutine rstencil(diff,n,k,rmin,rmax,dfac,kmin,kmax)

c     Computes the maximum and minimum indices of the segment
c     used in the polynomial interpolation between any two
c     points

      implicit none

c     Arguments: see text above, output args are kmin and kmax
      integer n,k,rmin,rmax,kmin,kmax
      real*8  diff(n,n),dfac

      integer  ref,i
      real*8   fl,fr

      kmin=k
      ref=min(k,n-k)
      if(ref .lt. rmin) ref=rmin
      if(ref .gt. rmax) ref=rmax
      
      do i=2,ref
         if(kmin .ne. 1) then
            fl=diff(kmin-1,kmin-1+i)
            fr=diff(kmin,kmin+i)
         
c           give more weight to the symmetrical 3rd order interpolation
c           according to dfac
            if(i .eq. 3) then
                 if(k .eq. kmin) then
                     if( abs(fl) .lt. dfac*abs(fr) ) kmin=kmin-1
                 else
                     if( dfac*abs(fl) .lt. abs(fr) ) kmin=kmin-1
                 endif
            else
                 if( abs(fl) .lt. abs(fr) ) kmin=kmin-1
            endif
         endif
      enddo


      kmax=kmin+ref
c kmin and kmax contain the output results, exit
      return
      end
C#########################################################################
      





C#########################################################################
      real*8 function poly(xgrid,diff,n,kmin,kmax,xint)

c     Given a 1D grid (xgrid), the value of the function at each
c     grid point (vgrid), the indices of the two extreme points
c     used in the polynomial interpolation (kmin,kmax) and a value
c     in the x coordinate, returns the value of the interpolated
c     function. The interpolant is simply a Newton polynomial of 
c     order kmax-kmin going thorough all the points of the segment.

      implicit none

      integer  n,kmin,kmax
      real*8   xgrid(n),diff(n,n),xint

      integer  i,j
      real*8   xroot

      poly=diff(kmin,kmin)
      
      do i=kmin+1,kmax
          xroot=1d0
          do j=kmin,i-1
              xroot=xroot*(xint-xgrid(j))
          enddo
          poly = poly + diff(kmin,i)*xroot
      enddo

      return
      end
C#########################################################################
         

      
      


C#########################################################################
      subroutine fdiff(xgrid,vgrid,n,diff)
      
c     Computes the numerical n-order finite difference between
c     all the points in the given segment. 
      
      implicit none

      integer  n,i,j
      real*8   xgrid(n),vgrid(n),diff(n,n)

      do i=1,n
         do j=i,1,-1
            if(i.eq.j) then
                diff(j,i)=vgrid(i)
            else
                diff(j,i)=(diff(j+1,i)-diff(j,i-1))/(xgrid(i)-xgrid(j))
            endif
         enddo
      enddo
            
      return
      end
C#########################################################################


