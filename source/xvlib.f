********************************************************************
C                                                                     
C                              XVLIB                                  
C                                                                     
C  Library module containing linear algebra routines that involve the 
C  multiplication of vectors or scalars.
C                                                                     
C  Nomenclature:
C    Each name has 6 basic characters:
C    First 2 characters denote the objects being multiplied:
C       q: quadratic matrix
C       m: general (rectangular) matrix
C       h: hermitian matrix
C       p: positive definite matrix
C       s: symmetric matrix
C       u: unitary matrix
C       d: diagonal matrix (only diagonal elements are supplied as a
C          vector)
C       t: tensor of third order
C       v: vector
C       x: scalar
C       e.g. 'xv' denotes the operation (scalar * vector)
C    Character 3 denotes how first object is used:
C       x: unchanged from input
C       t: transpose of input
C       a: adjoint of input
C       c: complex conjugate of input
C       i: inverse (real and complex scalar only)
C    Character 4 denotes how second object is used:
C       see Character 3 above.
C    Character 5, 6 denote data types of first, second object
C       respectively:
C       z: complex double precision (complex*16)
C       c: complex single precision (complex*8)
C       d: real double precision (real*8)
C       r: real single precision (real*4)
C       y: complex matrix stored as two double precision (real*8)
C          matrices
C     Further characters, if present, give more informaion:
C       a: the result is added to a further object
C       r: the result is subtracted (removed) from a further object
C       o: the result is stored (overwritten) in the second input object
C       c: the input matrices commute
C       h: the resulting matrix is hermitian
C       s: the resulting matrix is symmetric
C       1: the physical dimensions of the matrices differs from those
C          used.
C   
CC  xqxxzz (x,a,c,dim)
CC      Definition: x*a(j,i) = c(j,i)
CC      Dimensions: a(dim,dim),c(dim,dim),x
C
CC  xqxxzza (x,a,c,dim)
CC      Definition: x*a(j,i) + c(j,i) = c(j,i)
CC      Dimensions: a(dim,dim),c(dim,dim),x
C
CC  xmxxzz (x,a,c,dim1,dim2)
CC      Definition: x*a(j,i) = c(j,i)
CC      Dimensions: a(dim1,dim2),c(dim1,dim2),x
C
C   xmxxzza (x,a,c,dim1,dim2)
C       Definition: x*a(j,i) + c(j,i) = c(j,i)
C       Dimensions: a(dim1,dim2),c(dim1,dim2),x
C
C   xmxxdza (x,a,c,dim1,dim2)
C       Definition: x*a(j,i) + c(j,i) = c(j,i)
C       Dimensions: a(dim1,dim2),c(dim1,dim2),x
C
C   xmxtzza (x,a,c,dim1,dim2)
C       Definition: x*a(i,j) + c(j,i) = c(j,i)
C       Dimensions: a(dim2,dim1),c(dim1,dim2),x
C
C   xqxxdd (x,a,c,dim)
C       Definition: x*a(j,i) = c(j,i)
C       Dimensions: a(dim,dim),c(dim,dim),x
C
C   xqxxdda (x,a,c,dim)
C       Definition: x*a(j,i) + c(j,i) = c(j,i)
C       Dimensions: a(dim,dim),c(dim,dim),x
C
C   xqxxzza (x,a,c,dim)
C       Definition: x*a(j,i) + c(j,i) = c(j,i)
C       Dimensions: a(dim,dim),c(dim,dim),x
C
C   xqxxdda1 (x,a,c,phdim,dim)
C       Definition: x*a(j,i) + c(j,i) = c(j,i) ; 1 <= i,j <= dim
C       Dimensions: a(phdim,dim),c(phdim,dim),x
C
C   xqxxzd (x,a,c,dim)
C       Definition: x*a(j,i) = c(j,i)
C       Dimensions: a(dim,dim),c(dim,dim),x
C
C   xdxxdda (x,a,c,dim)
C       Definition: x*a(i) + c(i,i) = c(i,i)
C       Dimensions: a(dim),c(dim,dim),x
C
C   xdxsdda (x,a,c,dim)
C       Definition: x*a(i)**2 + c(i,i) = c(i,i)
C       Dimensions: a(dim),c(dim,dim),x
C
C   xuxxdda (x,c,dim)
C       Definition: x + c(i,i) = c(i,i)
C       Dimensions: c(dim,dim),x
C
C   xvixddo (x,v,dim)
C       Definition: inv(x)*v(i) = v(i)
C       Dimensions: v(dim),x
C
C   xvixzzo (x,v,dim)
C       Definition: inv(x)*v(i) = v(i)
C       Dimensions: v(dim),x
C
C   xvixdzo (x,v,dim)
C       Definition: inv(x)*v(i) = v(i)
C       Dimensions: v(dim),x
C
C   xvxxzza (x,v,w,dim)
C       Definition: x*v(i) + w(i) = w(i)
C       Dimensions: v(dim),w(dim),x
C
C   xvxxdzr (x,v,w,dim)
C       Definition: w(i) - x*v(i) = w(i)
C       Dimensions: v(dim),w(dim),x
C
C   xvxxdza (x,v,w,dim)
C       Definition: w(i) + x*v(i) = w(i)
C       Dimensions: v(dim),w(dim),x
C
C   xvxxzzo (x,v,dim)
C       Definition: x*v(i) = v(i)
C       Dimensions: v(dim),x
C
C   xvxxdzo (x,v,dim)
C       Definition: x*v(i) = v(i)
C       Dimensions: v(dim),x
C
C   xvxxddo (x,v,dim)
C       Definition: x*v(i) = v(i)
C       Dimensions: v(dim),x
C
C   xvxxdds (x,v,w,dim)
C       Definition: w(i) - x*v(i) = w(i)
C       Dimensions: v(dim),w(dim),x
C
C   xvxxzzs (x,v,w,dim)
C       Definition: w(i) - x*v(i) = w(i)
C       Dimensions: v(dim),w(dim),x
C
C   xqxxzzo (x,a,dim)
C       Definition: x*a(j,i) = a(j,i)
C       Dimensions: a(dim,dim),x
C
C   xqxxdzo (x,a,dim)
C       Definition: x*a(j,i) = a(j,i)
C       Dimensions: a(dim,dim),x
C
C   xqixdzo (x,a,dim)
C       Definition: inv(x)*a(j,i) = a(j,i)
C       Dimensions: a(dim,dim),x
C
C   xmxxzzo (x,a,dim1,dim2)
C       Definition: x*a(j,i) = a(j,i)
C       Dimensions: a(dim1,dim2),x
C
C   vvtxdd (u,v,x,dim)
C       Definition: u(i)*v(i) = x
C       Dimensions: u(dim),v(dim),x
C
C   vvtxzz (u,v,x,dim)
C       Definition: u(i)*v(i) = x
C       Dimensions: u(dim),v(dim),x
C
C   vvtxzza (u,v,x,dim)
C       Definition: u(i)*v(i) = x
C       Dimensions: u(dim),v(dim),x
C
C   vvaxzz (u,v,x,dim)
C       Definition: dconjg(u(i))*v(i) = x
C       Dimensions: u(dim),v(dim),x
C
C   vvxazz (u,v,m,dim)
C       Definition: u(i)*dconjg(v(j)) = m(i,j)
C       Dimensions: d(dim),v(dim),m(dim,dim)
C
C   qvxxdd (a,v,w,dim)
C       Definition: a(i,j)*v(j) = w(i) .
C       Dimensions: a(dim,dim),v(dim),w(dim)
C
C   qvxxzz (a,v,w,dim)
C       Definition: a(i,j)*v(j) = w(i) .
C       Dimensions: a(dim,dim),v(dim),w(dim)
C
C   qvxxdd1 (a,v,w,phdim,dim)
C       Definition: a(i,j)*v(j) = w(i) ; 1 <= i,j <= dim .
C       Dimensions: a(phdim,dim),v(dim),w(dim)
C
C   qvtxdd (a,v,w,dim)
C       Definition: a(j,i)*v(j) = w(i) .
C       Dimensions: a(dim,dim),v(dim),w(dim)
C
C   qvtxdd1 (a,v,w,phdim,dim)
C       Definition: a(j,i)*v(j) = w(i) ; 1 <= i,j <= dim .
C       Dimensions: a(phdim,dim),v(dim),w(dim)
C
C   mvtxdd (a,v,w,dim1,dim2)
C       Definition: a(j,i)*v(j) = w(i)
C       Dimensions: a(dim1,dim2),v(dim1),w(dim2)
C
C   mvtxdd1 (a,v,w,phdim,dim1,dim2)
C       Definition: a(j,i)*v(j) = w(i)
C       Dimensions: a(dim1,dim2),v(dim1),w(dim2) (= used dimensions)
C       Dimensions: a(phdim,phdim) (= allocated dimensions)
C
C   mvtxzz (a,v,w,dim1,dim2)
C       Definition: a(j,i)*v(j) = w(i) .
C       Dimensions: a(dim1,dim2),v(dim1),w(dim2)
C
C   mvxxzz (a,v,w,dim1,dim2)
C       Definition: a(i,j)*v(j) = w(i) .
C       Dimensions: a(dim1,dim2),v(dim2),w(dim1)
C
C   dvxxdd (a,v,w,dim)
C       Definition: a(i)*v(i) = w(i) .
C       Dimensions: a(dim),v(dim),w(dim)
C
C   dvxxdz (a,v,w,dim)
C       Definition: a(i)*v(i) = w(i) .
C       Dimensions: a(dim),v(dim),w(dim)
C
C   dvxxdzo (a,v,dim)
C       Definition: a(i)*v(i) = v(i) .
C       Dimensions: a(dim),v(dim)
C
C   dvxxddo (a,v,dim)
C       Definition: a(i)*v(i) = v(i) .
C       Dimensions: a(dim),v(dim)
C
C   qvxxdz (a,v,w,dim)
C       Definition: a(i,j)*v(j) = w(i) .
C       Dimensions: a(dim,dim),v(dim),w(dim)
C
C   qvxxdz1 (a,v,w,phdim,dim)
C       Definition: a(i,j)*v(j) = w(i) ; 1 <= i,j <= dim .
C       Dimensions: a(phdim,dim),v(dim),w(dim)
C
C   qvxxzz1 (a,v,w,phdim,dim)
C       Definition: a(i,j)*v(j) = w(i) ; 1 <= i,j <= dim .
C       Dimensions: a(phdim,dim),v(dim),w(dim)
C
C   qvxxyz (a,v,w,dim)
C       Definition: a(i,j)*v(j) = w(i) .
C       Dimensions: a(dim,dim),v(dim),w(dim)
C
C ----------------------------------------------------------------------


C-----------------------------------------------------------------------
C Library subroutine xqxxzz
C
C Multiplication of a complex scalar with a complex quadratic matrix:
C     x*a(j,i)=c(j,i)
C-----------------------------------------------------------------------

C     subroutine xqxxzz (x,a,c,dim)

C     implicit none

C     integer dim,i,j
C     complex*16  x,a(dim,dim),c(dim,dim)

C     do i = 1,dim
C        do j = 1,dim
C           c(j,i) = x*a(j,i)
C        enddo
C     enddo

C     return
C     end

C ----------------------------------------------------------------------
C Library subroutine xqxxzza
C
C Multiplication of a complex scalar with a complex quadratic matrix,
C the result of which is added to a further input complex quadratic
C matrix:
C     x*a(j,i) + c(j,i) = c(j,i)
C
C NB Input c matrix is overwritten on output
C-----------------------------------------------------------------------

C     subroutine xqxxzza (x,a,c,dim)

C     implicit none

C     integer dim,i,j
C     complex*16  x,a(dim,dim),c(dim,dim)

C     do i = 1,dim
C        do j = 1,dim
C           c(j,i) = c(j,i)+x*a(j,i)
C        enddo
C     enddo

C     return
C     end

C ----------------------------------------------------------------------
C Library subroutine xmxxzz
C
C Multiplication of a complex scalar with a complex rectangular matrix:
C     x*a(j,i)=c(j,i)
C-----------------------------------------------------------------------

C     subroutine xmxxzz (x,a,c,dim1,dim2)

C     implicit none

C     integer dim1,dim2,i,j
C     complex*16  x,a(dim1,dim2),c(dim1,dim2)

C     do i = 1,dim2
C        do j = 1,dim1
C           c(j,i) = x*a(j,i)
C        enddo
C     enddo

C     return
C     end

C ----------------------------------------------------------------------
C Library subroutine xmxxzza
C
C Multiplication of a complex scalar with a complex rectangular matrix,
C the result of which is added to a further input complex rectangular
C matrix:
C     x*a(j,i) + c(j,i) = c(j,i)
C
C NB Input c matrix is overwritten on output
C-----------------------------------------------------------------------

      subroutine xmxxzza (x,a,c,dim1,dim2)

      implicit none

      integer dim1,dim2,i,j
      complex*16  x,a(dim1,dim2),c(dim1,dim2)

      do i = 1,dim2
         do j = 1,dim1
            c(j,i) = c(j,i)+x*a(j,i)
         enddo
      enddo

      return
      end

C ----------------------------------------------------------------------
C Library subroutine xmxxdza
C
C Multiplication of a real scalar with a complex rectangular matrix,
C the result of which is added to a further input complex rectangular
C matrix:
C     x*a(j,i) + c(j,i) = c(j,i)
C
C NB Input c matrix is overwritten on output
C-----------------------------------------------------------------------

      subroutine xmxxdza (x,a,c,dim1,dim2)

      implicit none

      integer dim1,dim2,i,j
      real*8  x
      complex*16  a(dim1,dim2),c(dim1,dim2)

      do i = 1,dim2
         do j = 1,dim1
            c(j,i) = c(j,i)+x*a(j,i)
         enddo
      enddo

      return
      end

C ----------------------------------------------------------------------
C Library subroutine xmxtzza
C
C Multiplication of a complex scalar with a transposed complex
C rectangular matrix, the result of which is added to a further input
C complex rectangular matrix:
C     x*a(i,j) + c(j,i) = c(j,i)
C
C NB Input c matrix is overwritten on output
C-----------------------------------------------------------------------

      subroutine xmxtzza (x,a,c,dim1,dim2)

      implicit none

      integer dim1,dim2,i,j
      complex*16  x,a(dim2,dim1),c(dim1,dim2)

      do i = 1,dim2
         do j = 1,dim1
            c(j,i) = c(j,i)+x*a(i,j)
         enddo
      enddo

      return
      end

C ----------------------------------------------------------------------
C Library subroutine xqxxdd
C
C Multiplication of a real scalar with a real quadratic matrix:
C     x*a(j,i)=c(j,i)
C-----------------------------------------------------------------------

      subroutine xqxxdd (x,a,c,dim)

      implicit none

      integer dim,i,j
      real*8  x,a(dim,dim),c(dim,dim)

      do i = 1,dim
         do j = 1,dim
            c(j,i) = x*a(j,i)
         enddo
      enddo

      return
      end

C ----------------------------------------------------------------------
C Library subroutine xqxxdda
C
C Multiplication of a real scalar with a real quadratic matrix,
C the result of which is added to a further input real quadratic matrix:
C     x*a(j,i) + c(j,i) = c(j,i)
C
C NB Input c matrix is overwritten on output
C-----------------------------------------------------------------------

      subroutine xqxxdda (x,a,c,dim)

      implicit none

      integer dim,i,j
      real*8  x,a(dim,dim),c(dim,dim)

      do i = 1,dim
         do j = 1,dim
            c(j,i) = c(j,i)+x*a(j,i)
         enddo
      enddo

      return
      end

C ----------------------------------------------------------------------
C Library subroutine xqxxzza
C
C Multiplication of a complex scalar with a complex quadratic matrix,
C the result of which is added to a further input complex quadratic
C matrix:
C     x*a(j,i) + c(j,i) = c(j,i)
C
C NB Input c matrix is overwritten on output
C-----------------------------------------------------------------------

      subroutine xqxxzza (x,a,c,dim)

      implicit none

      integer    dim,i,j
      complex*16 x,a(dim,dim),c(dim,dim)

      do i = 1,dim
         do j = 1,dim
            c(j,i) = c(j,i)+x*a(j,i)
         enddo
      enddo

      return
      end

C ----------------------------------------------------------------------
C Library subroutine xqxxdda1
C
C Multiplication of a real scalar with a real quadratic matrix,
C the result of which is added to a further input real quadratic matrix:
C     x*a(j,i) + c(j,i) = c(j,i)
C
C NB Input c matrix is overwritten on output
C-----------------------------------------------------------------------

      subroutine xqxxdda1 (x,a,c,phdim,dim)

      implicit none

      integer phdim,dim,i,j
      real*8  x,a(phdim,dim),c(phdim,dim)

      do i = 1,dim
         do j = 1,dim
            c(j,i) = c(j,i)+x*a(j,i)
         enddo
      enddo

      return
      end

C ----------------------------------------------------------------------
C Library subroutine xqxxzd
C
C Multiplication of a complex scalar with a real quadratic matrix:
C     x*a(j,i)=c(j,i)
C-----------------------------------------------------------------------

      subroutine xqxxzd (x,a,c,dim)

      implicit none

      integer dim,i,j
      real*8  a(dim,dim)
      complex*16  x,c(dim,dim)

      do i = 1,dim
         do j = 1,dim
            c(j,i) = x*a(j,i)
         enddo
      enddo

      return
      end

C ----------------------------------------------------------------------
C Library subroutine xdxxdda
C
C Multiplication of a real scalar with a real diagonal matrix, where
C the result is added to a complex quadratic matrix
C     x*a(i,i) + c(i,i) = c(i,i)
C-----------------------------------------------------------------------

      subroutine xdxxdda (x,a,c,dim)

      implicit none

      integer dim,i
      real*8  x, a(dim)
      complex*16  c(dim,dim)

      do i = 1,dim
         c(i,i) = c(i,i) + x*a(i)
      enddo

      return
      end

C ----------------------------------------------------------------------
C Library subroutine xdxsdda
C
C Multiplication of a real scalar with a squared real diagonal matrix, 
C where the result is added to a complex quadratic matrix
C     x*a(i,i)**2 + c(i,i) = c(i,i)
C-----------------------------------------------------------------------

      subroutine xdxsdda (x,a,c,dim)

      implicit none

      integer dim,i
      real*8  x, a(dim)
      complex*16  c(dim,dim)

      do i = 1,dim
         c(i,i) = c(i,i) + x*a(i)**2
      enddo

      return
      end

C ----------------------------------------------------------------------
C Library subroutine xuxxdda
C
C Multiplication of a real scalar with the unity matrix, where
C the result is added to a complex quadratic matrix
C     x + c(i,i) = c(i,i)
C-----------------------------------------------------------------------

      subroutine xuxxdda (x,c,dim)

      implicit none

      integer dim,i
      real*8  x
      complex*16  c(dim,dim)

      do i = 1,dim
         c(i,i) = c(i,i) + x
      enddo

      return
      end

C ----------------------------------------------------------------------
C Library subroutine xvixddo
C
C Multiplication of the inverse of a real scalar with a real vector,
C where the result is stored in the initial array:
C     inv(x)*v(i) = v(i)
C
C NB Input v vector is overwritten on output
C    This routine can be used for the normalisation of a vector if
C    x is the sqrt of the vector scalar product.
C-----------------------------------------------------------------------

      subroutine xvixddo (x,v,dim)

      implicit none

      integer dim,i
      real*8  v(dim),x,inv

      inv=1.0d0/x
      do i = 1,dim
         v(i) = inv*v(i)
      enddo

      return
      end

C ----------------------------------------------------------------------
C Library subroutine xvixzzo
C
C Multiplication of the inverse of a complex scalar with a complex
C vector, where the result is stored in the initial array:
C     inv(x)*v(i) = v(i)
C
C NB Input v vector is overwritten on output
C    This routine can be used for the normalisation of a vector if
C    x is the sqrt of the vector's SYMMETRIC scalar product.
C-----------------------------------------------------------------------

      subroutine xvixzzo (x,v,dim)

      implicit none

      integer    dim,i
      complex*16 v(dim),x,inv

      inv=(1.0d0,0.0d0)/x
      do i = 1,dim
         v(i) = inv*v(i)
      enddo

      return
      end

C ----------------------------------------------------------------------
C Library subroutine xvixdzo
C
C Multiplication of the inverse of a real scalar with a complex vector,
C where the result is stored in the initial array:
C     inv(x)*v(i) = v(i)
C
C NB Input v vector is overwritten on output
C    This routine can be used for the normalisation of a vector if
C    x is the sqrt of the vector scalar product.
C-----------------------------------------------------------------------

      subroutine xvixdzo (x,v,dim)

      implicit none

      integer dim,i
      real*8  x,inv
      complex*16 v(dim)

      inv=1.0d0/x
      do i = 1,dim
         v(i) = inv*v(i)
      enddo

      return
      end

C ----------------------------------------------------------------------
C Library subroutine xvxxzza
C
C Multiplication of a complex scalar with a complex vector,
C where the result is added to a further vector
C     x*v(i) +w(i) = w(i)
C
C NB Input w vector is overwritten on output
C-----------------------------------------------------------------------

      subroutine xvxxzza (x,v,w,dim)

      implicit none

      integer dim,i
      complex*16 v(dim),w(dim),x

      do i = 1,dim
         w(i) = w(i)+x*v(i)
      enddo

      return
      end

C ----------------------------------------------------------------------
C Library subroutine xvxxdzr
C
C Multiplication of a real scalar with a complex vector,
C where the result is subtracted from a further vector
C     w(i) - x*v(i) = w(i)
C
C NB Input w vector is overwritten on output
C-----------------------------------------------------------------------

      subroutine xvxxdzr (x,v,w,dim)

      implicit none

      integer    dim,i
      real*8     x
      complex*16 v(dim),w(dim)

      do i = 1,dim
         w(i) = w(i)-x*v(i)
      enddo

      return
      end

C ----------------------------------------------------------------------
C Library subroutine xvxxdza
C
C Multiplication of a real scalar with a complex vector,
C where the result is added to a further vector
C     w(i) + x*v(i) = w(i)
C
C NB Input w vector is overwritten on output
C-----------------------------------------------------------------------

      subroutine xvxxdza (x,v,w,dim)

      implicit none

      integer    dim,i
      real*8     x
      complex*16 v(dim),w(dim)

      do i = 1,dim
         w(i) = w(i)+x*v(i)
      enddo

      return
      end

C ----------------------------------------------------------------------
C Library subroutine xvxxdds
C
C Multiplication of a real scalar with a real vector,
C where the result is subtracted from a further vector
C     w(i) - x*v(i) = w(i)
C
C NB Input w vector is overwritten on output
C-----------------------------------------------------------------------

      subroutine xvxxdds (x,v,w,dim)

      implicit none

      integer dim,i
      real*8  v(dim),w(dim),x

      do i = 1,dim
         w(i) = w(i)-x*v(i)
      enddo

      return
      end

C ----------------------------------------------------------------------
C Library subroutine xvxxzzo
C
C Multiplication of a complex scalar with a complex vector,
C where the result is stored in the input vector
C     x*v(i) = v(i)
C
C NB Input v vector is overwritten on output
C-----------------------------------------------------------------------

      subroutine xvxxzzo (x,v,dim)

      implicit none

      integer dim,i
      complex*16 v(dim),x

      do i = 1,dim
         v(i) = x*v(i)
      enddo

      return
      end

C ----------------------------------------------------------------------
C Library subroutine xvxxdzo
C
C Multiplication of a real scalar with a complex vector,
C where the result is stored in the input vector
C     x*v(i) = v(i)
C
C NB Input v vector is overwritten on output
C-----------------------------------------------------------------------

      subroutine xvxxdzo (x,v,dim)

      implicit none

      integer dim,i
      real*8  x
      complex*16 v(dim)

      do i = 1,dim
         v(i) = x*v(i)
      enddo

      return
      end

C ----------------------------------------------------------------------
C Library subroutine xvxxddo
C
C Multiplication of a real scalar with a real vector,
C where the result is stored in the input vector
C     x*v(i) = v(i)
C
C NB Input v vector is overwritten on output
C-----------------------------------------------------------------------

      subroutine xvxxddo (x,v,dim)

      implicit none

      integer dim,i
      real*8  v(dim),x

      do i = 1,dim
         v(i) = x*v(i)
      enddo

      return
      end

C ----------------------------------------------------------------------
C Library subroutine xvxxzzs
C
C Multiplication of a complex scalar with a complex vector,
C where the result is subtracted from a further vector
C     w(i) - x*v(i) = w(i)
C
C NB Input w vector is overwritten on output
C-----------------------------------------------------------------------

      subroutine xvxxzzs (x,v,w,dim)

      implicit none

      integer dim,i
      complex*16 v(dim),w(dim),x

      do i = 1,dim
         w(i) = w(i)-x*v(i)
      enddo

      return
      end

C ----------------------------------------------------------------------
C Library subroutine xqxxzzo
C
C Multiplication of a complex scalar with a quadratic complex matrix,
C where the result is stored in the input matrix
C     x*a(j,i) = a(j,i)
C
C NB Input a matrix is overwritten on output
C-----------------------------------------------------------------------

      subroutine xqxxzzo (x,a,dim)

      implicit none

      integer dim,i,j
      complex*16 a(dim,dim),x

      do i = 1,dim
         do j = 1,dim
            a(j,i) = x*a(j,i)
         enddo
      enddo

      return
      end

C ----------------------------------------------------------------------
C Library subroutine xqxxdzo
C
C Multiplication of a real scalar with a quadratic complex matrix,
C where the result is stored in the input matrix
C     x*a(j,i) = a(j,i)
C
C NB Input a matrix is overwritten on output
C-----------------------------------------------------------------------

      subroutine xqxxdzo (x,a,dim)

      implicit none

      integer dim,i,j
      real*8 x
      complex*16 a(dim,dim)

      do i = 1,dim
         do j = 1,dim
            a(j,i) = x*a(j,i)
         enddo
      enddo

      return
      end

C ----------------------------------------------------------------------
C Library subroutine xqixdzo
C
C Multiplication of the inverse of a real scalar with a quadratic 
C complex matrix, where the result is stored in the input matrix
C     x*a(j,i) = a(j,i)
C
C NB Input a matrix is overwritten on output
C-----------------------------------------------------------------------

      subroutine xqixdzo (x,a,dim)

      implicit none

      integer dim,i,j
      real*8 x,inv
      complex*16 a(dim,dim)

      inv=1.0d0/x
      do i = 1,dim
         do j = 1,dim
            a(j,i) = inv*a(j,i)
         enddo
      enddo

      return
      end

C ----------------------------------------------------------------------
C Library subroutine xmxxzzo
C
C Multiplication of a complex scalar with a rectangular complex matrix,
C where the result is stored in the input matrix
C     x*a(j,i) = a(j,i)
C
C NB Input a matrix is overwritten on output
C-----------------------------------------------------------------------

      subroutine xmxxzzo (x,a,dim1,dim2)

      implicit none

      integer dim1,dim2,i,j
      complex*16 a(dim1,dim2),x

      do i = 1,dim2
         do j = 1,dim1
            a(j,i) = x*a(j,i)
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine vvtxdd
C
C scalar product of two real vectors:
C     u(i)*v(i)=s
C-----------------------------------------------------------------------

      subroutine vvtxdd (u,v,s,dim)

      implicit none

      integer dim,i
      real*8  u(dim),v(dim),s

      s = u(1)*v(1)
      do i = 2,dim
         s = s+u(i)*v(i)
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine vvtxzz
C
C symmetric product of two complex vectors:
C     u(i)*v(i)=s
C-----------------------------------------------------------------------

      subroutine vvtxzz (u,v,s,dim)

      implicit none

      integer    dim,i
      complex*16 u(dim),v(dim),s

      s = u(1)*v(1)
      do i = 2,dim
         s = s+u(i)*v(i)
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine vvtxzz
C
C symmetric product of two complex vectors (result is added to s):
C     u(i)*v(i)=s
C-----------------------------------------------------------------------

      subroutine vvtxzza (u,v,s,dim)

      implicit none

      integer    dim,i
      complex*16 u(dim),v(dim),s

      do i = 1,dim
         s = s+u(i)*v(i)
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine vvaxzz
C
C scalar product of two complex vectors:
C     dconjg(u(i))*v(i)=s
C-----------------------------------------------------------------------
      subroutine vvaxzz (u,v,s,dim)

      implicit none

      integer     dim,i
      complex*16  u(dim),v(dim),s

      s = dconjg(u(1))*v(1)
      do i = 2,dim
         s = s+dconjg(u(i))*v(i)
      enddo

      return
      end


C-----------------------------------------------------------------------
C Library subroutine vvxazz 
C
C       u(i)*dconjg(v(j)) = m(i,j)
C       Dimensions: d(dim),v(dim),m(dim,dim)
C-----------------------------------------------------------------------
      subroutine vvxazz(u,v,m,dim)
 
      implicit none
      
      integer dim,i,j
      complex*16 u(dim),v(dim),m(dim,dim)

      do i=1,dim
         do j=1,dim
            m(i,j) = u(i)*dconjg(v(j))
         enddo
      enddo
      
      return
      end



C-----------------------------------------------------------------------
C Library subroutine qvxxdd
C
C Multiplication of a real quadratic matrix with a real vector
C     a(i,j)*v(j)=w(i)
C-----------------------------------------------------------------------
      subroutine qvxxdd (a,v,w,dim)

      implicit none

      integer dim,i,j
      real*8  a(dim,dim),v(dim),w(dim)

      do i = 1,dim
         w(i) = a(i,1)*v(1)
      enddo
      do j = 2,dim
         do i = 1,dim
            w(i) = w(i)+a(i,j)*v(j)
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine qvxxzz
C
C Multiplication of a complex quadratic matrix with a complex vector
C     a(i,j)*v(j)=w(i)
C-----------------------------------------------------------------------

      subroutine qvxxzz (a,v,w,dim)

      implicit none

      integer     dim,i,j
      complex*16  a(dim,dim),v(dim),w(dim)

      do i = 1,dim
         w(i) = a(i,1)*v(1)
      enddo
      do j = 2,dim
         do i = 1,dim
            w(i) = w(i)+a(i,j)*v(j)
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine qvxxdd1
C
C Multiplication of a real quadratic matrix with a real vector
C     a(i,j)*v(j)=w(i)
C
C NB phdim is physical (leading) dimension, dim is used dimension 
C-----------------------------------------------------------------------

      subroutine qvxxdd1 (a,v,w,phdim,dim)

      implicit none

      integer phdim,dim,i,j
      real*8  a(phdim,dim),v(dim),w(dim)

      do i = 1,dim
         w(i) = a(i,1)*v(1)
      enddo
      do j = 2,dim
         do i = 1,dim
            w(i) = w(i)+a(i,j)*v(j)
         enddo
      enddo

      return
      end

C ----------------------------------------------------------------------
C Library subroutine qvtxdd
C
C Multiplication of the transpose of a real quadratic matrix with a 
C real vector
C     a(j,i)*v(j)=w(i)
C-----------------------------------------------------------------------

      subroutine qvtxdd (a,v,w,dim)

      implicit none

      integer dim,i,j
      real*8  a(dim,dim),v(dim),w(dim)

      do i = 1,dim
         w(i) = a(1,i)*v(1)
         do j = 2,dim
            w(i) = w(i)+a(j,i)*v(j)
         enddo
      enddo

      return
      end

C ----------------------------------------------------------------------
C Library subroutine qvtxdd1
C
C Multiplication of the transpose of a real quadratic matrix with a 
C real vector
C     a(j,i)*v(j)=w(i)
C
C NB phdim is physical (leading) dimension, dim is used dimension 
C-----------------------------------------------------------------------

      subroutine qvtxdd1 (a,v,w,phdim,dim)

      implicit none

      integer phdim,dim,i,j
      real*8  a(phdim,dim),v(dim),w(dim)

      do i = 1,dim
         w(i) = a(1,i)*v(1)
         do j = 2,dim
            w(i) = w(i)+a(j,i)*v(j)
         enddo
      enddo

      return
      end

C ----------------------------------------------------------------------
C Library subroutine mvtxdd
C
C Multiplication of the transpose of a real rectangular matrix with a
C real vector
C     a(j,i)*v(j)=w(i)
C
C-----------------------------------------------------------------------

      subroutine mvtxdd (a,v,w,dim1,dim2)

      implicit none

      integer dim1,dim2,i,j
      real*8  a(dim1,dim2),v(dim1),w(dim2)

      do i = 1,dim2
         w(i) = a(1,i)*v(1)
         do j = 2,dim1
            w(i) = w(i)+a(j,i)*v(j)
         enddo
      enddo

      return
      end

C ----------------------------------------------------------------------
C Library subroutine mvtxdd1
C
C Multiplication of the transpose of a real rectangular matrix with a
C real vector
C     a(j,i)*v(j)=w(i)
C
C NB phdim is physical (leading) dimension, dim1-2 is used dimension 
C-----------------------------------------------------------------------

      subroutine mvtxdd1 (a,v,w,phdim,dim1,dim2)

      implicit none

      integer phdim,dim1,dim2,i,j
      real*8  a(phdim,phdim),v(dim1),w(dim2)

      do i = 1,dim2
         w(i) = a(1,i)*v(1)
         do j = 2,dim1
            w(i) = w(i)+a(j,i)*v(j)
         enddo
      enddo

      return
      end

C ----------------------------------------------------------------------
C Library subroutine mvtxzz
C
C Multiplication of the transpose of a complex rectangular matrix with a
C complex vector
C     a(j,i)*v(j)=w(i)
C-----------------------------------------------------------------------

      subroutine mvtxzz (a,v,w,dim1,dim2)

      implicit none

      integer dim1,dim2,i,j
      complex*16  a(dim1,dim2),v(dim1),w(dim2)

      do i = 1,dim2
         w(i) = a(1,i)*v(1)
         do j = 2,dim1
            w(i) = w(i)+a(j,i)*v(j)
         enddo
      enddo

      return
      end


C ----------------------------------------------------------------------
C Library subroutine mvxxzz
C
C Multiplication of a complex rectangular matrix with a complex vector.
C     a(i,j)*v(j)=w(i)
C-----------------------------------------------------------------------

      subroutine mvxxzz (a,v,w,dim1,dim2)

      implicit none

      integer dim1,dim2,i,j
      complex*16  a(dim1,dim2),v(dim2),w(dim1)

      do i = 1,dim1
         w(i) = a(i,1)*v(1)
      enddo
      do j = 2,dim2
         do i = 1,dim1
            w(i) = w(i)+a(i,j)*v(j)
         enddo
      enddo

      return
      end

C ----------------------------------------------------------------------
C Library subroutine dvxxdd
C
C Multiplication of a diagonal real matrix with a real vector
C     a(i)*v(i)=w(i)
C-----------------------------------------------------------------------

      subroutine dvxxdd (a,v,w,dim)

      implicit none

      integer dim,i
      real*8  a(dim),v(dim),w(dim)

      do i = 1,dim
         w(i) = a(i)*v(i)
      enddo

      return
      end

C ----------------------------------------------------------------------
C Library subroutine dvxxdz
C
C Multiplication of a diagonal real matrix with a complex vector
C     a(i)*v(i)=w(i)
C-----------------------------------------------------------------------

      subroutine dvxxdz (a,v,w,dim)

      implicit none

      integer dim,i
      real*8  a(dim)
      complex*16  v(dim),w(dim)

      do i = 1,dim
         w(i) = a(i)*v(i)
      enddo

      return
      end

C ----------------------------------------------------------------------
C Library subroutine dvxxdzo
C
C Multiplication of a diagonal real matrix with a complex vector
C     a(i)*v(i)=v(i)
C-----------------------------------------------------------------------

      subroutine dvxxdzo (a,v,dim)

      implicit none

      integer dim,i
      real*8  a(dim)
      complex*16 v(dim)

      do i = 1,dim
         v(i) = a(i)*v(i)
      enddo

      return
      end

C ----------------------------------------------------------------------
C Library subroutine dvxxddo
C
C Multiplication of a diagonal real matrix with a real vector
C     a(i)*v(i)=v(i)
C-----------------------------------------------------------------------

      subroutine dvxxddo (a,v,dim)

      implicit none

      integer dim,i
      real*8  a(dim),v(dim)

      do i = 1,dim
         v(i) = a(i)*v(i)
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine qvxxdz
C
C Multiplication of a real quadratic matrix with a complex vector
C     a(i,j)*v(j)=w(i)
C-----------------------------------------------------------------------

      subroutine qvxxdz (a,v,w,dim)

      implicit none

      integer dim,i,j
      real*8     a(dim,dim)
      complex*16 v(dim),w(dim)

      do i = 1,dim
         w(i) = a(i,1)*v(1)
      enddo
      do j = 2,dim
         do i = 1,dim
            w(i) = w(i)+a(i,j)*v(j)
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine qvxxdz1
C
C Multiplication of a real quadratic matrix with a complex vector
C     a(i,j)*v(j)=w(i)
C
C NB phdim is physical (leading) dimension, dim is used dimension 
C-----------------------------------------------------------------------

      subroutine qvxxdz1 (a,v,w,phdim,dim)

      implicit none

      integer    phdim,dim,i,j
      real*8     a(phdim,dim)
      complex*16 v(dim),w(dim)

      do i = 1,dim
         w(i) = a(i,1)*v(1)
      enddo
      do j = 2,dim
         do i = 1,dim
            w(i) = w(i)+a(i,j)*v(j)
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine qvxxzz1
C
C Multiplication of a complex quadratic matrix with a complex vector
C     a(i,j)*v(j)=w(i)
C
C NB phdim is physical (leading) dimension, dim is used dimension 
C-----------------------------------------------------------------------

      subroutine qvxxzz1 (a,v,w,phdim,dim)

      implicit none

      integer    phdim,dim,i,j
      complex*16 a(phdim,dim),v(dim),w(dim)

      do i = 1,dim
         w(i) = a(i,1)*v(1)
      enddo
      do j = 2,dim
         do i = 1,dim
            w(i) = w(i)+a(i,j)*v(j)
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine qvxxyz
C
C Multiplication of a complex quadratic matrix with a complex vector,
C where the complex matrix is stored as 2 real matrices.
C     a(i,j)*v(j)=w(i)
C-----------------------------------------------------------------------

      subroutine qvxxyz (a,v,w,dim)

      implicit none

      integer dim,i,j
      real*8     a(dim,dim,2)
      complex*16 v(dim),w(dim)

      do i = 1,dim
         w(i) = a(i,1,1)*v(1)+(0,1)*a(i,1,2)*v(1)
      enddo
      do j = 2,dim
         do i = 1,dim
            w(i) = w(i)+a(i,j,1)*v(j)+(0,1)*a(i,j,2)*v(j)
         enddo
      enddo

      return
      end


