C **********************************************************************
C                                                                     
C                              MMLIB                                  
C                                                                     
C  Library module containing linear algebra routines that involve the 
C  multiplication of matrices
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
C       e.g. 'qm' denotes the operation (quadratic matrix * rectangular
C       matrix)
C    Character 3 denotes how first object is used:
C       x: unchanged from input
C       t: transpose of input
C       a: adjoint of input
C       c: complex conjugate of input
C       v: as a vector
C    Character 4 denotes how second object is used:
C       see Character 3 above.
C    Character 5, 6 denote data types of first, second object
C       respectively:
C       z: complex double precision (complex*16)
C       c: complex single precision (complex*8)
C       d: real double precision (real*8)
C       r: real single precision (real*4)
C     Further characters, if present, give more informaion:
C       a: the result is added to a further object
C       r: the result is subtracted (removed) from a further object
C       c: the input matrices commute
C       h: the resulting matrix is hermitian
C       h1: the resulting matrix is anti-hermitian
C       s: the resulting matrix is symmetric
C       1: the physical dimensions of the matrices differs from those
C          used.
C   
C  Contents:                                                          
C  In the following list of available subroutines, matrices on the LHS
C  of the definition are input, that on the RHS output. The usual
C  summation convention is used i.e. a sum is made over repeated indices
C  on the LHS (NOTE: there is only elementwise multiplication and 
C  no subsequent summation if a diagonal matrix is involved !!!).
C
C   qqxxdd (a,b,c,dim)
C       Definition: a(j,k)*b(k,i) = c(j,i) .
C       Dimensions: a(dim,dim),b(dim,dim),c(dim,dim)
C
C   qqxxzd (a,b,c,dim)
C       Definition: a(j,k)*b(k,i) = c(j,i) .
C       Dimensions: a(dim,dim),b(dim,dim),c(dim,dim)
C
C   mmxxzz (a,b,c,dim1,dim2,dim3)
C       Definition: a(j,k)*b(k,i) = c(j,i) .
C       Dimensions: a(dim1,dim2),b(dim2,dim3),c(dim1,dim3)
C
C   mmvxzz (a,b,c,dim1,dim2,dim3,p)
C       Definition: a(p,k)*b(k,i) = c(i)
C       Dimensions: a(dim1,dim2),b(dim2,dim3),c(dim3)
C
C   mmxtzz (a,b,c,dim1,dim2,dim3)
C       Definition: a(j,k)*b(i,k) = c(j,i) .
C       Dimensions: a(dim1,dim2),b(dim3,dim2),c(dim1,dim3)
C
C   mmtczz (a,b,c,dim1,dim2,dim3)
C       Definition: a(k,j)*dconjg(b(k,i)) = c(j,i) .
C       Dimensions: a(dim1,dim2),b(dim1,dim3),c(dim2,dim3)
C
C   mmxxzza (a,b,c,dim1,dim2,dim3)
C       Definition: a(j,k)*b(k,i) + c(j,i) = c(j,i) .
C       Dimensions: a(dim1,dim2),b(dim2,dim3),c(dim1,dim3)
C
C   dmxxzz (a,b,c,dim1,dim2)
C       Definition: a(j)*b(j,i) = c(j,i) .
C       Dimensions: a(dim1),b(dim1,dim2),c(dim1,dim2)
C                                                                    
C   dmxxzza (a,b,c,dim1,dim2)
C       Definition: a(j)*b(j,i) + c(j,i) = c(j,i) .
C       Dimensions: a(dim1),b(dim1,dim2),c(dim1,dim2)
C                                                                    
C   dmxxdz (a,b,c,dim1,dim2)
C       Definition: a(j)*b(j,i) = c(j,i) .
C       Dimensions: a(dim1),b(dim1,dim2),c(dim1,dim2)
C
C   dmxxdd (a,b,c,dim1,dim2)
C       Definition: a(j)*b(j,i) = c(j,i) .
C       Dimensions: a(dim1),b(dim1,dim2),c(dim1,dim2)
C
C   ddxxdd (a,b,c,dim)
C       Definition: a(j)*b(j) = c(j) .
C       Dimensions: a(dim),b(dim),c(dim)
C
C   dqxxzz (a,b,c,dim)
C       Definition: a(j)*b(j,i) = c(j,i) .
C       Dimensions: a(dim),b(dim,dim),c(dim,dim)
C
C   dqxxdd (a,b,c,dim)
C       Definition: a(j)*b(j,i) = c(j,i) .
C       Dimensions: a(dim),b(dim,dim),c(dim,dim)
C
C   dqxxdd1 (a,b,c,phdim,dim)
C       Definition: a(j)*b(j,i) = c(j,i) .
C       Dimensions: a(dim),b(phdim,dim),c(phdim,dim)
C
C   dqxxdz (a,v,w,dim)
C       Definition: a(i)*v(i,j) = w(i,j) .
C       Dimensions: a(dim),v(dim,dim),w(dim,dim)
C
C   dqxxdz2 (a,v,w,dim)
C       Definition: a(j)*v(i,j) = w(i,j) .
C       Dimensions: a(dim),v(dim,dim),w(dim,dim)
C
CC  mmxxdz (a,b,c,dim1,dim2,dim3)
CC      Definition: a(j,k)*b(k,i) = c(j,i) .
CC      Dimensions: a(dim1,dim2),b(dim2,dim3),c(dim1,dim3)
C
C   mmtxdd (a,b,c,dim1,dim2,dim3)
C       Definition: a(k,j)*b(k,i) = c(j,i) .
C       Dimensions: a(dim2,dim1),b(dim2,dim3),c(dim1,dim3)
C
C   mmtxdd1 (a,b,c,phdim,dim1,dim2,dim3)
C       Definition: a(k,j)*b(k,i) = c(j,i) .
C       Dimensions: a(dim2,dim1),b(dim2,dim3),c(dim1,dim3) (= used dim.)
C       Dimensions: a,b,c(phdim,phdim) (= allocated dimension)
C
C   mmxtdd (a,b,c,dim1,dim2,dim3)
C       Definition: a(j,k)*b(i,k) = c(j,i) .
C       Dimensions: a(dim1,dim2),b(dim3,dim2),c(dim1,dim3)
C
C   mmxxdd (a,b,c,dim1,dim2,dim3)
C       Definition: a(j,k)*b(k,i) = c(j,i) .
C       Dimensions: a(dim1,dim2),b(dim2,dim3),c(dim1,dim3)
C
C   mmxxdd1 (a,b,c,phdim,dim1,dim2,dim3)
C       Definition: a(j,k)*b(k,i) = c(j,i) .
C       Dimensions: a(dim1,dim2),b(dim2,dim3),c(dim1,dim3) (= used dim.)
C       Dimensions: a,b,c(phdim,phdim) (= allocated dimension)
C
C   qqxtdd (a,b,c,dim)
C       Definition: a(j,k)*b(i,k) = c(j,i) .
C       Dimensions: a(dim,dim),b(dim,dim),c(dim,dim)
C
C   qqxtdd1 (a,b,c,phdim,dim)
C       Definition: a(j,k)*b(i,k) = c(j,i) ; 1 <= i,j,k <= dim .
C       Dimensions: a(phdim,dim),b(phdim,dim),c(phdim,dim)
C
C   qqtxdd (a,b,c,dim)
C       Definition: a(k,j)*b(k,i) = c(j,i) .
C       Dimensions: a(dim,dim),b(dim,dim),c(dim,dim)
C
C   qdxxzz (a,b,c,dim)
C       Definition: a(j,i)*b(i) = c(j,i) .
C       Dimensions: a(dim,dim),b(dim),c(dim,dim)
C
C   qdxxdd (a,b,c,dim)
C       Definition: a(j,i)*b(i) = c(j,i) .
C       Dimensions: a(dim,dim),b(dim),c(dim,dim)
C
C   qdxxdd1 (a,b,c,phdim,dim)
C       Definition: a(j,i)*b(i) = c(j,i) ; 1 <= i,j <= dim .
C       Dimensions: a(phdim,dim),b(dim),c(phdim,dim)
C
C   hhxtzzc (a,b,c,dim)
C       Definition: a(j,k)*b(i,k) = c(j,i) .
C       Dimensions: a(dim,dim),b(dim,dim),c(dim,dim)
C
C   mmaxzzh (a,b,c,dim1,dim2)
C       Definition: dconjg(a(k,j))*b(k,i) = c(j,i) .
C       Dimensions: a(dim1,dim2),b(dim1,dim2),c(dim2,dim2)
C
C   mmaxzzh1 (a,b,c,dim1,dim2)
C       Definition: dconjg(a(k,j))*b(k,i) = c(j,i) .
C       Dimensions: a(dim1,dim2),b(dim1,dim2),c(dim2,dim2)
C
C   mmtxzzs (a,b,c,dim1,dim2)
C       Definition: a(k,j)*b(k,i) = c(j,i) .
C       Dimensions: a(dim1,dim2),b(dim1,dim2),c(dim2,dim2)
C
C   qmxxzz (a,b,c,dim1,dim2)
C       Definition: a(j,k)*b(k,i) = c(j,i) .
C       Dimensions: a(dim1,dim1),b(dim1,dim2),c(dim1,dim2)
C
C   qqxxdd1 (a,b,c,phdim,dim)
C       Definition: a(j,k)*b(k,i) = c(j,i) ; 1 <= i,j,k <= dim .
C       Dimensions: a(phdim,dim),b(phdim,dim),c(phdim,dim)
C
CC  mmxxzz1 (a,b,c,phdim1,phdim2,phdim3,dim1,dim2,dim3)
CC      Definition: a(j,k)*b(k,i) = c(j,i) .
CC      Dimensions: a(phdim1,phdim2),b(phdim2,phdim3),c(phdim1,phdim3)
C
C   mmaxzz (a,b,c,dim1,dim2,dim3)
C       Definition: dconjg(a(k,j))*b(k,i) = c(j,i) .
C       Dimensions: a(dim1,dim2),b(dim1,dim3),c(dim2,dim3)
c
c   mmcxzz (a,b,c,dim1,dim2,dim3)
C       Definition: dconjg(a(j,k))*b(k,i) = c(j,i) .
C       Dimensions: a(dim1,dim2),b(dim2,dim3),c(dim1,dim3)
C
c   qqcxzz (a,b,c,dim)
C       Definition: dconjg(a(j,k))*b(k,i) = c(j,i) .
C       Dimensions: a(dim,dim),b(dim,dim),c(dim,dim)
C
C   mmxazz (a,b,c,dim1,dim2,dim3)
C       Definition: a(j,k)*dconjg(b(i,k)) = c(j,i) .
C       Dimensions: a(dim1,dim2),b(dim3,dim2),c(dim1,dim3)
C
C   mmxtzza (a,b,c,dim1,dim2,dim3)
C       Definition: a(j,k)*b(i,k) +c(j,i) = c(j,i) .
C       Dimensions: a(dim1,dim2),b(dim3,dim2),c(dim1,dim3)
C
C   mqxtzza (a,b,c,dim1,dim2)
C       Definition: a(j,k)*b(i,k) +c(j,i) = c(j,i) .
C       Dimensions: a(dim1,dim2),b(dim2,dim2),c(dim1,dim2)
C
C   mmcazza (a,b,c,dim1,dim2,dim3)
C       Definition: conjg(a(j,k))*conjg(b(i,k)) + c(j,i) = c(j,i)
C       Dimensions: a(dim1,dim2),b(dim3,dim2),c(dim1,dim3/)
C
C   mqxxzza (a,b,c,dim1,dim2)
C       Definition: a(j,k)*b(k,i) +c(j,i) = c(j,i) .
C       Dimensions: a(dim1,dim2),b(dim2,dim2),c(dim1,dim2)
C
C   mqxtzd (a,b,c,dim1,dim2)
C       Definition: a(j,k)*b(i,k) = c(j,i) .
C       Dimensions: a(dim1,dim2),b(dim2,dim2),c(dim1,dim2)
C
C   qqtxzz (a,b,c,dim)
C       Definition: a(k,j)*b(k,i) = c(j,i) .
C       Dimensions: a(dim,dim),b(dim,dim),c(dim,dim)
C
C   mqxazz (a,b,c,dim1,dim2)
C       Definition: a(j,k)*dconjg(b(i,k)) = c(j,i) .
C       Dimensions: a(dim1,dim2),b(dim2,dim2),c(dim1,dim2)
C
C   qmxxdz (a,b,c,dim1,dim2)
C       Definition: a(j,k)*b(k,i) = c(j,i) .
C       Dimensions: a(dim1,dim1),b(dim1,dim2),c(dim1,dim2)
C
C   qmtxdz (a,b,c,dim1,dim2)
C       Definition: a(k,j)*b(k,i) = c(j,i) .
C       Dimensions: a(dim1,dim1),b(dim1,dim2),c(dim1,dim2)
C
C   qqxxzz (a,b,c,dim)
C       Definition: a(j,k)*b(k,i) = c(j,i) .
C       Dimensions: a(dim,dim),b(dim,dim),c(dim,dim)
C
C   mqxtzz (a,b,c,dim1,dim2)
C       Definition: a(j,k)*b(i,k) = c(j,i) .
C       Dimensions: a(dim1,dim2),b(dim2,dim2),c(dim1,dim2)
C
C   mqxxzz (a,b,c,dim1,dim2)
C       Definition: a(j,k)*b(k,i) = c(j,i) .
C       Dimensions: a(dim1,dim2),b(dim2,dim2),c(dim1,dim2)
C
C   mmtxzz (a,b,c,dim1,dim2,dim3)
C       Definition: a(k,j)*b(k,i) = c(j,i) .
C       Dimensions: a(dim1,dim2),b(dim1,dim3),c(dim2,dim3)
C
C   mmtxzza (a,b,c,dim1,dim2,dim3)
C       Definition: a(k,j)*b(k,i) + c(j,i) = c(j,i) .
C       Dimensions: a(dim1,dim2),b(dim1,dim3),c(dim2,dim3)
C                                                                  
C   mmxczz (a,b,c,dim1,dim2,dim3)
C       Definition: a(j,k)*dconjg(b(k,i)) = c(j,i) .
C       Dimensions: a(dim1,dim2),b(dim2,dim3),c(dim1,dim3)
C
C   mqxczz (a,b,c,dim1,dim2)
C       Definition: a(j,k)*dconjg(b(k,i)) = c(j,i) .
C       Dimensions: a(dim1,dim2),b(dim2,dim2),c(dim1,dim2)
C
C   qqaxzz (a,b,c,dim)
C       Definition: dconjg(a(k,j))*b(k,i) = c(j,i) .
C       Dimensions: a(dim,dim),b(dim,dim),c(dim,dim)      
C
C***********************************************************************


C-----------------------------------------------------------------------
C Library subroutine qqxxdd
C
C Multiplication of two real quadratic matrices:
C     a(j,k)*b(k,i) = c(j,i)
C-----------------------------------------------------------------------

      subroutine qqxxdd (a,b,c,dim)

      implicit none

      integer dim,i,j,k
      real*8  a(dim,dim),b(dim,dim),c(dim,dim)

      do i=1,dim
         do j=1,dim
            c(j,i) = a(j,1)*b(1,i)
         enddo
      enddo
      do i = 1,dim
         do k = 2,dim
            do j = 1,dim
               c(j,i) = c(j,i)+a(j,k)*b(k,i)
            enddo
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine qqxxzd
C
C Multiplication of a complex quadratic matrix with a real quadratic
C matrix:
C     a(j,k)*b(k,i) = c(j,i)
C-----------------------------------------------------------------------

      subroutine qqxxzd (a,b,c,dim)

      implicit none

      integer dim,i,j,k
      real*8  b(dim,dim)
      complex*16 a(dim,dim),c(dim,dim)

      do i=1,dim
         do j=1,dim
            c(j,i) = a(j,1)*b(1,i)
         enddo
      enddo
      do i = 1,dim
         do k = 2,dim
            do j = 1,dim
               c(j,i) = c(j,i)+a(j,k)*b(k,i)
            enddo
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine mmxxzz
C
C Multiplication of two complex rectangular matrices:
C     a(j,k)*b(k,i) = c(j,i)
C-----------------------------------------------------------------------

      subroutine mmxxzz (a,b,c,dim1,dim2,dim3)

      implicit none

      integer     dim1,dim2,dim3,i,j,k
      complex*16  a(dim1,dim2),b(dim2,dim3),c(dim1,dim3)

      do i=1,dim3
         do j=1,dim1
            c(j,i) = a(j,1)*b(1,i)
         enddo
      enddo
      do i = 1,dim3
         do k = 2,dim2
            do j = 1,dim1
               c(j,i) = c(j,i)+a(j,k)*b(k,i)
            enddo
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C subroutine mmvxzz
C
C Multiplication of two complex rectangular matrices with p fixed (uses
C first matrix as a vector):
C 
C     a(p,k)*b(k,i) = c(i)
C-----------------------------------------------------------------------

      subroutine mmvxzz (a,b,c,dim1,dim2,dim3,p)

      implicit none

      integer     dim1,dim2,dim3,p,i,k
      complex*16  a(dim1,dim2),b(dim2,dim3),c(dim3)

      do i=1,dim3
         c(i) = a(p,1)*b(1,i)
      enddo
      do i = 1,dim3
         do k = 2,dim2
            c(i) = c(i)+a(p,k)*b(k,i)
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine mmxtzz
C
C Multiplication of a complex rectangular matrix with the transpose of
C a rectangular complex matrix
C     a(j,k)*b(i,k) = c(j,i)
C-----------------------------------------------------------------------

      subroutine mmxtzz (a,b,c,dim1,dim2,dim3)

      implicit none

      integer     dim1,dim2,dim3,i,j,k
      complex*16  a(dim1,dim2),b(dim3,dim2),c(dim1,dim3)

      do i=1,dim3
         do j=1,dim1
            c(j,i) = a(j,1)*b(i,1)
         enddo
      enddo

      do i = 1,dim3
         do k = 2,dim2
            do j = 1,dim1
               c(j,i) = c(j,i)+a(j,k)*b(i,k)
            enddo
         enddo
      enddo

      return
      end


C-----------------------------------------------------------------------
C Library subroutine mmxxzza
C
C Multiplication of a complex rectangular matrix with a rectangular 
C complex matrix, the result of which is added to a further matrix.
C     a(j,k)*b(k,i) +c(j,i) = c(j,i)
C-----------------------------------------------------------------------

      subroutine mmxxzza (a,b,c,dim1,dim2,dim3)

      implicit none

      integer     dim1,dim2,dim3,i,j,k
      complex*16  a(dim1,dim2),b(dim2,dim3),c(dim1,dim3)

      do k = 1,dim2
         do i = 1,dim3
            do j = 1,dim1
               c(j,i) = c(j,i)+a(j,k)*b(k,i)
            enddo
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine dmxxzz
C
C Multiplication of a diagonal complex matrix with a complex rectangular
C matrix:
C     a(j)*b(j,i) = c(j,i)
C-----------------------------------------------------------------------

      subroutine dmxxzz (a,b,c,dim1,dim2)

      implicit none

      integer     dim1,dim2,i,j
      complex*16  a(dim1),b(dim1,dim2),c(dim1,dim2)

      do i = 1,dim2
         do j = 1,dim1
            c(j,i) = a(j)*b(j,i)
         enddo
      enddo

      return
      end


C-----------------------------------------------------------------------
C Library subroutine dmxxzza
C
C Multiplication of a diagonal complex matrix with a complex rectangular
C matrix:
C     a(j)*b(j,i) + c(j,i) = c(j,i)
C-----------------------------------------------------------------------

      subroutine dmxxzza (a,b,c,dim1,dim2)

      implicit none

      integer     dim1,dim2,i,j
      complex*16  a(dim1),b(dim1,dim2),c(dim1,dim2)

      do i = 1,dim2
         do j = 1,dim1
            c(j,i) = c(j,i) + a(j)*b(j,i)
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine dmxxdz
C
C Multiplication of a diagonal real matrix with a real rectangular
C matrix:
C     a(j)*b(j,i) = c(j,i)
C-----------------------------------------------------------------------

      subroutine dmxxdz (a,b,c,dim1,dim2)

      implicit none

      integer     dim1,dim2,i,j
      real*8      a(dim1)
      complex*16  b(dim1,dim2),c(dim1,dim2)

      do i = 1,dim2
         do j = 1,dim1
            c(j,i) = a(j)*b(j,i)
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine dmxxdd
C
C Multiplication of a diagonal real matrix with a real rectangular
C matrix:
C     a(j)*b(j,i) = c(j,i)
C-----------------------------------------------------------------------

      subroutine dmxxdd (a,b,c,dim1,dim2)

      implicit none

      integer     dim1,dim2,i,j
      real*8      a(dim1),b(dim1,dim2),c(dim1,dim2)

      do i = 1,dim2
         do j = 1,dim1
            c(j,i) = a(j)*b(j,i)
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine dqxxzz
C
C Multiplication of a diagonal complex matrix with a complex quadratic 
C matrix:
C     a(j)*b(j,i) = c(j,i)
C-----------------------------------------------------------------------

      subroutine dqxxzz (a,b,c,dim)

      implicit none

      integer     dim,i,j
      complex*16  a(dim),b(dim,dim),c(dim,dim)

      do i = 1,dim
         do j = 1,dim
            c(j,i) = a(j)*b(j,i)
         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
C Library subroutine ddxxdd
C
C Multiplication of a diagonal real matrix with a diagonal real matrix:
C     a(j)*b(j) = c(j)
C-----------------------------------------------------------------------

      subroutine ddxxdd (a,b,c,dim)

      implicit none

      integer     dim,j
      real*8      a(dim),b(dim),c(dim)

      do j = 1,dim
         c(j) = a(j)*b(j)
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine dqxxdd
C
C Multiplication of a diagonal real matrix with a real quadratic matrix:
C     a(j)*b(j,i) = c(j,i)
C-----------------------------------------------------------------------

      subroutine dqxxdd (a,b,c,dim)

      implicit none

      integer     dim,i,j
      real*8      a(dim),b(dim,dim),c(dim,dim)

      do i = 1,dim
         do j = 1,dim
            c(j,i) = a(j)*b(j,i)
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine dqxxdd1
C
C Multiplication of a diagonal real matrix with a real quadratic matrix:
C     a(j)*b(j,i) = c(j,i)
C
C NB phdim is physical (leading) dimension, dim is used dimension 
C-----------------------------------------------------------------------

      subroutine dqxxdd1 (a,b,c,phdim,dim)

      implicit none

      integer phdim,dim,i,j
      real*8  a(dim),b(phdim,dim),c(phdim,dim)

      do i = 1,dim
         do j = 1,dim
            c(j,i) = a(j)*b(j,i)
         enddo
      enddo

      return
      end

C ----------------------------------------------------------------------
C Library subroutine dqxxdz
C
C Multiplication of a diagonal real matrix with a complex matrix
C     a(i)*v(i,j)=w(i,j)
C-----------------------------------------------------------------------

      subroutine dqxxdz (a,v,w,dim)

      implicit none

      integer dim,i,j
      real*8  a(dim)
      complex*16  v(dim,dim),w(dim,dim)

      do j=1,dim
         do i=1,dim
            w(i,j) = a(i)*v(i,j)
         enddo
      enddo

      return
      end


C ----------------------------------------------------------------------
C Library subroutine dqxxdz2
C
C Multiplication of a diagonal real matrix with a complex matrix 
C     a(j)*v(i,j)=w(i,j)
C-----------------------------------------------------------------------

      subroutine dqxxdz2 (a,v,w,dim)

      implicit none

      integer dim,i,j
      real*8  a(dim)
      complex*16  v(dim,dim),w(dim,dim)

      do j=1,dim
         do i=1,dim
            w(i,j) = a(j)*v(i,j)
         enddo
      enddo

      return
      end


C-----------------------------------------------------------------------
C Library subroutine mmxxdz
C
C Multiplication of a rectangular real matrix with a rectangular complex
C matrix:
C     a(j,k)*b(k,i) = c(j,i)
C-----------------------------------------------------------------------

C     subroutine mmxxdz (a,b,c,dim1,dim2,dim3)

C     implicit none

C     integer     dim1,dim2,dim3,i,j,k
C     real*8      a(dim1,dim2)
C     complex*16  b(dim2,dim3),c(dim1,dim3)

C     do i=1,dim3
C        do j=1,dim1
C           c(j,i) = a(j,1)*b(1,i)
C        enddo
C     enddo
C     do i = 1,dim3
C        do k = 2,dim2
C           do j = 1,dim1
C              c(j,i) = c(j,i)+a(j,k)*b(k,i)
C           enddo
C        enddo
C     enddo

C     return
C     end

C-----------------------------------------------------------------------
C Library subroutine mmtxdd
C
C Multiplication of the transpose of a rectangular real matrix with a 
C rectangular real matrix:
C     a(k,j)*b(k,i) = c(j,i)
C-----------------------------------------------------------------------

      subroutine mmtxdd (a,b,c,dim1,dim2,dim3)

      implicit none

      integer dim1,dim2,dim3,i,j,k
      real*8  a(dim2,dim1),b(dim2,dim3),c(dim1,dim3)

      do i = 1,dim3
         do j = 1,dim1
            c(j,i) = a(1,j)*b(1,i)
            do k = 2,dim2
               c(j,i) = c(j,i)+a(k,j)*b(k,i)
            enddo
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine mmxtdd
C
C Multiplication of the rectangular real matrix with the transpose of a
C rectangular real matrix:
C       a(j,k)*b(i,k) = c(j,i)
C-----------------------------------------------------------------------
 
      subroutine mmxtdd(a,b,c,dim1,dim2,dim3)
 
      implicit none
 
      integer dim1,dim2,dim3,i,j,k
      real*8  a(dim1,dim2),b(dim3,dim2),c(dim1,dim3)
 
      do i = 1,dim3
         do j = 1,dim1
            c(j,i) = a(j,1)*b(i,1)
            do k = 2,dim2
               c(j,i) = c(j,i)+a(j,k)*b(i,k)
            enddo
         enddo
      enddo
 
      return
      end

C-----------------------------------------------------------------------
C Library subroutine mmtxdd1
C
C Multiplication of the transpose of a rectangular real matrix with a 
C rectangular real matrix:
C     a(k,j)*b(k,i) = c(j,i)
C
C NB phdim is physical (leading) dimension, dim1-3 is used dimension 
C-----------------------------------------------------------------------

      subroutine mmtxdd1 (a,b,c,phdim,dim1,dim2,dim3)

      implicit none

      integer phdim,dim1,dim2,dim3,i,j,k
      real*8  a(phdim,phdim),b(phdim,phdim),c(phdim,phdim)

      do i = 1,dim3
         do j = 1,dim1
            c(j,i) = a(1,j)*b(1,i)
            do k = 2,dim2
               c(j,i) = c(j,i)+a(k,j)*b(k,i)
            enddo
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine mmxxdd
C
C Multiplication of two real rectangular matrices:
C     a(j,k)*b(k,i) = c(j,i)
C-----------------------------------------------------------------------

      subroutine mmxxdd (a,b,c,dim1,dim2,dim3)

      implicit none

      integer dim1,dim2,dim3,i,j,k
      real*8  a(dim1,dim2),b(dim2,dim3),c(dim1,dim3)

      do i=1,dim3
         do j=1,dim1
            c(j,i) = a(j,1)*b(1,i)
         enddo
      enddo
      do i = 1,dim3
         do k = 2,dim2
            do j = 1,dim1
               c(j,i) = c(j,i)+a(j,k)*b(k,i)
            enddo
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine mmxxdd1
C
C Multiplication of two real rectangular matrices:
C     a(j,k)*b(k,i) = c(j,i)
C
C NB phdim is physical (leading) dimension, dim1-3 is used dimension 
C-----------------------------------------------------------------------

      subroutine mmxxdd1 (a,b,c,phdim,dim1,dim2,dim3)

      implicit none

      integer phdim,dim1,dim2,dim3,i,j,k
      real*8  a(phdim,phdim),b(phdim,phdim),c(phdim,phdim)

      do i=1,dim3
         do j=1,dim1
            c(j,i) = a(j,1)*b(1,i)
         enddo
      enddo
      do i = 1,dim3
         do k = 2,dim2
            do j = 1,dim1
               c(j,i) = c(j,i)+a(j,k)*b(k,i)
            enddo
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine qqxtdd
C
C Multiplication of a real quadratic matrix with the transpose of a real
C quadratic matrix:
C     a(j,k)*b(i,k) = c(j,i)
C-----------------------------------------------------------------------

      subroutine qqxtdd (a,b,c,dim)

      implicit none

      integer dim,i,j,k
      real*8  a(dim,dim),b(dim,dim),c(dim,dim)

      do i = 1,dim
         do j = 1,dim
            c(j,i) = a(j,1)*b(i,1)
         enddo
      enddo
      do k = 2,dim
         do i = 1,dim
            do j = 1,dim
               c(j,i) = c(j,i)+a(j,k)*b(i,k)
            enddo
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine qqxtdd1
C
C Multiplication of a real quadratic matrix with the transpose of a real
C quadratic matrix:
C     a(j,k)*b(i,k) = c(j,i)
C
C NB phdim is physical (leading) dimension, dim is used dimension 
C-----------------------------------------------------------------------

      subroutine qqxtdd1 (a,b,c,phdim,dim)

      implicit none

      integer phdim,dim,i,j,k
      real*8  a(phdim,dim),b(phdim,dim),c(phdim,dim)

      do i = 1,dim
         do j = 1,dim
            c(j,i) = a(j,1)*b(i,1)
         enddo
      enddo
      do k = 2,dim
         do i = 1,dim
            do j = 1,dim
               c(j,i) = c(j,i)+a(j,k)*b(i,k)
            enddo
         enddo
      enddo

      return
      end

C ----------------------------------------------------------------------
C Library subroutine qqtxdd
C
C Multiplication of the transpose of a real quadratic matrix with a
C real quadratic matrix:
C     a(k,j)*b(k,i) = c(j,i)
C-----------------------------------------------------------------------

      subroutine qqtxdd (a,b,c,dim)

      implicit none

      integer dim,i,j,k
      real*8  a(dim,dim),b(dim,dim),c(dim,dim)

      do i = 1,dim
         do j = 1,dim
            c(j,i) = a(1,j)*b(1,i)
            do k = 2,dim
               c(j,i) = c(j,i)+a(k,j)*b(k,i)
            enddo
         enddo
      enddo

      return
      end

C----------------------------------------------------------------------
C Library subroutine qdxxzz
C
C Multiplication of a complex quadratic matrix with a diagonal complex 
C matrix:
C     a(j,i)*b(i) = c(j,i)
C-----------------------------------------------------------------------

      subroutine qdxxzz (a,b,c,dim)

      implicit none

      integer dim,i,j
      complex*16 a(dim,dim),b(dim),c(dim,dim)

      do i = 1,dim
         do j = 1,dim
            c(j,i) = a(j,i)*b(i)
         enddo
      enddo

      return
      end

C----------------------------------------------------------------------
C Library subroutine qdxxdz
C
C Multiplication of a complex quadratic matrix with a diagonal real
C matrix:
C     a(j,i)*b(i) = c(j,i)
C-----------------------------------------------------------------------
 
      subroutine qdxxdz (a,b,c,dim)
 
      implicit none
 
      integer dim,i,j
      real*8  b(dim)
      complex*16 a(dim,dim),c(dim,dim)
 
      do i = 1,dim
         do j = 1,dim
            c(j,i) = a(j,i)* b(i)
         enddo
      enddo
 
      return
      end

C----------------------------------------------------------------------
C Library subroutine qdxxdd
C
C Multiplication of a real quadratic matrix with a diagonal real matrix:
C     a(j,i)*b(i) = c(j,i)
C-----------------------------------------------------------------------

      subroutine qdxxdd (a,b,c,dim)

      implicit none

      integer dim,i,j
      real*8  a(dim,dim),b(dim),c(dim,dim)

      do i = 1,dim
         do j = 1,dim
            c(j,i) = a(j,i)*b(i)
         enddo
      enddo

      return
      end

C----------------------------------------------------------------------
C Library subroutine qdxxdd1
C
C Multiplication of a real quadratic matrix with a diagonal real matrix:
C     a(j,i)*b(i) = c(j,i)
C
C NB phdim is physical (leading) dimension, dim is used dimension 
C-----------------------------------------------------------------------

      subroutine qdxxdd1 (a,b,c,phdim,dim)

      implicit none

      integer phdim,dim,i,j
      real*8  a(phdim,dim),b(dim),c(phdim,dim)

      do i = 1,dim
         do j = 1,dim
            c(j,i) = a(j,i)*b(i)
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine hhxtzzc
C
C Multiplication of a complex hermitian matrix with the transpose of a 
C complex hermitian matrix, where the two matrices commute:
C     a(j,k)*b(i,k) = c(j,i)
C
C NB The fact that the two matrices commute means that the result of
C    the multiplication is also hermitian
C-----------------------------------------------------------------------

      subroutine hhxtzzc (a,b,c,dim)

      implicit none

      integer    dim,i,j,k
      complex*16 a(dim,dim),b(dim,dim),c(dim,dim)

      do i = 1,dim
         do j = i,dim
            c(j,i) = a(j,1)*b(i,1)
         enddo
      enddo
      do k = 2,dim
         do i = 1,dim
            do j = i,dim
               c(j,i) = c(j,i)+a(j,k)*b(i,k)
            enddo
         enddo
      enddo
      do i = 1,dim
         do j = 1,i-1
            c(j,i) = dconjg(c(i,j))
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine mmaxzzh
C
C Multiplication of the adjoint of a complex rectangular matrix with
C a rectangular complex matrix, where the result is a hermitian matrix
C     dconjg(a(k,j))*b(k,i) = c(j,i)
C
C NB this routine can be used for the overlap of two sets of "spfs" in
C    the same basis
C-----------------------------------------------------------------------

      subroutine mmaxzzh (a,b,c,dim1,dim2)

      implicit none

      integer     dim1,dim2,i,j,k
      complex*16  a(dim1,dim2),b(dim1,dim2),c(dim2,dim2)

      do i = 1,dim2
         do j = i,dim2
            c(j,i) = dconjg(a(1,j))*b(1,i)
            do k = 2,dim1
               c(j,i) = c(j,i)+dconjg(a(k,j))*b(k,i)
            enddo
         enddo
      enddo
C
C now fill in other half of hermitian matrix
C
      do i=1,dim2
         c(i,i)=dble(c(i,i))
      enddo
      do i=1,dim2
         do j=1,i-1
            c(j,i)=dconjg(c(i,j))
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine mmaxzzh1
C
C Multiplication of the adjoint of a complex rectangular matrix with
C a rectangular complex matrix, where the result is an anti- hermitian 
C matrix
C     dconjg(a(k,j))*b(k,i) = c(j,i)
C
C-----------------------------------------------------------------------

      subroutine mmaxzzh1 (a,b,c,dim1,dim2)

      implicit none

      integer     dim1,dim2,i,j,k
      real*8      x
      complex*16  a(dim1,dim2),b(dim1,dim2),c(dim2,dim2)

      do i = 1,dim2
         do j = i,dim2
            c(j,i) = dconjg(a(1,j))*b(1,i)
            do k = 2,dim1
               c(j,i) = c(j,i)+dconjg(a(k,j))*b(k,i)
            enddo
         enddo
      enddo
C
C now fill in other half of anti-hermitian matrix
C
      do i=1,dim2
         x=dimag(c(i,i))
         c(i,i)=dcmplx(0.0d0,x)
      enddo
      do i=1,dim2
         do j=1,i-1
            c(j,i)=-dconjg(c(i,j))
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine mmtxzzs
C
C Multiplication of a transposed complex rectangular matrix with
C a rectangular complex matrix, where the result is a symmetric matrix
C     a(k,j)*b(k,i) = c(j,i)
C-----------------------------------------------------------------------

      subroutine mmtxzzs (a,b,c,dim1,dim2)

      implicit none

      integer     dim1,dim2,i,j,k
      complex*16  a(dim1,dim2),b(dim1,dim2),c(dim2,dim2)

      do i = 1,dim2
         do j = i,dim2
            c(j,i) = a(1,j)*b(1,i)
            do k = 2,dim1
               c(j,i) = c(j,i)+a(k,j)*b(k,i)
            enddo
         enddo
      enddo
C
C now fill in other half of symmetric matrix
C
      do i=1,dim2
         do j=1,i-1
            c(j,i)=c(i,j)
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine qmxxzz
C
C Multiplication of a complex quadratic matrix with a complex 
C rectangular matrix:
C     a(j,k)*b(k,i) = c(j,i)
C-----------------------------------------------------------------------

      subroutine qmxxzz (a,b,c,dim1,dim2)

      implicit none

      integer     dim1,dim2,i,j,k
      complex*16  a(dim1,dim1),b(dim1,dim2),c(dim1,dim2)

      do i = 1,dim2
         do j = 1,dim1
            c(j,i) = a(j,1)*b(1,i)
         enddo
      enddo
      do i = 1,dim2
         do k = 2,dim1
            do j = 1,dim1
               c(j,i) = c(j,i)+a(j,k)*b(k,i)
            enddo
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine qqxxdd1
C
C Multiplication of two real quadratic matrices:
C     a(j,k)*b(k,i) = c(j,i)
C
C NB phdim is physical (leading) dimension, dim is used dimension 
C-----------------------------------------------------------------------

      subroutine qqxxdd1 (a,b,c,phdim,dim)

      implicit none

      integer phdim,dim,i,j,k
      real*8  a(phdim,dim),b(phdim,dim),c(phdim,dim)

      do i = 1,dim
         do j = 1,dim
            c(j,i) = a(j,1)*b(1,i)
         enddo
      enddo
      do i = 1,dim
         do k = 2,dim
            do j = 1,dim
               c(j,i) = c(j,i)+a(j,k)*b(k,i)
            enddo
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine mmxxzz1
C
C Multiplication of two complex rectangular matrices:
C     a(j,k)*b(k,i) = c(j,i)
C
C NB phdims are physical dimensions, dims are used dimensions
C-----------------------------------------------------------------------

C     subroutine mmxxzz1 (phdim1,phdim2,phdim3,dim1,dim2,dim3,a,b,c)

C     implicit none

C     integer     phdim1,phdim2,phdim3,dim1,dim2,dim3,i,j,k
C     complex*16  a(phdim1,phdim2),b(phdim2,phdim3),c(phdim1,phdim3)

C     do i = 1,dim3
C        do j = i,dim1
C           c(j,i) = a(j,1)*b(1,i)
C        enddo
C     enddo
C     do i = 1,dim3
C        do k = 2,dim2
C           do j = 1,dim1
C              c(j,i) = c(j,i)+a(j,k)*b(k,i)
C           enddo
C        enddo
C     enddo

C     return
C     end

C-----------------------------------------------------------------------
C Library subroutine mmaxzz
C
C Multiplication of the adjoint of a complex rectangular matrix with
C a rectangular complex matrix
C     dconjg(a(k,j))*b(k,i) = c(j,i)
C
C NB this routine can be used for the overlap of two sets of vectors in
C    different spf bases
C-----------------------------------------------------------------------
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C If you do some changes in mmaxzz, please do them in dwtmmaxzz too!
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mmaxzz (a,b,c,dim1,dim2,dim3)

      implicit none

      integer     dim1,dim2,dim3,i,j,k
      complex*16  a(dim1,dim2),b(dim1,dim3),c(dim2,dim3)

      do i = 1,dim3
         do j = 1,dim2
            c(j,i) = dconjg(a(1,j))*b(1,i)
            do k = 2,dim1
               c(j,i) = c(j,i)+dconjg(a(k,j))*b(k,i)
            enddo
         enddo
      enddo

      return
      end


C-----------------------------------------------------------------------
C Library subroutine mmtczz
C
C Multiplication of the transpose of a complex rectangular matrix with
C the complex conjugate of a rectangular complex matrix
C     a(k,j)*dconjg(b(k,i)) = c(j,i)
C
C-----------------------------------------------------------------------
      subroutine mmtczz (a,b,c,dim1,dim2,dim3)

      implicit none

      integer     dim1,dim2,dim3,i,j,k
      complex*16  a(dim1,dim2),b(dim1,dim3),c(dim2,dim3)

      do i = 1,dim3
         do j = 1,dim2
            c(j,i) = a(1,j)*dconjg(b(1,i))
            do k = 2,dim1
               c(j,i) = c(j,i)+a(k,j)*dconjg(b(k,i))
            enddo
         enddo
      enddo

      return
      end


C-----------------------------------------------------------------------
C Library subroutine mmaczz
C
C Multiplication of the adjoint of a complex rectangular matrix with
C a adjoint of a rectangular complex matrix
C     dconjg(a(k,j))*dconjg(b(k,i)) = c(j,i)
C
C NB this routine can be used for the overlap of two sets of vectors in
C    different spf bases
C-----------------------------------------------------------------------
      subroutine mmaczz (a,b,c,dim1,dim2,dim3)

      implicit none

      integer     dim1,dim2,dim3,i,j,k
      complex*16  a(dim1,dim2),b(dim1,dim3),c(dim2,dim3)

      do i = 1,dim3
         do j = 1,dim2
            c(j,i) = dconjg(a(1,j))*dconjg(b(1,i))
            do k = 2,dim1
               c(j,i) = c(j,i)+dconjg(a(k,j))*dconjg(b(k,i))
            enddo
         enddo
      enddo

      return
      end


C-----------------------------------------------------------------------
C Library subroutine mmcxzz
C
C Multiplication of the adjoint of a complex rectangular matrix with
C a rectangular complex matrix
C     dconjg(a(j,k))*b(k,i) = c(j,i)
C
C NB this routine can be used for the overlap of two sets of vectors in
C    different spf bases
C-----------------------------------------------------------------------
      subroutine mmcxzz (a,b,c,dim1,dim2,dim3)

      implicit none

      integer     dim1,dim2,dim3,i,j,k
      complex*16  a(dim1,dim2),b(dim2,dim3),c(dim1,dim3)

      do i = 1,dim3
         do j = 1,dim1
            c(j,i) = dconjg(a(j,1))*b(1,i)
            do k = 2,dim2
               c(j,i) = c(j,i)+dconjg(a(j,k))*b(k,i)
            enddo
         enddo
      enddo

      return
      end



C-----------------------------------------------------------------------
C Library subroutine qqcxzz
C
C Multiplication of the adjoint of a complex quadratic matrix with
C a quadratic complex matrix
C     dconjg(a(j,k))*b(k,i) = c(j,i)
C
C-----------------------------------------------------------------------
      subroutine qqcxzz (a,b,c,dim)

      implicit none

      integer     dim,i,j,k
      complex*16  a(dim,dim),b(dim,dim),c(dim,dim)

      do i = 1,dim
         do j = 1,dim
            c(j,i) = dconjg(a(j,1))*b(1,i)
            do k = 2,dim
               c(j,i) = c(j,i)+dconjg(a(j,k))*b(k,i)
            enddo
         enddo
      enddo

      return
      end





C-----------------------------------------------------------------------
C Library subroutine mmaxzz
C
C Multiplication of the adjoint of a complex rectangular matrix with
C a rectangular complex matrix
C     a(j,k)*dconjg(b(i,k)) = c(j,i)
C
C NB this routine can be used for the overlap of two sets of vectors in
C    different spf bases
C-----------------------------------------------------------------------

      subroutine mmxazz (a,b,c,dim1,dim2,dim3)

      implicit none

      integer     dim1,dim2,dim3,i,j,k
      complex*16  a(dim1,dim2),b(dim3,dim2),c(dim1,dim3)

      do i = 1,dim3
         do j = 1,dim1
            c(j,i) = 0.0d0
         enddo
      enddo

      do i = 1,dim3
         do j = 1,dim1
            do k = 1,dim2
               c(j,i) = c(j,i)+a(j,k)*dconjg(b(i,k))
            enddo
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine mmxtzza
C
C Multiplication of a complex rectangular matrix with the transpose of
C a rectangular complex matrix, the result of which is added to a
C further matrix.
C     a(j,k)*b(i,k) + c(j,i) = c(j,i)
C-----------------------------------------------------------------------

      subroutine mmxtzza (a,b,c,dim1,dim2,dim3)

      implicit none

      integer     dim1,dim2,dim3,i,j,k
      complex*16  a(dim1,dim2),b(dim3,dim2),c(dim1,dim3)

      do k = 1,dim2
         do i = 1,dim3
            do j = 1,dim1
               c(j,i) = c(j,i)+a(j,k)*b(i,k)
            enddo
         enddo
      enddo

      return
      end
      
      
C-----------------------------------------------------------------------
C Library subroutine mmcazza
C
C Multiplication of a conjugated complex rectangular matrix with the 
C adjoint of a quadratic complex matrix, the result of which is added 
C to a further rectangular matrix.
C     conjg(a(j,k))*conjg(b(i,k)) + c(j,i) = c(j,i)
C-----------------------------------------------------------------------
      subroutine mmcazza (a,b,c,dim1,dim2,dim3)

      implicit none

      integer     dim1,dim2,dim3,i,j,k
      complex*16  a(dim1,dim2),b(dim3,dim2),c(dim1,dim3)

      do k = 1,dim2
         do i = 1,dim3
            do j = 1,dim1
               c(j,i) = c(j,i)+conjg(a(j,k))*conjg(b(i,k))
            enddo
         enddo
      enddo

      return
      end


C-----------------------------------------------------------------------
C Library subroutine mqxxzza
C
C Multiplication of a complex rectangular matrix with a quadratic
C complex matrix, the result of which is added to a further rectangular
C matrix.
C     a(j,k)*b(k,i) + c(j,i) = c(j,i)
C-----------------------------------------------------------------------
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C If you do some changes in mqxxzza, please do them in
C dwtmqxxzzaandaddmxxzopar too!
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mqxxzza (a,b,c,dim1,dim2)

      implicit none

      integer     dim1,dim2,i,j,k
      complex*16  a(dim1,dim2),b(dim2,dim2),c(dim1,dim2)

      do i = 1,dim2
         do k = 1,dim2
            do j = 1,dim1
               c(j,i) = c(j,i)+a(j,k)*b(k,i)
            enddo
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine mqxtzd
C
C Multiplication of a rectangular complex matrix with the transpose of 
C a quadratic real matrix:
C     a(j,k)*b(i,k) = c(j,i)
C-----------------------------------------------------------------------

      subroutine mqxtzd (a,b,c,dim1,dim2)

      implicit none

      integer     dim1,dim2,i,j,k
      real*8      b(dim2,dim2)
      complex*16  a(dim1,dim2),c(dim1,dim2)

      do j=1,dim1
         do i=1,dim2
            c(j,i) = a(j,1)*b(i,1)
         enddo
      enddo
      do k=2,dim2
         do j=1,dim1
            do i=1,dim2
               c(j,i) = c(j,i)+a(j,k)*b(i,k)
            enddo
         enddo
      enddo

      return
      end

C----------------------------------------------------------------------
C Library subroutine qqtxzz
C
C Multiplication of the transpose of a complex quadratic matrix with a
C complex quadratic matrix:
C     a(k,j)*b(k,i) = c(j,i)
C-----------------------------------------------------------------------

      subroutine qqtxzz (a,b,c,dim)

      implicit none

      integer dim,i,j,k
      complex*16  a(dim,dim),b(dim,dim),c(dim,dim)

      do i = 1,dim
         do j = 1,dim
            c(j,i) = a(1,j)*b(1,i)
            do k = 2,dim
               c(j,i) = c(j,i)+a(k,j)*b(k,i)
            enddo
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine mqxazz
C
C Multiplication of a complex rectangular matrix with the adjoint of 
C a quadratic complex matrix:
C     a(j,k)*dconjg(b(i,k)) = c(j,i)
C-----------------------------------------------------------------------

      subroutine mqxazz (a,b,c,dim1,dim2)

      implicit none

      integer     dim1,dim2,i,j,k
      complex*16  a(dim1,dim2),b(dim2,dim2),c(dim1,dim2)

      do i = 1,dim2
         do j = 1,dim1
            c(j,i) = a(j,1)*dconjg(b(i,1))
         enddo
      enddo
      do k = 2,dim2
         do i = 1,dim2
            do j = 1,dim1
               c(j,i) = c(j,i)+a(j,k)*dconjg(b(i,k))
            enddo
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine qmxxdz
C
C Multiplication of a quadratic real matrix with a rectangular complex
C matrix:
C     a(j,k)*b(k,i) = c(j,i)
C-----------------------------------------------------------------------

      subroutine qmxxdz (a,b,c,dim1,dim2)

      implicit none

      integer     dim1,dim2,i,j,k
      real*8      a(dim1,dim1)
      complex*16  b(dim1,dim2),c(dim1,dim2)

      do i=1,dim2
         do j=1,dim1
            c(j,i) = a(j,1)*b(1,i)
         enddo
      enddo
      do i = 1,dim2
         do k = 2,dim1
            do j = 1,dim1
               c(j,i) = c(j,i)+a(j,k)*b(k,i)
            enddo
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine qmtxdz
C
C Multiplication of the transpose of a quadratic real matrix with a 
C rectangular complex matrix:
C     a(k,j)*b(k,i) = c(j,i)
C-----------------------------------------------------------------------

      subroutine qmtxdz (a,b,c,dim1,dim2)

      implicit none

      integer     dim1,dim2,i,j,k
      real*8      a(dim1,dim1)
      complex*16  b(dim1,dim2),c(dim1,dim2)

      do i = 1,dim2
         do j = 1,dim1
            c(j,i) = a(1,j)*b(1,i)
            do k = 2,dim1
               c(j,i) = c(j,i)+a(k,j)*b(k,i)
            enddo
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine qqxxzz
C
C Multiplication of two complex quadratic matrices:
C     a(j,k)*b(k,i) = c(j,i)
C-----------------------------------------------------------------------
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C If you do some changes in qqxxzz, please do them in dwtqqxxzz too!
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine qqxxzz (a,b,c,dim)

      implicit none

      integer dim,i,j,k
      complex*16   a(dim,dim),b(dim,dim),c(dim,dim)

      do i=1,dim
         do j=1,dim
            c(j,i) = a(j,1)*b(1,i)
         enddo
      enddo
      do i = 1,dim
         do k = 2,dim
            do j = 1,dim
               c(j,i) = c(j,i)+a(j,k)*b(k,i)
            enddo
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine mqxtzz
C
C Multiplication of a complex rectangular matrix with the transpose of
C a quadratic complex matrix:
C     a(j,k)*b(i,k) = c(j,i)
C-----------------------------------------------------------------------
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C If you do some changes in mqxtzz, please do them in dwtmqxtzzpar too!
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mqxtzz (a,b,c,dim1,dim2)

      implicit none

      integer     dim1,dim2,i,j,k
      complex*16  a(dim1,dim2),b(dim2,dim2),c(dim1,dim2)

      do i = 1,dim2
         do j = 1,dim1
            c(j,i) = a(j,1)*b(i,1)
         enddo
      enddo
      do k = 2,dim2
         do i = 1,dim2
            do j = 1,dim1
               c(j,i) = c(j,i)+a(j,k)*b(i,k)
            enddo
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine mqxxzz
C
C Multiplication of a complex rectangular matrix with a quadratic
C complex matrix:
C     a(j,k)*b(k,i) = c(j,i)
C-----------------------------------------------------------------------
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C If you do some changes in mqxxzz, please do them in dwtmqxxzz too!
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine mqxxzz (a,b,c,dim1,dim2)

      implicit none

      integer     dim1,dim2,i,j,k
      complex*16  a(dim1,dim2),b(dim2,dim2),c(dim1,dim2)

      do i = 1,dim2
         do j = 1,dim1
            c(j,i) = a(j,1)*b(1,i)
         enddo
      enddo
      do i = 1,dim2
         do k = 2,dim2
            do j = 1,dim1
               c(j,i) = c(j,i)+a(j,k)*b(k,i)
            enddo
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine mmtxzz
C
C Multiplication of a transposed complex rectangular matrix with
C a rectangular complex matrix
C     a(k,j)*b(k,i) = c(j,i)
C-----------------------------------------------------------------------

      subroutine mmtxzz (a,b,c,dim1,dim2,dim3)

      implicit none

      integer     dim1,dim2,dim3,i,j,k
      complex*16  a(dim1,dim2),b(dim1,dim3),c(dim2,dim3)

      do i = 1,dim3
         do j = 1,dim2
            c(j,i) = a(1,j)*b(1,i)
            do k = 2,dim1
               c(j,i) = c(j,i)+a(k,j)*b(k,i)
            enddo
         enddo
      enddo

      return
      end
C-----------------------------------------------------------------------
C Library subroutine mmtxzza
C
C Multiplication of a transposed complex rectangular matrix with
C a rectangular complex matrix
C     a(k,j)*b(k,i) + c(j,i) = c(j,i)
C-----------------------------------------------------------------------

      subroutine mmtxzza (a,b,c,dim1,dim2,dim3)

      implicit none

      integer     dim1,dim2,dim3,i,j,k
      complex*16  a(dim1,dim2),b(dim1,dim3),c(dim2,dim3)

      do i = 1,dim3
         do j = 1,dim2
            do k = 1,dim1
               c(j,i) = c(j,i)+a(k,j)*b(k,i)
            enddo
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine mmxczz
C
C Multiplication of a complex rectangular matrices with the complex
C conjugate of a rectangular matrix:
C     a(j,k)*dconjg(b(k,i)) = c(j,i)
C-----------------------------------------------------------------------

      subroutine mmxczz (a,b,c,dim1,dim2,dim3)

      implicit none

      integer     dim1,dim2,dim3,i,j,k
      complex*16  a(dim1,dim2),b(dim2,dim3),c(dim1,dim3)

      do i=1,dim3
         do j=1,dim1
            c(j,i) = a(j,1)*dconjg(b(1,i))
         enddo
      enddo
      do i = 1,dim3
         do k = 2,dim2
            do j = 1,dim1
               c(j,i) = c(j,i)+a(j,k)*dconjg(b(k,i))
            enddo
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine mqxczz
C
C Multiplication of a complex rectangular matrices with the complex
C conjugate of a quadratic matrix:
C     a(j,k)*dconjg(b(k,i)) = c(j,i)
C-----------------------------------------------------------------------

      subroutine mqxczz (a,b,c,dim1,dim2)

      implicit none

      integer     dim1,dim2,i,j,k
      complex*16  a(dim1,dim2),b(dim2,dim2),c(dim1,dim2)

      do i=1,dim2
         do j=1,dim1
            c(j,i) = a(j,1)*dconjg(b(1,i))
         enddo
      enddo
      do i = 1,dim2
         do k = 2,dim2
            do j = 1,dim1
               c(j,i) = c(j,i)+a(j,k)*dconjg(b(k,i))
            enddo
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C     Library subroutine qqaxzz                                         
C                                                                       
C Multiplication of the adjoint of a complex quadratic matrix with      
C a another quadratic complex matrix                                    
C     dconjg(a(k,j))*b(k,i) = c(j,i)                                    
C                                                                       
C-----------------------------------------------------------------------

      subroutine qqaxzz (a,b,c,dim)

      implicit none

      integer     dim,i,j,k
      complex*16  a(dim,dim),b(dim,dim),c(dim,dim)

      do i = 1,dim
         do j = 1,dim
            c(j,i) = dconjg(a(1,j))*b(1,i)
            do k = 2,dim
               c(j,i) = c(j,i)+dconjg(a(k,j))*b(k,i)
            enddo
         enddo
      enddo

      return
      end
