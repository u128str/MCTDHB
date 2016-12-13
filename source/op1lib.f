C-----------------------------------------------------------------------
C                    OP1LIB
C       SUBROUTINES ACTING ON A SINGLE OBJECT (LEVEL 1 ROUTINES)
C
C NOMENCLATURE:
C    Each of the following routines starts with:
C       sum2: the squares of elements are summed
C       norm: the norm (sqrt of sum of squares) is calculated
C       tr:   the trace is calculated
C       zero: all elements are set to zero
C       unit: sets up a unit matrix
C       cp:   input matrix/vector is copied to output matrix/vector.
C       init: initialises an object i.e. sets all elements to 1.
C       over: calculates an overlap matrix between two sets of vectors
C       tran: transposes input object
c       cut:  calculate a cut through the source object (the target has
c             one dimension less than the source object)
C    This is followed by four chracters:
C    Character 1 denotes the object on which the action is performed:
C       v: vector
C       q: quadratic matrix
C       m: general (rectangular) matrix
C       h: hermitian matrix
c       t: tensor of third order
C    Character 2 defines how the object is used:
C       t: in sum2, the product of the element with its transpose is
C          summed
C       a: in cp, the adjoint of the object is copied
C       x: this is a blank
C    Character 3 defines the data type of the object:
C       s: real single precision (real*4)
C       d: real double precision (real*8)
C       c: complex single precision (complex*8)
C       z: complex double precision (complex*16)
C       i: integer
C       l: logical
C     Further characters, if present, give more informaion:
C       1: the physical dimensions of the matrices differs from those
C          used.
C       n: (= 1 or 2 or ...) in a cut the dimension to be held constant
c     
C Contents:
C In the following list of available subroutines, objects on the LHS
C of the definition are input, that on the RHS output. The usual
C summation convention is used i.e. a sum is made over repeated indices
C on the LHS
C
C    sum2vxd (vec,sum2,dim)
C        Definition: vec(i)*vec(i) = sum2
C        Dimensions: vec(dim),sum2
C
C    sum2vxz (vec,sum2,dim)
C        Definition: dconjg(vec(i))*vec(i) = sum2
C        Dimensions: vec(dim),sum2
C
CC   sum2qxd (mat,sum2,dim)
CC       Definition: mat(i,j)*mat(i,j) = sum2
CC       Dimensions: mat(dim,dim),sum2
C
C    sum2mxd (mat,sum2,dim1,dim2)
C        Definition: mat(i,j)*mat(i,j) = sum2
C        Dimensions: mat(dim1,dim2),sum2
C
C    sum2qtz (mat,sum2,dim)
C        Definition: mat(i,j)*mat(j,i) = sum2
C        Dimensions: mat(dim,dim),sum2
C
C    normvxd (vec,norm,dim)
C        Definition: vec(i)*vec(i) = norm**2
C        Dimensions: vec(dim),norm
C
C    normvxz (vec,norm,dim)
C        Definition: dconjg(vec(i))*vec(i) = norm**2
C        Dimensions: vec(dim),norm
C
C    trhxz (mat,trace,dim)
C        Definition: dble(mat(i,i)) = trace
C        Dimensions: mat(dim,dim),trace
C
C    trqxd (mat,trace,dim)
C        Definition: mat(i,i) = trace
C        Dimensions: mat(dim,dim),trace
C
C    trqxz (mat,trace,dim)
C        Definition: mat(i,i) = trace
C        Dimensions: mat(dim,dim),trace
C
C    trmmcxzz(a,b,trace,dim1,dim2)
C        Definiton: conjg(a(i,j))*b(j,i) = trace
C        Dimension: a(dim1,dim2), b(dim2,dim1)
C
C    trmmaxzz(a,b,trace,dim1,dim2)
C      Definiton: conjg(a(j,i))*b(j,i) = trace
C      Dimension: a(dim1,dim2), b(dim1,dim2)
C
C    trmmaxzza(a,b,scal,dim1,dim2)
C        Definiton: conjg(a(i,j))*b(j,i)) = scal
C        Dimension: a(dim1,dim2), b(dim1,dim2)
C
C    trtxz(ten,vec,dim1,dim2)
C        Definiton: ten(i,j,j)) = vec(i)
C        Dimension: ten(dim1,dim2,dim2), vec(dim1)
C
C    zeromxz (mat,dim1,dim2)
C        Definition: mat(i,j) = (0.0d0,0.0d0)
C        Dimensions: mat(dim1,dim2)
C
C    zeromxd (mat,dim1,dim2)
C        Definition: mat(i,j) = 0.0d0
C        Dimensions: mat(dim1,dim2)
C
CC   zeromxs (mat,dim1,dim2)
CC       Definition: mat(i,j) = 0.0
CC       Dimensions: mat(dim1,dim2)
C
C    zeromxi (mat,dim1,dim2)
C        Definition: mat(i,j) = 0
C        Dimensions: mat(dim1,dim2)
C
C    zeromxl (mat,dim1,dim2)
C        Definition: mat(i,j) = .false.
C        Dimensions: mat(dim1,dim2)
C
C    zerovxz (vec,dim)
C        Definition: vec(i) = (0.0d0,0.0d0)
C        Dimensions: vec(dim)
C
C    zerovxd (vec,dim)
C        Definition: vec(i) = 0.0d0
C        Dimensions: vec(dim)
C
C    zerovxs (vec,dim)
C        Definition: vec(i) = 0.0d0
C        Dimensions: vec(dim)
C
C    zerovxi (vec,dim)
C        Definition: vec(i) = 0
C        Dimensions: vec(dim)
C
C    zerovxl (vec,dim)
C        Definition: vec(i) = .false.
C        Dimensions: vec(dim)
C
C    unitqxz (mat,dim)
C        Definition: mat(i,j) = d(i,j)   (d(i,j) is the kronecker delta)
C        Dimensions: mat(dim,dim)
C
C    unitqxd (mat,dim)
C        Definition: mat(i,j) = d(i,j)   (d(i,j) is the kronecker delta)
C        Dimensions: mat(dim,dim)
C
C    cpqxd (a,c,dim)
C        Definition: c(j,i) = a(j,i)
C        Dimensions: a(dim,dim),c(dim,dim)
C
C    cpqxd1 (a,c,phdim,dim)
C        Definition: c(j,i) = a(j,i) ; 1 <= i,j <= dim
C        Dimensions: a(phdim,dim),c(phdim,dim)
C
C    cpqxz (a,c,dim)
C        Definition: c(j,i) = a(j,i)
C        Dimensions: a(dim,dim),c(dim,dim)
C
C    cpqxdz (a,c,dim)
C        Definition: c(j,i) = a(j,i)
C        Dimensions: a(dim,dim),c(dim,dim)
C
C    cpqaz (a,c,dim)
C       Definition: c(i,j) = dconj(a(j,i))
C       Dimensions: a(dim,dim),c(dim,dim)
C
C    cpqtz (a,c,dim)
C        Definition: c(i,j) = a(j,i)
C        Dimensions: a(dim,dim),c(dim,dim)
C
C    cpmaz (a,c,dim1,dim2) 
C        Definition: c(i,j) = dconj(a(j,i))
C        Dimensions: a(dim1,dim2),c(dim2,dim1)
C
C    cpmtz (a,c,dim1,dim2)
C        Definition: c(i,j) = a(j,i)
C        Dimensions: a(dim1,dim2),c(dim2,dim1)
C
C    cpmxz (a,c,dim1,dim2)
C        Definition: c(j,i) = a(j,i)
C        Dimensions: a(dim1,dim2),c(dim1,dim2)
C
C    cpmxd (a,c,dim1,dim2)
C        Definition: c(j,i) = a(j,i)
C        Dimensions: a(dim1,dim2),c(dim1,dim2)
C
C    cpmxi (a,c,dim1,dim2)
C        Definition: c(j,i) = a(j,i)
C        Dimensions: a(dim1,dim2),c(dim1,dim2)
C
C    cpvxd (v,w,dim)
C        Definition: w(i) = v(i)
C        Dimensions: v(dim),w(dim)
C
C    cpvxz (v,w,dim)
C        Definition: w(i) = v(i)
C        Dimensions: v(dim),w(dim)
C
C    cpvcz (v,w,dim)
C        Definition: w(i) = dconjg(v(i))
C        Dimensions: v(dim), w(dim)
C
C    cpvxdz (v,w,dim)
C        Definition: w(i) = v(i)
C        Dimensions: v(dim),w(dim)
C
C    cpvxzd (v,w,dim)
C        Definition: w(i) = v(i)
C        Dimensions: v(dim),w(dim)
C
C    cpvxi (v,w,dim)
C        Definition: w(i) = v(i)
C        Dimensions: v(dim),w(dim)
C
C    cpvxi2(v,w,dim)
C        Definition: w(i) = v(i) &&  v(i) = w(i)
C        Dimensions: v(dim),w(dim)
C
C    cpvxl (v,w,dim)
C        Definition: w(i) = v(i)
C        Dimensions: v(dim),w(dim)
C
C    initvxz (vec,dim)
C        Definition: vec(i) = (1.0d0,0.0d0)
C        Dimensions: vec(dim)
C
C    initvxd (vec,dim)
C        Definition: vec(i) = 1.0d0
C        Dimensions: vec(dim)
C
C    initvxi (vec,dim)
C        Definition: vec(i) = 1
C        Dimensions: vec(dim)
C
C    initvxl (vec,dim)
C        Definition: vec(i) = .true.
C        Dimensions: vec(dim)
C
C    initmxz (a,dim1,dim2)
C        Definition: a(j,i) = (1.0d0,0.0d0)
C        Dimensions: a(dim1,dim2)
C
C    initmxd (a,dim1,dim2)
C        Definition: a(j,i) = (1.0d0,0.0d0)
C        Dimensions: a(dim1,dim2)
C
C    initmxl (a,dim1,dim2)
C        Definition: a(j,i) = .true.
C        Dimensions: a(dim1,dim2)
C
C    overmxz (a,c,dim1,dim2)
C        Definition: dconjg(a(k,j))*a(k,i) = c(j,i)
C        Dimensions: a(dim1,dim2),c(dim2,dim2)
C
C    overmcz (a,c,dim1,dim2)
C        Definition: a(k,j)*a(k,i) = c(j,i)
C        Dimensions: a(dim1,dim2),c(dim2,dim2)
C
C    tranmxz (a,c,dim1,dim2)
C        Definition: a(k,j) = c(j,k)
C        Dimensions: a(dim1,dim2),c(dim2,dim1)
C
C    tranmxd (a,c,dim1,dim2)
C        Definition: a(k,j) = c(j,k)
C        Dimensions: a(dim1,dim2),c(dim2,dim1)
C
C    tranqxz (a,dim)
C        Definition: a(k,j) = a(j,k)
C        Dimensions: a(dim,dim)
C
C    tranqxd (a,dim)
C        Definition: a(k,j) = a(j,k)
C        Dimensions: a(dim,dim)
C
C    cpvxz_s (v,w,dim1,dim2,index1,index2)
C        Definition: index1(j)=i
C                    index2(k)=i
C                    v(i) = w(i)
C        Dimensions: v(dim1),w(dim2),index1(dim1),index2(dim2)
c     
c    cuttxd2 (t,m,dim1,dim2,dim3,jcut)
c        Definition: m(i,k)=t(i,jcut,k)
c        Dimensions: t(dim1,dim2,dim3),m(dim1,dim3)
C
C-----------------------------------------------------------------------


c-----------------------------------------------------------------------
c Library subroutine sum2vxd
c
C sums the squares of the elements of a real vector
C        vec(i)*vec(i) = sum2
C i.e. the square of the norm of the vector
c-----------------------------------------------------------------------

      subroutine sum2vxd(vec,sum2,dim)

      implicit none

      integer i,dim
      real*8  vec(dim),sum2
      
      sum2=0.
      do i=1,dim
         sum2=sum2+vec(i)*vec(i)
      enddo

      return
      end

c-----------------------------------------------------------------------
c Library subroutine sum2vxz
c
C sums the squares of the elements of a complex vector
C        dconjg(vec(i))*vec(i) = sum2
C i.e. the square of the norm of the vector
c-----------------------------------------------------------------------

      subroutine sum2vxz(vec,sum2,dim)

      implicit none

      integer     i,dim
      complex*16  vec(dim),sum2
      
      sum2=(0.0d0,0.0d0)
      do i=1,dim
         sum2=sum2+dconjg(vec(i))*vec(i)
      enddo

      return
      end

c-----------------------------------------------------------------------
c Library subroutine sum2qxd
c
C sums the squares of the elements of a real quadratic matrix
C        mat(i,j)*mat(i,j) = sum2
c-----------------------------------------------------------------------
     
C     subroutine sum2qxd(mat,sum2,dim)
C     
C     implicit none

C     integer i,j,dim
C     real*8  mat(dim,dim),sum2

C     sum2=0.
C     do j=1,dim
C        do i=1,dim
C           sum2=sum2+mat(i,j)*mat(i,j)
C        enddo
C     enddo

C     return
C     end

C----------------------------------------------------------------------- 
C Library subroutine sum2mxd
C
C sums the squares of the elements of a real rectangular matrix
C        mat(i,j)*mat(i,j) = sum2
C-----------------------------------------------------------------------
    
      subroutine sum2mxd(mat,sum2,dim1,dim2)
      
      implicit none

      integer i,j,dim1,dim2
      real*8  mat(dim1,dim2),sum2

      sum2=0.
      do j=1,dim2
         do i=1,dim1
            sum2=sum2+mat(i,j)*mat(i,j)
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine sum2qtz
C
C sums the products of the elements of a real quadratic matrix
C with its transpose
C        mat(i,j)*mat(j,i) = sum2
C-----------------------------------------------------------------------
     
      subroutine sum2qtz(mat,sum2,dim)
      
      implicit none

      integer     i,j,dim
      complex*16  mat(dim,dim),sum2

      sum2=0.
      do j=1,dim
         do i=1,dim
            sum2=sum2+mat(i,j)*mat(j,i)
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine normvxd
C
C calculates the norm of a real vector
C        vec(i)*vec(i) = norm**2
C-----------------------------------------------------------------------

      subroutine normvxd(vec,norm,dim)
     
      implicit none

      integer i,dim
      real*8  vec(dim),norm
     
      norm=vec(1)*vec(1)
      do i=2,dim
         norm=norm+vec(i)*vec(i)
      enddo
      norm=sqrt(norm)

      return
      end

C-----------------------------------------------------------------------
c Library subroutine normvxz
c
C calculates the norm of a complex vector
C        dconjg(vec(i))*vec(i) = norm**2
c-----------------------------------------------------------------------

      subroutine normvxz(vec,norm,dim)
      
      implicit none

      integer     i,dim
      real*8      norm
      complex*16  vec(dim)
      
      norm=0.
      do i=1,dim
         norm=norm+dconjg(vec(i))*vec(i)
      enddo
      norm=sqrt(norm)

      return
      end


C-----------------------------------------------------------------------
C Library subroutine trhxz
C
C calculates the trace of a hermitian matrix
C      dble(mat(i,i)) = trace
C
C NB diagonal elements of complex matrix are real due to hermiticity
C-----------------------------------------------------------------------

      subroutine trhxz (mat,trace,dim)
      
      implicit none

      integer    dim,i
      real*8     trace
      complex*16 mat(dim,dim)
      
      trace=dble(mat(1,1))
      do i=2,dim
         trace=trace+dble(mat(i,i))
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine trqxd
C
C calculates the trace of a quadratic matrix
C      mat(i,i) = trace
C-----------------------------------------------------------------------

      subroutine trqxd (mat,trace,dim)
      
      implicit none

      integer    dim,i
      real*8     mat(dim,dim),trace
      
      trace=mat(1,1)
      do i=2,dim
         trace=trace+mat(i,i)
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine trqxz
C
C calculates the trace of a quadratic matrix
C      mat(i,i) = trace
C-----------------------------------------------------------------------

      subroutine trqxz (mat,trace,dim)
      
      implicit none

      integer    dim,i
      complex*16 mat(dim,dim),trace
      
      trace=mat(1,1)
      do i=2,dim
         trace=trace+mat(i,i)
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine trmmcxzz
C
C calculates the trace of the hermitian product of two complex matrices
C
C      Definiton: conjg(a(i,j))*b(j,i) = trace
C      Dimension: a(dim1,dim2), b(dim2,dim1)
C
C NB diagonal elements of complex matrix are real due to hermiticity
C-----------------------------------------------------------------------

      subroutine trmmcxzz(a,b,trace,dim1,dim2)

      implicit none

      integer    dim1, dim2, i, j
      complex*16 a(dim1,dim2), b(dim2,dim1), trace

      trace = conjg(a(1,1))*b(1,1)
      do j=2,dim2
         trace = trace + conjg(a(1,j))*b(j,1)
      enddo
      do i=2,dim1
         do j=1,dim2
            trace = trace + conjg(a(i,j))*b(j,i) 
         enddo
      enddo

      return
      end


C-----------------------------------------------------------------------
C Library subroutine trmmcxzz
C
C calculates the trace of the hermitian product of two complex matrices
C
C      Definiton: conjg(a(j,i))*b(j,i) = trace
C      Dimension: a(dim1,dim2), b(dim1,dim2
C
C-----------------------------------------------------------------------
      subroutine trmmaxzz(a,b,trace,dim1,dim2)

      implicit none

      integer    dim1, dim2, i, j
      complex*16 a(dim1,dim2), b(dim1,dim2), trace

      trace = conjg(a(1,1))*b(1,1)
      do j=2,dim1
         trace = trace + conjg(a(j,1))*b(j,1)
      enddo
      do i=2,dim2
         do j=1,dim1
            trace = trace + conjg(a(j,i))*b(j,i) 
         enddo
      enddo

      return
      end



C-----------------------------------------------------------------------
C Library subroutine trmmaxzza
C
C calculates the trace of the product of two complex matrices and adds 
C it to scal
C
C         Definiton: conjg(a(j,i))*b(j,i)) + scal = scal
C         Dimension: a(dim1,dim2), b(dim1,dim2)
C
C-----------------------------------------------------------------------

      subroutine trmmaxzza(a,b,scal,dim1,dim2)

      implicit none

      integer    dim1, dim2, i, j
      complex*16 a(dim1,dim2), b(dim1,dim2), scal

      do j=1,dim2
         do i=1,dim1
            scal = scal + dconjg(a(i,j)) * b(i,j) 
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine trtxz
C
C calculates the trace of tensor, reducing it to a vector.
C
C         Definiton: ten(i,j,j) = vec(i)
C         Dimension: ten(dim1,dim2,dim2), vec(dim1)
C
C-----------------------------------------------------------------------

      subroutine trtxz(ten,vec,dim1,dim2)

      implicit none

      integer    dim1, dim2, i, j
      complex*16 ten(dim1,dim2,dim2), vec(dim1)

      call zerovxz(vec,dim1)

      do j=1,dim2
         do i=1,dim1
            vec(i) = vec(i) + ten(i,j,j) 
         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
c Library subroutine zeromxz
c
C makes all the elements of a complex rectangular matrix zero
C      mat(i,j) = (0.0d0,0.0d0)
c-----------------------------------------------------------------------

      subroutine zeromxz(mat,dim1,dim2)

      implicit none

      integer    dim1,dim2,i,j
      complex*16 mat(dim1,dim2)

      do j=1,dim2
         do i=1,dim1
            mat(i,j)=(0.0d0,0.0d0)
         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
c Library subroutine zeromxd
c
C makes all the elements of a real*8 rectangular matrix zero
C      mat(i,j) = 0.0d0
c-----------------------------------------------------------------------

      subroutine zeromxd(mat,dim1,dim2)

      implicit none

      integer    dim1,dim2,i,j
      real*8     mat(dim1,dim2)

      do j=1,dim2
         do i=1,dim1
            mat(i,j)=0.0d0
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine zeromxs
C
C makes all the elements of a real*4 rectangular matrix zero
C      mat(i,j) = 0.0
C-----------------------------------------------------------------------

C     subroutine zeromxs(mat,dim1,dim2)

C     implicit none

C     integer    dim1,dim2,i,j
C     real*4     mat(dim1,dim2)

C     do j=1,dim2
C        do i=1,dim1
C           mat(i,j)=0.0
C        enddo
C     enddo

C     return
C     end

c-----------------------------------------------------------------------
c Library subroutine zeromxi
c
C makes all the elements of a rectangular integer matrix zero
C      mat(i,j) = 0
c-----------------------------------------------------------------------

      subroutine zeromxi(mat,dim1,dim2)

      implicit none

      integer    dim1,dim2,i,j
      integer    mat(dim1,dim2)

      do j=1,dim2
         do i=1,dim1
            mat(i,j)=0
         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
c Library subroutine zeromxl
c
C makes all the elements of a rectangular matrix of logicals .false.
C     mat(i,j) = .false.
c-----------------------------------------------------------------------

      subroutine zeromxl(mat,dim1,dim2)

      implicit none

      integer    dim1,dim2,i,j
      logical    mat(dim1,dim2)

      do j=1,dim2
         do i=1,dim1
            mat(i,j)=.false.
         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
c Library subroutine zerovxz
c
C makes all the elements of a complex vector zero
C     vec(i)=(0,0d0,0.0d0)
c-----------------------------------------------------------------------

      subroutine zerovxz(vec,dim)

      implicit none

      integer    dim,i
      complex*16 vec(dim)

      do i=1,dim
         vec(i)=(0.0d0,0.0d0)
      enddo

      return
      end

c-----------------------------------------------------------------------
c Library subroutine zerovxd
c
C makes all the elements of a real*8 vector zero
C     vec(i)=0.0d0
c-----------------------------------------------------------------------
      subroutine zerovxd(vec,dim)

      implicit none

      integer    dim,i
      real*8     vec(dim)

      do i=1,dim
         vec(i)=0.0d0
      enddo

      return
      end

c-----------------------------------------------------------------------
c Library subroutine zerovxs
c
C makes all the elements of a real*4 vector zero
C     vec(i)=0.0d0
c-----------------------------------------------------------------------
      subroutine zerovxs(vec,dim)

      implicit none

      integer    dim,i
      real*4     vec(dim)

      do i=1,dim
         vec(i)=0.0d0
      enddo

      return
      end

c-----------------------------------------------------------------------
c Library subroutine zerovxi
c
C makes all the elements of an integer vector zero
C     vec(i)=0
c-----------------------------------------------------------------------
      subroutine zerovxi(vec,dim)

      implicit none

      integer    dim,i
      integer    vec(dim)

      do i=1,dim
         vec(i)=0
      enddo

      return
      end

c-----------------------------------------------------------------------
c Library subroutine zerovxl
c
C makes all the elements of a logical vector .false.
C     vec(i)=.false.
c-----------------------------------------------------------------------
      subroutine zerovxl(vec,dim)

      implicit none

      integer    dim,i
      logical    vec(dim)

      do i=1,dim
         vec(i)=.false.
      enddo

      return
      end

c-----------------------------------------------------------------------
c Library subroutine unitqxz
c
C sets up a unit matrix
C     mat(i,j) = 0, if i.ne.j
C     mat(i,j) = 1, if i.eq.j
c-----------------------------------------------------------------------
      subroutine unitqxz(mat,dim)

      implicit none

      integer    dim,i,j
      complex*16 mat(dim,dim)

      do j=1,dim
         do i=1,dim
            mat(i,j)=(0.0d0,0.0d0)
         enddo
      enddo
      do i=1,dim
         mat(i,i)=(1.0d0,0.0d0)
      enddo

      return
      end

c-----------------------------------------------------------------------
c Library subroutine unitqxd
c
C sets up a unit matrix
C     mat(i,j) = 0, if i.ne.j
C     mat(i,j) = 1, if i.eq.j
c-----------------------------------------------------------------------
      subroutine unitqxd(mat,dim)

      implicit none

      integer    dim,i,j
      real*8 mat(dim,dim)

      do j=1,dim
         do i=1,dim
            mat(i,j)=(0.0d0,0.0d0)
         enddo
      enddo
      do i=1,dim
         mat(i,i)=(1.0d0,0.0d0)
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine cpqxd
C
C copies a real quadratic matrix to another real quadratic matrix
C     c(i,j) = a(i,j)
C-----------------------------------------------------------------------

      subroutine cpqxd (a,c,dim)

      implicit none

      integer dim,i,j
      real*8  a(dim,dim),c(dim,dim)

      do i = 1,dim
         do j = 1,dim
            c(j,i) = a(j,i)
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine cpqxd1
C
C copies a real quadratic matrix to another real quadratic matrix
C     c(i,j) = a(i,j)
C
C NB phdim is physical (leading) dimension, dim is used dimension 
C-----------------------------------------------------------------------

      subroutine cpqxd1 (a,c,phdim,dim)

      implicit none

      integer phdim,dim,i,j
      real*8  a(phdim,dim),c(phdim,dim)

      do i = 1,dim
         do j = 1,dim
            c(j,i) = a(j,i)
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine cpqxz
C
C copies a complex quadratic matrix to another complex quadratic matrix
C     c(i,j) = a(i,j)
C-----------------------------------------------------------------------

      subroutine cpqxz (a,c,dim)

      implicit none

      integer dim,i,j
      complex*16  a(dim,dim),c(dim,dim)

      do i = 1,dim
         do j = 1,dim
            c(j,i) = a(j,i)
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine cpqxdz
C
C copies a real quadratic matrix to a complex quadratic matrix
C     c(i,j) = a(i,j)
C-----------------------------------------------------------------------

      subroutine cpqxdz (a,c,dim)

      implicit none

      integer dim,i,j
      real*8  a(dim,dim)
      complex*16  c(dim,dim)

      do i = 1,dim
         do j = 1,dim
            c(j,i) = a(j,i)
         enddo
      enddo

      return
      end

C ----------------------------------------------------------------------
C Library subroutine cpqaz
C
C copies the adjoint of a complex quadratic matrix to a complex 
C quadratic matrix
C     c(j,i) = dconjg(a(i,j)) 
C-----------------------------------------------------------------------

      subroutine cpqaz (a,c,dim)

      implicit none

      integer     dim,i,j
      complex*16  a(dim,dim),c(dim,dim)

      do i = 1,dim
         do j = 1,dim
            c(j,i) = dconjg(a(i,j))
         enddo
      enddo

      return
      end

C ----------------------------------------------------------------------
C Library subroutine cpqtz
C
C copies the transpose of a complex quadratic matrix to a complex 
C quadratic matrix
C     c(j,i) = a(i,j) 
C-----------------------------------------------------------------------

      subroutine cpqtz (a,c,dim)

      implicit none
      
      integer     dim,i,j
      complex*16  a(dim,dim),c(dim,dim)

      do i = 1,dim
         do j = 1,dim
            c(j,i) = a(i,j)
         enddo
      enddo
      
      return
      end


C----------------------------------------------------------------------
C Library subroutine cpmtz
C
C copies the transpose of a complex rectangular atrix to a complex
C rectangular matrix
C     c(j,i) = a(i,j)
C-----------------------------------------------------------------------

      subroutine cpmtz (a,c,dim1,dim2)

      implicit none

      integer     dim1,dim2,i,j
      complex*16  a(dim1,dim2),c(dim2,dim1)

      do i = 1,dim1
         do j = 1,dim2
            c(j,i) = a(i,j)
         enddo
      enddo

      return
      end


C----------------------------------------------------------------------
C Library subroutine cpmaz
C
C copies the adjoint of a rectangular quadratic matrix to a complex
C rectangular matrix
C     c(j,i) = dconjg(a(i,j))
C-----------------------------------------------------------------------

      subroutine cpmaz (a,c,dim1,dim2)

      implicit none

      integer     dim1,dim2,i,j
      complex*16  a(dim1,dim2),c(dim2,dim1)

      do i = 1,dim1
         do j = 1,dim2
            c(j,i) = dconjg(a(i,j))
         enddo
      enddo

      return
      end



C ----------------------------------------------------------------------
C Library subroutine cpmxz
c
C copies a complex rectangular matrix to a different complex 
C rectangular matrix
C     c(j,i) = a(j,i)
c-----------------------------------------------------------------------

      subroutine cpmxz (a,c,dim1,dim2)

      implicit none

      integer     dim1,dim2,i,j
      complex*16  a(dim1,dim2),c(dim1,dim2)

      do i = 1,dim2
         do j = 1,dim1
            c(j,i) = a(j,i)
         enddo
      enddo

      return
      end

C ----------------------------------------------------------------------
c Library subroutine cpmxd
c
C copies a real rectangular matrix to a different real rectangular matrix
C     c(j,i) = a(j,i)
c-----------------------------------------------------------------------

      subroutine cpmxd (a,c,dim1,dim2)

      implicit none

      integer     dim1,dim2,i,j
      real*8      a(dim1,dim2),c(dim1,dim2)

      do i = 1,dim2
         do j = 1,dim1
            c(j,i) = a(j,i)
         enddo
      enddo

      return
      end

C ----------------------------------------------------------------------
C Library subroutine cpmxi
c
C copies an integer rectangular matrix to a different integer
C rectangular matrix
C     c(j,i) = a(j,i)
c-----------------------------------------------------------------------

      subroutine cpmxi (a,c,dim1,dim2)

      implicit none

      integer     dim1,dim2,i,j
      integer  a(dim1,dim2),c(dim1,dim2)

      do i = 1,dim2
         do j = 1,dim1
            c(j,i) = a(j,i)
         enddo
      enddo

      return
      end

C ----------------------------------------------------------------------
c Library subroutine cpvxd
c
C copies a real vector to a different real vector
C     w(i) = v(i)
c-----------------------------------------------------------------------

      subroutine cpvxd (v,w,dim)

      integer dim,i
      real*8  v(dim),w(dim)

      do i = 1,dim
         w(i) = v(i)
      enddo

      return
      end

C ----------------------------------------------------------------------
c Library subroutine cpvxz
c
C copies a complex vector to a different complex vector
C     w(i) = v(i)
c-----------------------------------------------------------------------

      subroutine cpvxz (v,w,dim)

      integer dim,i
      complex*16  v(dim),w(dim)

      do i = 1,dim
         w(i) = v(i)
      enddo

      return
      end

C ----------------------------------------------------------------------
c Libraray subroutine cpvcz
c
c copies a conjugate of a complex vector to a different complex vector
C     w(i) = dconjg(vi))
C ----------------------------------------------------------------------

      subroutine cpvcz (v,w,dim)

      implicit none

      integer dim,i
      complex*16 v(dim),w(dim)
      
      do i=1,dim
         w(i) = dconjg(v(i))
      enddo
      
      return
      end


C ----------------------------------------------------------------------
c Library subroutine cpvxdz
c
C copies a real vector to a  complex vector
C     w(i) = v(i)
c-----------------------------------------------------------------------

      subroutine cpvxdz (v,w,dim)

      integer dim,i
      real*8  v(dim)
      complex*16  w(dim)

      do i = 1,dim
         w(i) = v(i)
      enddo

      return
      end

C ----------------------------------------------------------------------
c Library subroutine cpvxzd
c
C copies the real part of a complex vector to a  real vector
C     w(i) = v(i)
c-----------------------------------------------------------------------

      subroutine cpvxzd (v,w,dim)

      integer dim,i
      real*8  w(dim)
      complex*16  v(dim)

      do i = 1,dim
         w(i) = dble(v(i))
      enddo

      return
      end

C ----------------------------------------------------------------------
c Library subroutine cpvxi
c
C copies a integer vector to a different integer vector
C     w(i) = v(i)
c-----------------------------------------------------------------------

      subroutine cpvxi (v,w,dim)

      integer dim,i
      integer  v(dim),w(dim)

      do i = 1,dim
         w(i) = v(i)
      enddo

      return
      end

C ----------------------------------------------------------------------
c Library subroutine cpvxi2
c
C replaces 2 integer vectors
C     w(i) = v(i)  &&  v(i) = w(i)
c-----------------------------------------------------------------------
 
      subroutine cpvxi2 (v,w,dim)
 
      integer dim,i,a
      integer  v(dim),w(dim)
 
      do i = 1,dim
         a = w(i)
         w(i) = v(i)
         v(i) = a
      enddo
 
      return
      end

C ----------------------------------------------------------------------
c Library subroutine cpvxl
c
C copies a logical vector to a different logical vector
C     w(i) = v(i)
c-----------------------------------------------------------------------

      subroutine cpvxl (v,w,dim)

      integer dim,i
      logical  v(dim),w(dim)

      do i = 1,dim
         w(i) = v(i)
      enddo

      return
      end

c-----------------------------------------------------------------------
c Library subroutine initvxz
c
C initialises a vector
C    vec(i)=(1.0d0,0.0d0)
c-----------------------------------------------------------------------

      subroutine initvxz(vec,dim)

      implicit none

      integer    dim,i
      complex*16 vec(dim)

      do i=1,dim
         vec(i)=(1.0d0,0.0d0)
      enddo

      return
      end

c-----------------------------------------------------------------------
c Library subroutine initvxd
c
C initialises a vector
C    vec(i)=1.0d0
c-----------------------------------------------------------------------

      subroutine initvxd(vec,dim)

      implicit none

      integer    dim,i
      real*8     vec(dim)

      do i=1,dim
         vec(i)=1.0d0
      enddo

      return
      end

c-----------------------------------------------------------------------
c Library subroutine initvxi
c
C initialises a vector
C    vec(i)=1
c-----------------------------------------------------------------------

      subroutine initvxi(vec,dim)

      implicit none

      integer    dim,i,vec(dim)

      do i=1,dim
         vec(i)=1
      enddo

      return
      end

c-----------------------------------------------------------------------
c Library subroutine initvxl
c
C initialises a vector
C    vec(i)=.true.
c-----------------------------------------------------------------------

      subroutine initvxl(vec,dim)

      implicit none

      integer    dim,i
      logical    vec(dim)

      do i=1,dim
         vec(i)=.true.
      enddo

      return
      end

c-----------------------------------------------------------------------
c Library subroutine initmxz
c
C initialises a matrix
C    a(j,i)=(1.0d0,0.0d0)
c-----------------------------------------------------------------------

      subroutine initmxz(a,dim1,dim2)

      implicit none

      integer    dim1,dim2,i,j
      complex*16 a(dim1,dim2)

      do i=1,dim2
         do j=1,dim1
            a(j,i)=(1.0d0,0.0d0)
         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
c Library subroutine initmxd
c
C initialises a matrix
C    a(j,i)=(1.0d0,0.0d0)
c-----------------------------------------------------------------------

      subroutine initmxd(a,dim1,dim2)

      implicit none

      integer    dim1,dim2,i,j
      real*8     a(dim1,dim2)

      do i=1,dim2
         do j=1,dim1
            a(j,i)=1.0d0
         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
c Library subroutine initmxl
c
C initialises a logical matrix
C    a(j,i)=.true.
c-----------------------------------------------------------------------

      subroutine initmxl(a,dim1,dim2)

      implicit none

      integer    dim1,dim2,i,j
      logical    a(dim1,dim2)

      do i=1,dim2
         do j=1,dim1
            a(j,i)=.true.
         enddo
      enddo

      return
      end


C-----------------------------------------------------------------------
C Library subroutine overmxz
C
C Overlap of a complex set of vectors with itself (i.e. multiplication
C of a matrix with its adjoint)
C     dconjg(a(k,j))*a(k,i) = c(j,i)
C NB resultant matrix is hermitian
C-----------------------------------------------------------------------
      subroutine overmxz (a,c,dim1,dim2)

      implicit none

      integer     dim1,dim2,i,j,k
      complex*16  a(dim1,dim2),c(dim2,dim2)

      do i = 1,dim2
         do j = i,dim2
            c(j,i) = dconjg(a(1,j))*a(1,i)
            do k = 2,dim1
               c(j,i) = c(j,i)+dconjg(a(k,j))*a(k,i)
            enddo
         enddo
      enddo
C
C now fill in other half of hermitian matrix
C
      do i=1,dim2
         c(i,i)=dble(c(i,i))
         do j=1,i-1
            c(j,i)=dconjg(c(i,j))
         enddo
      enddo

      return
      end




C-----------------------------------------------------------------------
C Library subroutine overmcz
C
C Overlap of a complex set of vectors with itself (i.e. multiplication
C of a matrix with its transpose)
C     a(k,j)*a(k,i) = c(j,i)
C NB standard definition takes c.c. of first matrix elements.
C NB resultant matrix is symmetric
C-----------------------------------------------------------------------

      subroutine overmcz (a,c,dim1,dim2)

      implicit none

      integer     dim1,dim2,i,j,k
      complex*16  a(dim1,dim2),c(dim2,dim2)

      do i = 1,dim2
         do j = i,dim2
            c(j,i) = a(1,j)*a(1,i)
            do k = 2,dim1
               c(j,i) = c(j,i)+a(k,j)*a(k,i)
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
C Library subroutine tranmxz
C
C transposes a complex matrix:
C   a(i,j) = c(j,i)
C-----------------------------------------------------------------------

      subroutine tranmxz(a,c,dim1,dim2)

      implicit none

      integer dim1,dim2,i,j
      complex*16 a(dim1,dim2),c(dim2,dim1)

      do i=1,dim1
         do j=1,dim2
            c(j,i)=a(i,j)
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine tranmxd
C
C transposes a real matrix:
C   a(i,j) = c(j,i)
C-----------------------------------------------------------------------

      subroutine tranmxd(a,c,dim1,dim2)

      implicit none

      integer dim1,dim2,i,j
      real*8 a(dim1,dim2),c(dim2,dim1)

      do i=1,dim1
         do j=1,dim2
            c(j,i)=a(i,j)
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine tranqxz
C
C transposes a complex quadratic matrix to itself:
C   a(i,j) = a(j,i)
C-----------------------------------------------------------------------

      subroutine tranqxz(a,dim)

      implicit none

      integer dim,i,j
      complex*16 a(dim,dim),b

      do i=2,dim
         do j=1,i-1
            b=a(i,j)
            a(i,j)=a(j,i)
            a(j,i)=b
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine tranqxd
C
C transposes a real quadratic matrix to itself:
C   a(i,j) = a(j,i)
C-----------------------------------------------------------------------

      subroutine tranqxd(a,dim)

      implicit none

      integer dim,i,j
      real*8 a(dim,dim),b

      do i=2,dim
         do j=1,i-1
            b=a(i,j)
            a(i,j)=a(j,i)
            a(j,i)=b
         enddo
      enddo

      return
      end

C-----------------------------------------------------------------------
C Library subroutine cpvxz_s
C
C copies a complex vector to a different complex vector, where the 
C vectors are not completeky stored, but managed by index arrays
C     w(i) = v(i)
C
C-----------------------------------------------------------------------


      subroutine cpvxz_s(v,w,index1,index2,dim1,dim2)

      implicit none

      integer dim1,dim2,index1(dim1),index2(dim2),
     +        b1,b2
      complex*16 v(dim1),w(dim2)

      b1=1
      do b2=1,dim2
  100    continue
         if (index2(b2) .eq. index1(b1)) then
            w(b2)=v(b1)
         else if (index2(b2) .gt. index1(b1) .and.
     +        b1 .lt. dim1) then
            b1=b1+1
            go to 100
         else
            w(b2)=0.0d0
         endif
      enddo

      return
      end


C-----------------------------------------------------------------------
C Library subroutine cuttxd2
C
C calculates a cut through a real*8 tensor of third order for a fixed given
c value of the second index jcut. the result is written to a real*8 matrix m.
c     m(i,k) = t(i,jcut,k)
C
C-----------------------------------------------------------------------


      subroutine cuttxd2(t,m,dim1,dim2,dim3,jcut)

      implicit none

      integer dim1,dim2,dim3,jcut,i,k
      real*8 t(dim1,dim2,dim3),m(dim1,dim3)

      if (dim1.eq.1) then
         do k=1,dim3
            m(1,k)=t(1,jcut,k)
         enddo
      else if (dim3.eq.1) then
         do i=1,dim1
            m(i,1)=t(i,jcut,1)
         enddo
      else
         do k=1,dim3
            do i=1,dim1
               m(i,k)=t(i,jcut,k)
            enddo
         enddo
      endif

      return
      end
