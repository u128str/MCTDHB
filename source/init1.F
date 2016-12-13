C***********************************************************************
C 
C                Module INIT1
C
C Contents: initho, initrho, initsin, initexp, initfft, initsphfbr 
C
C calls lapack-routine DSTEQR to diagonalise tri-diagonal matrices.
C
C Initialises the arrays trafo, dif2mat and dif1mat for a DVR basis set
C
C variables passed:
C  trafo   : DVR/FBR transformation matrices 
C  ort     : DVR grid points
C  dif1mat : DVR first  derivative matrix
C  dif2mat : DVR second derivative matrix
C  gdim    : number of DVR grid points for this mode
C  xi      : Value of first DVR grid point
C  xf      : Value of last  DVR grid point
C  hofm    : mass * vibrational frequency (homass*hofreq)
C  hoxeq   : equilibrium position of oscillator 
C  workr   : work array 
C
C
C***********************************************************************
C#######################################################################

C-----------------------------------------------------------------------
C                Subroutine INITHO
C             Harmonic Oscillator DVR
C-----------------------------------------------------------------------

      subroutine initho(trafo,ort,dif2mat,dif1mat,gdim
     +                   ,hofm,hoxeq,dx,hotyp,message)

      implicit none

      integer g,g1,g2,gdim,error,hotyp
      real*8 trafo(gdim,gdim),ort(gdim),dif2mat(gdim,gdim)
     +     ,dif1mat(gdim,gdim),hofm,hoxeq,x,dx
      character*(*) message

C-----------------------------------------------------------------------
C Calculate symmetric tridiagonal matrix representation X_mn of the
C coordinate operator (note that the diagonal elements are always zero).
C-----------------------------------------------------------------------
C --- Set diagonal elements to 0  

      call zerovxd(ort,gdim)

C --- Set super-diagonal elements for HO-FBR basis.(g=n+1)

      do g=1,gdim-1
         dif2mat(g,1)=dsqrt(g/(2.d0*hofm))
      enddo


C --- Diagonalize coordinate-matrix to obtain "ort" and "trafo".
      call DSTEQR('I',gdim,ort,dif2mat,trafo,gdim,dif1mat,error)
      if (error.ne.0) then
         write(6,*) 'Error in diagonalising HO-coordinate-matrix'
c str         call errormsg
      endif

C --- If xi-xf flag is set
 
      if(hotyp .eq. 1) then
         hofm = dx/(ort(gdim)-ort(1))
         do  g = 1, gdim
            ort(g) = ort(g)*hofm
         end do
         hofm = hofm**(-2)
      end if

C-----------------------------------------------------------------------
C Change the sign of the eigenvector if necessary, and thus
C  ensure that the DVR-weights are all positive. (See genphi1.f)
C-----------------------------------------------------------------------
      do g=1,gdim
         if(trafo(1,g) .lt. 0.d0 ) then
            do g1=1,gdim
               trafo(g1,g) = -trafo(g1,g) 
            enddo
         end if
      enddo

C-----------------------------------------------------------------------
C Calculate second derivative matrix in corresponding DVR basis
C
C Harmonic oscillator: H=T+V in FBR basis and V in DVR basis are
C analytically known --> define 2*mass*H in FBR basis, transform to
C DVR basis and subtract 2*mass*V (diagonally) given in DVR basis to
C get  d^2/dQ^2 = -2*mass*T = -2*mass*(H-V) in DVR basis.
C-----------------------------------------------------------------------

      call zeromxd(dif2mat,gdim,gdim)

C --- Total energy in FBR basis (given analytically and diagonal)

      do g=1,gdim
         x=(g-0.5d0)*hofm*2.d0

C --- Transform to DVR basis 

         do g1=1,gdim
            do g2=1,gdim
               dif2mat(g2,g1)=dif2mat(g2,g1)-trafo(g,g2)*x*trafo(g,g1)
            enddo
         enddo
      enddo

C --- Subtract potential contribution for harmonic oscillator

      do g=1,gdim
         dif2mat(g,g) = dif2mat(g,g) + (hofm*ort(g))*(hofm*ort(g))
      enddo
 

C-----------------------------------------------------------------------
C Calculate DVR matrix representation for first derivative.
C
C  The FBR representation of d/dx is antisymmetric tridiagonal with 
C  zero diagonal elements and analytically given as (assuming hbar=1): 
C
C       <ho_n|d/dx|ho_m> = -sqrt(mass*freq/2)*
C                         [ sqrt(m+1) delta_nm+1 - sqrt(m) delta_nm-1]
C
C-----------------------------------------------------------------------

      call zeromxd(dif1mat,gdim,gdim)
      
      do g=1,gdim-1

C --- Get super-diagonal matrix element in FBR basis

         x = -dsqrt(hofm*g/2.d0)

C --- Transform from FBR to DVR basis

         do g1=1,gdim
            do g2=1,gdim
               dif1mat(g2,g1)=dif1mat(g2,g1)+trafo(g+1,g2)
     +              *trafo(g,g1)*x-trafo(g,g2)*trafo(g+1,g1)*x
            enddo
         enddo
      enddo


C --- Shift grid points

      do g=1,gdim
         ort(g)=ort(g)+hoxeq
      enddo
      return
      end


C#######################################################################

C-----------------------------------------------------------------------
C                Subroutine INITRHO
C            Radial Harmonic Oscillator DVR
C-----------------------------------------------------------------------

      subroutine initrho(trafo,ort,dif2mat,dif1mat,gdim
     +                   ,hofm,hoxeq,xi,dx,hotyp,message)

      implicit none

      integer g,g1,g2,gdim,error,hotyp
      real*8 trafo(gdim,gdim),ort(gdim),dif2mat(gdim,gdim)
     +     ,dif1mat(gdim,gdim),hofm,hoxeq,x,xi,dx
      character*(*) message

C-----------------------------------------------------------------------
C phi_j = sqrt(2)*phi^{HO}_2j-1 ; odd harmonic oscillator functions.
C Calculate symmetric tridiagonal matrix representation X_mn of the
C transformed coordinate operator z = (x-x_0)^2 .
C-----------------------------------------------------------------------

C --- Set superdiagonal elements for RHO-FBR basis.

         do g=1,gdim-1
            dif2mat(g,1)=dsqrt(dble(2*g*(2*g+1)))/(2.d0*hofm)
         enddo

C --- Set diagonal elements for RHO-FBR basis.

         do g=1,gdim
            ort(g) = (4*g-1)/(2.d0*hofm)
         enddo

C --- Diagonalize coordinate-matrix to obtain "ort" and "trafo".

      call DSTEQR('I',gdim,ort,dif2mat,trafo,gdim,dif1mat,error)

      if (error.ne.0) then
         message = 'Error in diagonalising HO-coordinate-matrix'
c          call errormsg
      endif

C --- Transform ort
      
      do g=1,gdim
         ort(g) = dsqrt(ort(g))
      enddo

C --- If xi-xf flag is set
 
      if(hotyp .eq. 1) then
         hofm  = dx/(ort(gdim)-ort(1))
         do  g = 1, gdim
            ort(g) = ort(g)*hofm 
         end do
         hoxeq = xi - ort(1)
         hofm  = hofm**(-2)
      end if

C-----------------------------------------------------------------------
C Change the sign of the eigenvector if necessary, and thus
C ensure that the DVR-weights are all positive. (See genphi1.f)
C-----------------------------------------------------------------------
      do g=1,gdim
         if(trafo(1,g) .lt. 0.d0 ) then
            do g1=1,gdim
               trafo(g1,g) = -trafo(g1,g) 
            enddo
         end if
      enddo

C-----------------------------------------------------------------------
C Calculate second derivative matrix in corresponding DVR basis
C
C RHarmonic oscillator: H=T+V in FBR basis and V in DVR basis are
C analytically known --> define 2*mass*H in FBR basis, transform to
C DVR basis and subtract 2*mass*V (diagonally) given in DVR basis to
C get  d^2/dQ^2 = -2*mass*T = -2*mass*(H-V) in DVR basis.
C-----------------------------------------------------------------------

      call zeromxd(dif2mat,gdim,gdim)

C --- Total energy in FBR basis (given analytically and diagonal)

      do g=1,gdim
         x=(2*g-0.5d0)*hofm*2.d0

C --- Transform to DVR basis 

         do g1=1,gdim
            do g2=1,gdim
               dif2mat(g2,g1)=dif2mat(g2,g1)-trafo(g,g2)*x*trafo(g,g1)
            enddo
         enddo
      enddo

C --- Subtract potential contribution for harmonic oscillator

      do g=1,gdim
         dif2mat(g,g) = dif2mat(g,g) + (hofm*ort(g))**2
      enddo
 

C-----------------------------------------------------------------------
C Calculate DVR matrix representation for first derivative.
C
C  The FBR representation of 0.5*(x*d/dx+d/dx*x) is antisymmetric 
C  tridiagonal with zero diagonal elements and analytically given as: 
C
C       <ho_{2n-1}|d/dx|ho_{2m-1}> = 
C           sqrt(2m*(2m+1)) delta_nm+1 - sqrt(2m-1)*(2m-2)) delta_nm-1
C
C-----------------------------------------------------------------------

      call zeromxd(dif1mat,gdim,gdim)
      
      do g=1,gdim-1

C --- Get super-diagonal matrix element in FBR basis

         x = -dsqrt(dble(2*g*(2*g+1)))

C --- Transform from FBR to DVR basis

         do g1=1,gdim
            do g2=1,gdim
               dif1mat(g2,g1)=dif1mat(g2,g1)+trafo(g+1,g2)
     +              *trafo(g,g1)*x-trafo(g,g2)*trafo(g+1,g1)*x
            enddo
         enddo
      enddo


C --- Shift grid points

      do g=1,gdim
         ort(g)=ort(g)+hoxeq
      enddo

      return
      end



C-----------------------------------------------------------------------
C                Subroutine INITSIN
C                   Sine DVR
C-----------------------------------------------------------------------

      subroutine initsin(trafo,ort,dif2mat,dif1mat,gdim,xi,xf,workr)

      implicit none

      integer g1,g2,gdim
      real*8 trafo(gdim,gdim),ort(gdim),dif2mat(gdim,gdim)
     +     ,dif1mat(gdim,gdim),workr(gdim,gdim),
     +     x,xi,xf,deltax,fac1,Pi

      Pi     = 4.d0*atan(1.d0)
      deltax = (xf-xi)/dble(gdim-1)
      x      = gdim+1.d0
      fac1   = Pi/x

C --- Grid points are analytically given

      do g1=1,gdim
         ort(g1)=xi+(g1-1)*deltax
      enddo         

C --- DVR/FBR transformation matrix is analytically given

      do g1=1,gdim
         do g2=1,gdim
            trafo(g2,g1)= dsqrt(2d0/x)*sin(g2*g1*fac1)
         enddo
      enddo

C --- Second derivative in Sine-DVR basis (given analytically)

      do g1=1,gdim
         do g2=1,g1-1
            dif2mat(g2,g1)=-(Pi/deltax)**2
     +           *2d0*(-1d0)**(g2-g1)/x**2
     +           *sin(g2*fac1)*sin(g1*fac1)
     +           /(cos(g2*fac1)-cos(g1*fac1))**2
            dif2mat(g1,g2)=dif2mat(g2,g1)
         enddo
         dif2mat(g1,g1)=-(Pi/deltax)**2
     +        *(1d0/3d0+1d0/(6d0*x**2)
     +        -1d0/(2d0*x**2*sin(g1*fac1)**2))
      enddo        

C-----------------------------------------------------------------------
C Derivative in FBR basis is analytically given for SIN-DVR
C but is not tridiagonal:
C   i)  all diagonal elements and all elements for which the row
C       and the column index are simultaneously odd or even
C       are zero
C   ii) for all other elements it holds:
C
C           d_ab = 4ab/(a**2-b**2) 1/(gdim+1) 1/deltax
C
C   iii) from ii) it follows that the matrix is still antisymmetric
C-----------------------------------------------------------------------
         

      fac1 = 4.d0/((gdim+1)*deltax)
      do g1=1,gdim 
         do g2=1,g1-1
            if (mod(g2+g1,2).eq.0) then
               dif1mat(g2,g1) = 0.d0
            else
               dif1mat(g2,g1) = fac1*g2*g1/(g2**2-g1**2)
            endif
            dif1mat(g1,g2) = -dif1mat(g2,g1)
         enddo
         dif1mat(g1,g1) = 0.d0 
      enddo
      
C --- use special spin form, or transform from FBR to DVR basis
      if(gdim.eq.2 .and. abs(xi+0.5d0).lt.1.0d-9 .and. 
     +                   abs(xf-0.5d0).lt.1.0d-9 ) then
         dif1mat(1,1) = 0.d0
         dif1mat(1,2) = -0.5d0
         dif1mat(2,1) = 0.5d0
         dif1mat(2,2) = 0.d0
         dif2mat(1,1) = 0.d0
         dif2mat(1,2) = 0.5d0
         dif2mat(2,1) = 0.5d0
         dif2mat(2,2) = 0.d0
         PRINT * , 'SPIN Form'
      else
         call qqxxdd(dif1mat,trafo,workr,gdim)
         call qqtxdd(trafo,workr,dif1mat,gdim)
      endif

      return
      end



C#######################################################################

C-----------------------------------------------------------------------
C                Subroutine INITEXP
C                  Exponential DVR
C-----------------------------------------------------------------------

      subroutine initexp(trafo,ort,dif2mat,dif1mat,exphin,exprueck,gdim,
     +                   xi,xf,kmin,kdim)

      implicit none

      integer g,g1,g2,gdim,m,kmin,kdim
      real*8 trafo(gdim,gdim),ort(gdim),dif2mat(gdim,gdim)
     +     ,dif1mat(gdim,gdim),xi,xf,deltax,fac1,fac2,Pi
      complex*16 exphin(gdim,kdim),exprueck(kdim,gdim)


      if( mod(gdim,2) .eq. 0 )  then
         write(6,*) ' ERROR in init1.f : '
         write(6,*) 'Number of grid points must be odd for exp-DVR'
         stop
      end if

      Pi     = 4.d0*atan(1.d0)
      deltax = (xf-xi)/dble(gdim-1)
      fac1   = Pi/(gdim*deltax)
      fac2   = Pi/dble(gdim)


C --- Grid points are analytically given

      do g=1,gdim
         ort(g)=xi+(g-1)*deltax
      enddo         

CC --- DVR/FBR transformation matrix is analytically given, however,
CC --- it is complex. Thus the correct trafo is outcommented and 
CC --- trafo is zeroed.

      call zeromxd(trafo,gdim,gdim)

c     do g1=1,gdim
c        do g2=1,gdim
C           trafo(g2,g1)= exp(fac2*(0.d0,2.d0)*g1*g2)/sqrt(gdim)
c        enddo
c     enddo

      do g1=1,kdim
         m=g1-1+kmin
         do g2=1,gdim
            exphin(g2,g1) = exp(-dcmplx(0.d0,m*ort(g2)))/
     +                      dsqrt(dble(gdim))
            exprueck(g1,g2) = exp(dcmplx(0.d0,m*ort(g2)))/
     +                        dsqrt(dble(gdim))
         enddo
      enddo

C --- Second derivative in Exponential-DVR basis (given analytically)

      do g1=1,gdim
         do g2=1,g1-1
            g = g1-g2
            dif2mat(g1,g2) = -2.d0*fac1**2 *( (-1)**g ) 
     $                   * cos(fac2*g)/(sin(fac2*g))**2
            dif2mat(g2,g1) =  dif2mat(g1,g2)   
         enddo
         dif2mat(g1,g1) = -fac1**2 *(gdim**2-1)/3.d0
      enddo


C --- First derivative in Exponential-DVR basis (given analytically)

      do g1=1,gdim
         do g2=1,g1-1
            g = g1-g2
            dif1mat(g1,g2) =  fac1*( (-1)**g ) /sin(fac2*g)
            dif1mat(g2,g1) =  -dif1mat(g1,g2)   
         enddo
         dif1mat(g1,g1) = 0.d0
      enddo

      return
      end




C#######################################################################

C***********************************************************************
C
C                   Subroutine INITFFT
C
C calls fftinit
C
C Initialises the arrays ort,hin,rueck,fftfak,fftp for FFT 
C
C variables passed for one mode:
C  ort   : FFT grid points
C  hin   : array for position --> momentum FFT transformation
C  rueck : array for momentum --> position FFT transformation
C  fftfak: FFT decomposition factors 
C  fftp  : FFT momentum array
C  gdim  : number of FFT grid points
C  xi    : First point of FFT grid
C  xf    : last point  of FFT grid
C
C   11/96 AJ
C
C***********************************************************************

      subroutine initfft(ort,fftp,gdim,xi,xf)

      implicit none

      integer gdim,g
      real*8  xi,xf,deltax,length,deltap,pmax,Pi
      real*8  ort(0:gdim-1), fftp(0:gdim-1)

      Pi=4.d0*atan(1.d0)

C-----------------------------------------------------------------------
C Define grid
C-----------------------------------------------------------------------
      deltax = (xf-xi)/(gdim-1)
      do g=0,gdim-1
         ort(g)=xi+g*deltax
      enddo

C-----------------------------------------------------------------------
C Calculate momentum grid
C-----------------------------------------------------------------------
      length = gdim*deltax
      deltap = 2d0*Pi/length
      pmax   = Pi/deltax 
      if (mod(gdim,2).eq.0) then
         do g=0,gdim/2-1
            fftp(g)=g*deltap
            fftp(g+gdim/2)=-pmax+g*deltap
         enddo
      else
         do g=0,gdim/2
            fftp(g)=g*deltap
         enddo
         do g=0,gdim/2-1
            fftp(g+gdim/2+1)=-pmax+(g+0.5d0)*deltap
         enddo
      endif      

      return
      end        

C#######################################################################

C-----------------------------------------------------------------------
C
C                   Subroutine initsphfbr
C                 spherical  harmonics  FBR
C
C
C Initialises the arrays kinsph,jsph,msph for a FBR representation
C in a spherical harmonics basis set.
C
C variables passed for one mode:
C  kinsph : FBR kinetic energy matrix (diagonal matrix)
C  jsph   : azimutal quantum numbers for the basis functions
C  msph   : magnetic quantum numbers for the basis functions
C  gdim   : number of FBR basis functions
C  jmax   : maximum azimutal quantum number for spherical harmonics
C  sphtyp : type of spherical harmonics basis :  1: j,m=all (abs(m)<= jmax)
C           2: j=selected and m=all (j=0,2.. or 1,3..)
C           3: j=selected and m=truncated
C                   (j=0,2.. or 1,3..) and (abs(m) <= mmax)
C           4: j=all  and m=truncated  (abs(m) <= mmax)
C  mmax   : maximum magnetic quantum number for truncated basis
C  mincr  : increment of m's, starting from mmax, if m=truncated
C  f      : fth mode
C
C   08/97 MCH
C
C-----------------------------------------------------------------------

      subroutine initsphfbr(gdim,kinsph,jzsph,jzzsph,jsph,msph,
     +                      jmax,sphtyp,mmax,mincr,mmin,ilog)

      implicit none

      integer gdim,jsph(gdim),msph(gdim),jmax,sphtyp,ilog
      integer m,j0,j,g,xx,mmax,mincr,mmin
      real*8  kinsph(gdim),jzsph(gdim),jzzsph(gdim)


C  quantum numbers arrays jsph,msph are defined
C  according to the type of basis.
      g=0
      if (sphtyp.le.2)  then
         if (sphtyp.eq.1) then
            j0=0
         else
            j0=mod(jmax,2)
         endif

         do j=j0,jmax,sphtyp
            do m=-j,j,1
               g=g+1
               jsph(g)=j
               msph(g)=m
            enddo
         enddo

      else

         if (sphtyp.eq.3) then
           xx=2
           j0=mod(jmax,2)
         else
           xx=1
           j0=0
         endif

         do j=j0,jmax,xx
            do m=mmin,mmax,mincr
               if (iabs(m).le.min(mmax,j)) then
                  g=g+1
                  jsph(g)=j
                  msph(g)=m
               endif
            enddo
         enddo

      endif

c to debug
      if (gdim.ne.g) then
         write(ilog,'(a,2i5)') '+++ problem in initsphfbr ',gdim,g
         stop
      endif

C Construction of the j^2, j_z and j_z^2 matrices in the FBR-representation
C diagonal matrices stored as an one-dimensional arrays.

      do g=1,gdim
         kinsph(g)=jsph(g)*(jsph(g)+1)
         jzsph (g)=msph(g)
         jzzsph(g)=msph(g)**2
      enddo

      return
      end

C#######################################################################

C-----------------------------------------------------------------------
C                Subroutine INITLAGU
C                  Laguerre DVR
C phi_n(x) = (x/n)**(1/2) * exp(-x/2) * L^{1}_{n-1}(x); (alpha=1)
C phi_n(x) = A_n * x**(alpha/2) * exp(-x/(2*breite)) * 
C                                 L^{alpha}_{n-1}(x/breite)
C-----------------------------------------------------------------------

      subroutine initlagu(trafo,ort,dif2mat,dif1mat,workr,gdim,
     +                    x0,breite,xi,xf,typ,icut,message,basis)

      implicit none

      integer g,g1,gdim,typ,icut,error,basis,alpha
      real*8  trafo(gdim,gdim),ort(gdim),workr(gdim,gdim),
     +        dif1mat(gdim,gdim),dif2mat(gdim,gdim),
     +        x0,breite,xi,xf
      character*(*) message
      character*5  clgu


C-----------------------------------------------------------------------
C  Calculate second derivative matrix in corresponding DVR basis
C  (d^2/dx^2 - c/x^2) in Lagu-FBR basis is analytically evaluated. 
C  c = -1/4, 0, 3/4, 2  for alpha = 1,2,3,4 .
C  The centrifugal term c/x^2 is finally removed in DVR representation. 
C-----------------------------------------------------------------------
      if(basis.eq.15) then
         alpha = 1
         clgu  = 'Lagu1'
         do g=1,gdim
            trafo(g,g) = 0.25d0*(1-2*g)/breite**2
            do g1=1,g-1
               trafo(g,g1) = -0.5d0*dsqrt(dble(g1)/dble(g))*g1/breite**2
               trafo(g1,g) = trafo(g,g1)
            enddo
         enddo
      elseif(basis.eq.16) then
         alpha = 2
         clgu  = 'Lagu2'
         do g=1,gdim
            trafo(g,g) = dble(1-4*g)/(12.d0*breite**2)
            do g1=1,g-1
               trafo(g,g1) = -dsqrt(dble(g1*(g1+1))/dble(g*(g+1)))
     +             * dble(2*g1+1)/(6.d0*breite**2)
               trafo(g1,g) = trafo(g,g1)
            enddo
         enddo
      elseif(basis.eq.17) then
         alpha = 3
         clgu  = 'Lagu3'
         do g=1,gdim
            trafo(g,g) = -g/(4.d0*breite**2)
            do g1=1,g-1
               trafo(g,g1) = dsqrt(dble(g1*(g1+1)*(g1+2))/
     +                            dble(g*(g+1)*(g+2)))
     +             * dble(-g1-1)/(4.d0*breite**2)
               trafo(g1,g) = trafo(g,g1)
            enddo
         enddo
      elseif(basis.eq.18) then
         alpha = 4
         clgu  = 'Lagu4'
         do g=1,gdim
            trafo(g,g) = -dble(4*g+1)/(20.d0*breite**2)
            do g1=1,g-1
               trafo(g,g1) = dsqrt((dble(g1*(g1+1))/dble(g*(g+1)))*
     +                    (dble((g1+2)*(g1+3))/dble((g+2)*(g+3))))
     +             * dble(-2*g1-3)/(10.d0*breite**2)
               trafo(g1,g) = trafo(g,g1)
            enddo
         enddo
      endif

C-----------------------------------------------------------------------
C In order to remove extremly large eigenvalues of dif2mat, this
C matrix is diagonalized, the largest (negative) eigenvalues are 
C modified and the diagonal matrix is transformed back.
C-----------------------------------------------------------------------
C get all eigenvalues, stored in "ort" and the corresponding normalized
C eigenvectors of the second derivatives matrix "trafo". 
C On output, the eigenvectors are store in "trafo".
C-----------------------------------------------------------------------
      g = gdim*gdim   
      call dsyev('V','U',gdim,trafo,gdim,ort,workr,g,error)

      if (error.ne.0) then
         message = 'Error in diagonalising dif2mat (Lagu)'
c         call errormsg
      endif

C---- Modify eigenvalues.
      write(2,'(2a,i3,3(a,e12.4))') clgu,': icut =',icut,
     +        ', EV :',-ort(1), ', cut :', -ort(icut+1),
     +        ', EV/cut :',ort(1)/ort(icut+1)
      do g = 1, icut
         ort(g) = ort(icut+1)
      enddo

C-----------------------------------------------------------------------
C dif2mat = trafo * ort * trafo-transposed.
C-----------------------------------------------------------------------
c      call qdxxdd (trafo,ort,workr,gdim)
c      call qqxtdd (trafo,workr,dif2mat,gdim)
      
C-----------------------------------------------------------------------
C Calculate symmetric tridiagonal matrix representation X_mn of the
C x operator (note that the diagonal elements are not zero).
C-----------------------------------------------------------------------

C --- Set diagonal elements of X_mn 

      do g=1,gdim
         ort(g) = 2*g + alpha -1
      enddo

C --- Set subdiagonal elements of X_mn.

      do g=1,gdim-1
         workr(g,1) = -dsqrt(dble(g*(g+alpha)))
      enddo

C --- Diagonalize x-matrix to obtain "ort" and "trafo".

      call DSTEQR('I',gdim,ort,workr,trafo,gdim,dif1mat,error)

      if (error.ne.0) then
         message = 'Error in diagonalising X-matrix (Lagu)'
c         call errormsg
      endif

C-----------------------------------------------------------------------
C Change the sign of the eigenvector if necessary, and thus
C  ensure that the DVR-weights are all positive. (See genphi1.f)
C-----------------------------------------------------------------------
      do g=1,gdim
         if(trafo(1,g) .lt. 0.d0 ) then
            do g1=1,gdim
               trafo(g1,g) = -trafo(g1,g)
            enddo
         end if
      enddo

C-------------------------------------------------------------------------
C   Transform to DVR basis
C-------------------------------------------------------------------------
c      call qqxxdd(dif2mat,trafo,workr,gdim)
c      call qqtxdd(trafo,workr,dif2mat,gdim)

C-------------------------------------------------------------------------
C Remove the artifical centrifugal potential c/x^2  (c = -1/4, 0, 3/4, 2).
C-------------------------------------------------------------------------
      if(basis.eq.15) then
         do g=1,gdim
            dif2mat(g,g) = dif2mat(g,g) - 0.25d0/(ort(g)*breite)**2
         enddo
      elseif(basis.eq.17) then
         do g=1,gdim
            dif2mat(g,g) = dif2mat(g,g) + 0.75d0/(ort(g)*breite)**2
         enddo
      elseif(basis.eq.18) then
         do g=1,gdim
            dif2mat(g,g) = dif2mat(g,g) + 2.0d0/(ort(g)*breite)**2
         enddo
      endif

C-------------------------------------------------------------------------
C  Calculate first derivative matrix in corresponding DVR basis
C  d/dx in Lagu-FBR basis is analytically known  
C   --> transform to DVR basis
C-------------------------------------------------------------------------
      if(basis.eq.15) then
         do g=1,gdim
            dif1mat(g,g) = 0.d0
            do g1=1,g-1
               dif1mat(g,g1) = 0.5d0*dsqrt(dble(g1)/dble(g))/breite
               dif1mat(g1,g) = -dif1mat(g,g1)
            enddo
         enddo
      elseif(basis.eq.16) then
         do g=1,gdim
            dif1mat(g,g) = 0.d0
            do g1=1,g-1
              dif1mat(g,g1) = 0.5d0*dsqrt(dble(g1*(g1+1))/dble(g*(g+1)))
     +                         /breite
              dif1mat(g1,g) = -dif1mat(g,g1)
            enddo
         enddo
      elseif(basis.eq.17) then
         do g=1,gdim
            dif1mat(g,g) = 0.d0
            do g1=1,g-1
               dif1mat(g,g1) = 0.5d0*dsqrt(dble(g1*(g1+1)*(g1+2))/
     +                                dble(g*(g+1)*(g+2)))/breite
               dif1mat(g1,g) = -dif1mat(g,g1)
            enddo
         enddo
      elseif(basis.eq.18) then
         do g=1,gdim
            dif1mat(g,g) = 0.d0
            do g1=1,g-1
              dif1mat(g,g1) =0.5d0*dsqrt((dble(g1*(g1+1))/dble(g*(g+1)))
     +                * (dble((g1+2)*(g1+3))/dble((g+2)*(g+3))))/breite
               dif1mat(g1,g) = -dif1mat(g,g1)
            enddo
         enddo
      endif

c      call qqxxdd(dif1mat,trafo,workr,gdim)
c      call qqtxdd(trafo,workr,dif1mat,gdim)

C-----------------------------------------------------------------------
C If xi-xf or x0-xf is set
C-----------------------------------------------------------------------
      if(typ.gt.0) then
         if(typ.eq.1) then
            breite = (xf-xi)/(ort(gdim)-ort(1))
            x0 = xi - ort(1)*breite
         elseif(typ.eq.2) then
            breite = (xf-x0)/ort(gdim)
         endif
         do g=1,gdim
            do g1=1,gdim
               dif1mat(g,g1) = dif1mat(g,g1)/breite
               dif2mat(g,g1) = dif2mat(g,g1)/breite**2
            enddo
         enddo
      endif

C-----------------------------------------------------------------------
C  Shift grid points 
C-----------------------------------------------------------------------
      do g=1,gdim
         ort(g) = ort(g)*breite + x0
      enddo

      xi = ort(1)
      xf = ort(gdim)

      return
      end
C-----------------------------------------------------------------------

