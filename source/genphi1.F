C***********************************************************************
C 
C            MCTDH MODULE GENPHI1
C                
C***********************************************************************
C 
C            MCTDH subroutine GENDISC
C                
C 
C sets up single particle functions for a discrete basis set. 
C The spfs are a cmoplete set of orthonormal vectors. 
C
C variables passed for one mode:
C  psi: single particle functions (states)
C  dim: number of spfs (states)
C  gdim: as dim
C
C-----------------------------------------------------------------------
C basis(f) : basis type for dof f. (dof=degree of freedom)
C(phifbr), el, HO, Leg, sin, FFT, exp, sphfbr, kleg, k, pleg, rHO, Leg/r
C   -1      0   1    2    3    4    5      6    7    8    9    11   12  
C***********************************************************************

      subroutine gendisc (psi,dim,gdim)

      implicit none

      integer dim,gdim,e,g
      complex*16 psi(gdim,dim)
      
      do e=1,dim
         do g=1,gdim
            if (e.ne.g) then
               psi(g,e)=(0.d0,0.d0)
            else
               psi(g,e)=(1.d0,0.d0)
            endif
        enddo
      enddo

      return
      end


C***********************************************************************
C
C             MCTDH  subroutine  GENGAUSS
C
C  11/96  AJ
C
C***********************************************************************

      subroutine gengauss(psi,gdim,dim,weight,trafo,ort,basis,rpbaspar,
     +           ipbaspar,ipbaspar1,rsbaspar,isbaspar,ggdim)

      implicit none

      integer ione 
      parameter (ione=1)

      integer e,g,gdim,dim,basis,next,ipbaspar(*),blz,ibleg,dim1,
     +        ipbaspar1(*),typ,power,isbaspar(*),ierr,ggdim
      real*8  trafo(ggdim,gdim),ort(gdim),weight(gdim),
     +        gxeq,gmom,norm,homass,hofreq,hoxeq,ortg,shift,
     +        rpbaspar(*),rsbaspar(*)
      complex*16 psi(gdim,dim),gwidth

C-----------------------------------------------------------------------
C get primitive basis (DVR) parameters
C-----------------------------------------------------------------------
      if (basis .eq. 1 .or. basis .eq. 11) then
         hoxeq  = rpbaspar(1)
         hofreq = rpbaspar(2)
         homass = rpbaspar(3)
      else if (basis .eq. 2 .or. basis .eq. 12) then
         blz=ipbaspar(1)
         ibleg=ipbaspar(2)
      else if (basis .eq. 7) then
         blz=ipbaspar1(3)
         ibleg=ipbaspar(1)
      else if (basis .eq. 9) then
         blz=ipbaspar(5)
         ibleg=ipbaspar(1)
      endif

C-----------------------------------------------------------------------
C get DVR weigts
C-----------------------------------------------------------------------
      call dvrweights(trafo,ort,gdim,ggdim,weight,basis,
     +                rpbaspar,ipbaspar,ipbaspar1,.false.,.false.)

C-----------------------------------------------------------------------
C get type (i.e. all or only odd / even functions)
C-----------------------------------------------------------------------
      typ=isbaspar(1)

C-----------------------------------------------------------------------
C get spf basis parameters
C-----------------------------------------------------------------------
      gxeq   = rsbaspar(1)
      gmom   = rsbaspar(2)
      gwidth = dcmplx(rsbaspar(3),rsbaspar(4))
      shift  = 0.d0 
      if(isbaspar(5).eq.1) shift = ort(2) + ort(gdim) - 2.d0*ort(1)

C-----------------------------------------------------------------------
C Multiply initial Gaussian by DVR weights
C-----------------------------------------------------------------------
      do g=1,gdim
         ortg = ort(g)-gxeq
         if( isbaspar(5) .eq. 1 ) then
            if(abs(ortg) .gt. abs(ortg-shift)) ortg = ortg-shift
            if(abs(ortg) .gt. abs(ortg+shift)) ortg = ortg+shift
            if(abs(ortg) .gt. abs(ortg-shift)) ortg = ortg-shift
         endif

         psi(g,1) = weight(g)*exp(-0.25d0*(ortg*gwidth)**2)

         if (typ .eq. 1) then
C.....if odd functions requested, set up x*psi
            psi(g,1)=ortg*psi(g,1)
         endif
      enddo

C-----------------------------------------------------------------------
C Normalise initial gaussian
C-----------------------------------------------------------------------
      call normvxz(psi,norm,gdim)
      call xvixdzo(norm,psi,gdim)

C-----------------------------------------------------------------------
C Set up orthonormal higher initial single-particle functions:
C  next is a parameter to say that 2nd spf is first function to be
C  generated
C-----------------------------------------------------------------------
      next=2
      dim1=min(dim,gdim)
      power=1
      if (typ .gt. 0) power=2
      call genhigher(psi,gdim,dim1,next,ort,basis,ipbaspar,power,
     +                gxeq,shift)
      call schmidtortho(psi,gdim,dim1,ierr)
     
C-----------------------------------------------------------------------
C     If the gauss is too localized, it is impossible to create a large
C     number of linearly independent spf's. schmidtortho then returns
C     with ierr=number of function that could to be orthonormalized.
C-----------------------------------------------------------------------
      if( ierr .gt. 0 ) then
         next = ierr+1
         do g=1,gdim
            psi(g,ierr)=(ort(g)-gxeq)**typ
         enddo
         call normvxz(psi(1,ierr),norm,gdim)
         call xvixdzo(norm,psi(1,ierr),gdim)
         call genhigher(psi,gdim,dim1,next,ort,basis,ipbaspar,
     +                  power,gxeq,shift)
         call schmidtortho(psi,gdim,dim1,ierr)
         if( ierr .gt. 0 ) then
            write(6,'(/a,/a,i4,a,/a,i4,a,i4,a,i4)')
     +     '############# ERROR in gengauss (genphi1.F) #############', 
     +              'Cannot Schmidt-orthogonalize the',
     +               ierr,'-th spf (build).',
     +              ' basis =',basis,',  gdim =', gdim,', dim =',dim 
 
            write(2,'(/a,/a,i4,a,/a,i4,a,i4,a,i4)') 
     +     '############# ERROR in gengauss (genphi1.F) #############', 
     +              'Cannot Schmidt-orthogonalize the',
     +               ierr,'-th spf (build).',
     +              ' basis =',basis,',  gdim =', gdim,', dim =',dim 
            stop 1
         end if
      end if



C-----------------------------------------------------------------------
C Multiply by factor for initial momentum 
C-----------------------------------------------------------------------
      do e=1,dim
         do g=1,gdim
            psi(g,e)=psi(g,e)
     +           *exp((0,1)*gmom*(ort(g)-gxeq))
         enddo
      enddo
c str      write(6,*)"HOFUNK FROM GETGAUSS",psi
      return
      end

C***********************************************************************
C
C             MCTDH  subroutine  GENALEG
C
C
C    Generates a set of the first 'dim' normalized associated
C    Legendre polynomials P_l^n for fixed slz
C
C    10/96 AJ
C
C***********************************************************************

      subroutine genaleg(psi,gdim,dim,weight,trafo,ort,basis,ggdim,
     +                   rpbaspar,ipbaspar,ipbaspar1,isbaspar,message)

      implicit none

      integer e,g,gdim,dim,basis,slz,sl0,isleg,ibleg,j,ipbaspar(*),
     +        ipbaspar1(*),isbaspar(*),isleg1,ggdim
      real*8  ort(gdim),trafo(ggdim,gdim),weight(gdim),rpbaspar(*)
      real*8  dlegendre
      complex*16 psi(gdim,dim)
      character*(*) message
      logical lerr

c      external dlegendre

C-----------------------------------------------------------------------
C get primitive basis (DVR) parameters
C-----------------------------------------------------------------------
      if (basis .eq. 2  .or.  basis .eq. 12) then
         ibleg = ipbaspar(2)
      elseif (basis .eq. 7) then
         ibleg = ipbaspar(1)
      elseif (basis .eq. 9) then
         ibleg = ipbaspar(1)
      endif

C-----------------------------------------------------------------------
C get DVR weigts
C-----------------------------------------------------------------------
      call dvrweights(trafo,ort,gdim,ggdim,weight,basis,
     +                rpbaspar,ipbaspar,ipbaspar1,.false.,.false.)

C-----------------------------------------------------------------------
C get spf basis parameters
C-----------------------------------------------------------------------
      slz   = isbaspar(1)
      sl0   = isbaspar(2)
      isleg = isbaspar(3)

C-----------------------------------------------------------------------
C Check parameters
C-----------------------------------------------------------------------
      lerr = .false.
      if (isleg.eq.0) then
         if (sl0-slz+1.gt.dim)    lerr = .true.
      elseif (isleg.eq.1) then
         if ((sl0-slz)/2+1.gt.dim)  lerr = .true.
      elseif (isleg.eq.2) then
         if ((sl0-slz)/2+1.gt.dim)  lerr = .true.
      else
         lerr = .true.
      endif
      if (dim.gt.gdim) lerr = .true.
      if (ibleg.gt.0 .and. ibleg.ne.isleg ) lerr = .true.
      if (basis.ne.2 .and. basis.ne.7 .and. slz.ne.0) lerr = .true.
      if ((basis.le.0 .or. basis.gt.5) .and. (basis.ne.7 .or. 
     +     basis.ne.9 .or. basis.ne.12) ) lerr = .true.
      if ((basis.eq.2 .or. basis.eq.12) .and. ibleg.eq.0 .and.
     +     isleg.gt.0  .and. (2*dim+1-isleg).gt.gdim  )  lerr = .true.
      message='Illegal input parameters for generation of Leg. func.'
      if (lerr) write(6,*) message 


      if(basis.eq.2.or.basis.eq.7.or.basis.eq.9.or.basis.eq.12) then

C-----------------------------------------------------------------------
C Take the initial spf from trafo (DVR.eq.Leg)
C-----------------------------------------------------------------------
         do e=1,dim
            j = e 
            if(ibleg.eq.0 .and. isleg.gt.0 ) j = 2*e-1+mod(slz+ibleg,2)
            do g=1,gdim
               psi(g,e) = trafo(j,g)
            end do
         end do

      else

C-----------------------------------------------------------------------
C Get Legendre polynomials multiplied by their DVR weights (DVR.ne.Leg)
C-----------------------------------------------------------------------
         do g=1,gdim
            do e=1,dim
               if (isleg.eq.0) then
c                  psi(g,e)=weight(g)*dlegendre(e-1,ort(g))
               else 
                  isleg1 = 2-mod(isleg+slz,2)
c                  psi(g,e)=weight(g)*dlegendre(2*e-isleg1,ort(g))
               endif
            enddo
         enddo
      
C-----------------------------------------------------------------------
C Normalize single-particle functions
C-----------------------------------------------------------------------
c         call normspf(psi,gdim,dim)
      end if

      return
      end

C***********************************************************************
C
C             MCTDH  subroutine  GENKLEG
C
C
C    Generates a set of the normalized associated
C    Legendre polynomials P_l^m for the given set of (l,m) pairs.
C
C    Calls:
C           Subroutine excitations3:  Calculates the set of the
C           most probable (l,m) pairs.
C
C    03/99 SS
C***********************************************************************

      subroutine genkleg(psi,gdim,kdim,dim,trafo,basis,ort,ipbaspar1,
     +                    ipbaspar2,isbaspar1,isbaspar2,ilog,message)

      implicit none

      integer gdim,dim,basis,slz,sl0,isleg,ibleg,ipbaspar1(*),
     +        ipbaspar2(*),isbaspar1(*),isbaspar2(*),kmin,kdim
      integer e,g,j,dj,j1,m,m1,dm,mmin,mmax,ilog,jmin
      real*8  ort(gdim)
      character*(*) message
      logical lerr,ok

      real*8 trafo(gdim,gdim,kdim)
      complex*16 psi(gdim,kdim,dim)

C-----------------------------------------------------------------------
C get primitive basis (DVR) parameters
C-----------------------------------------------------------------------
      ibleg = ipbaspar1(1)
      if (basis.eq.7) then
         kmin  = ipbaspar2(1)
      elseif (basis.eq.9) then
         kmin  = ipbaspar1(2)
      endif

C-----------------------------------------------------------------------
C get spf basis parameters
C-----------------------------------------------------------------------
      sl0   = isbaspar1(1)
      isleg = isbaspar1(2)
      slz   = isbaspar2(1)
      mmin  = isbaspar2(2)
      mmax  = isbaspar2(3)
      dm    = isbaspar2(4)

C-----------------------------------------------------------------------
C Check parameters
C-----------------------------------------------------------------------
      lerr = .false.

      if (dim.gt.gdim*kdim) lerr = .true.
      if (ibleg.gt.0 .and. ibleg.ne.isleg ) lerr = .true.
      if (basis.ne.7 .and. basis.ne.9) lerr = .true.
      if (ibleg.eq.0 .and. isleg.gt.0 .and. (2*dim+1-isleg).gt.gdim)
     +                lerr=.true.
      message='Illegal input parameters for generation of KLeg. func.'
      if (lerr) write(6,*) message 

C-----------------------------------------------------------------------
C Take the initial spf from trafo (DVR.eq.KLeg)
C-----------------------------------------------------------------------

      dj=1
      j=sl0
      m=slz

      m1=m+1-kmin
      if (ibleg.eq.0) then
         j1=j-abs(m)+1
      else
         jmin = abs(m) + mod(abs(m)+ibleg,2)
         j1 = (j-jmin)/2 + 1
      endif

      write (ilog,*)
      write (ilog,*) '*** The (j,m) pairs for the KLEG IWF ***'
      write (ilog,*)
      write (ilog,'(a,i3,a,i3)') '    j=',j,'  m=',m

         do g=1,gdim
            psi(g,m1,1) = trafo(j1,g,m1)
         enddo

      do e=2,dim
 30      continue
         call excite3(j,m,sl0,slz,dj,dm)
         call check(j,m,mmin,mmax,ibleg,isleg,gdim,ok)
         if (.not.ok) goto 30

         m1=m+1-kmin
         if (ibleg.eq.0) then
            j1=j-abs(m)+1
         else
            jmin = abs(m) + mod(abs(m)+ibleg,2)
            j1 = (j-jmin)/2 + 1
         endif

         write (ilog,'(a,i3,a,i3)') '    j=',j,'  m=',m

         do g=1,gdim
            psi(g,m1,e) = trafo(j1,g,m1)
         enddo
      enddo
      write (ilog,*)

      return
      end

C***********************************************************************

      subroutine check(j,m,mmin,mmax,ibleg,isleg,gdim,ok)

      implicit none

      integer j, m, mmin, mmax, ibleg, isleg, gdim, jmax
      logical ok

      ok = .true.
C........j >= 0
      if (j .lt. 0) then
         ok = .false.
         return
      endif
C........jmin <= j <= jmax
      if (isleg.eq.0) then
          jmax=abs(m)+gdim-1
      else
          jmax=abs(m)+2*(gdim-1)
      endif
      if (j.lt.abs(m) .or. j.gt.jmax) then
          ok = .false.
          return
      endif
C........mmax <= m <= mmin
      if (m.lt.mmin .or. m.gt.mmax) then
         ok = .false.
         return
      endif
C......---symmetry (j is even or odd)---
      if (ibleg .ne. 0) then
         if (mod(j+ibleg,2).ne.0) then
            ok = .false.
            return
         elseif (isleg.gt.0) then
            if (mod(j+isleg,2).ne.0) then
               ok = .false.
               return
            endif 
         endif
      endif


      return
      end

C***********************************************************************
      
      subroutine excite3(j,m,j0,m0,dj,dm)

      implicit none

      integer j, m, j0, m0, dj, dm, del, delj, delm, sj, sm


C     Calculate auxiliary values
      sj = 1
      if (j .lt. j0)  sj = -1
      sm = 1
      if (m .lt. m0)  sm = -1
      delj = abs(j - j0)/dj
      delm = abs(m - m0)/dm
      del = delm + delj

C     Determine next configuration

C     Try to flip spin
      if ((delm .gt. 0) .and. (sm .eq. 1)) then
         sm = -1
         goto 10
      endif

C     Try to flip 'logical spin' for j
      if ((delj .gt. 0) .and. (sj .eq. -1)) then
         sm = 1
         sj = 1
         goto 10
      endif

C     Try to move excitation from j to m
      if (delm .lt. del) then
         sm = 1
         sj = -1
         delm = delm + 1
         delj = delj -1
         goto 10
      endif

C     Choose 'logical' initial configuration for del+1
      del = del + 1
      sm = 1
      sj = -1
      delm = 0
      delj = del

 10   continue

C     Calculate new configuration
      j = j0 + delj*sj*dj
      m = m0 + delm*sm*dm

      return
      end


C***********************************************************************
C
C            MCTDH  subroutine GENK
C
C     Generates the IWF for the combined-mode containing "k" inwf.
C     The IWF of the previous dof can be a gauss- or leg-inwf.
C
C  10/99 SS
C***********************************************************************

      subroutine genk(psi,workc,isbaspar,gdim,kdim,dim)

      implicit none

      integer gdim,kdim,dim,isbaspar(*),slz,mmin,m1,g,g1
      complex*16 psi(gdim,kdim,dim),workc(gdim,dim)


C-----------------------------------------------------------------------
C get spf basis parameters
C-----------------------------------------------------------------------
      slz  = isbaspar(1)
      mmin = isbaspar(2)
      m1   = slz-mmin+1

C-----------------------------------------------------------------------
C Build the IWF from the WF of the previous dof of the combined-mode, 
C with fixed value of k=m1.
C-----------------------------------------------------------------------

      do g=1,dim
         call zeromxz(psi(1,1,g),gdim,kdim)
         do g1=1,gdim
            psi(g1,m1,g) = workc(g1,g)
         enddo
      enddo

      return   
      end


C***********************************************************************
C
C            MCTDH  subroutine PLEGRUECK
C
C   Transforms KLeg INWF from (Theta,K)-grid to (Theta,Phi)-grid for
C   the combined-mode with (PLeg,Exp)-DVRs.
C
C***********************************************************************

      subroutine plegrueck(psik,exprueck,psi,gdim1,gdim2,kdim,dim)

      integer gdim1,gdim2,kdim,dim,e,ierr
      complex*16 psi(gdim1,gdim2,dim),psik(gdim1,kdim,dim),
     +           exprueck(kdim,gdim2)

      do e=1,dim
c         call mmxxzz(psik(1,1,e),exprueck,psi(1,1,e),gdim1,kdim,gdim2)
      enddo

      call schmidtortho(psi,gdim1*gdim2,dim,ierr)

      return
      end

C***********************************************************************
C
C            MCTDH  subroutine PLEGHIN
C
C   Transforms INWF from (Theta,Phi)-grid to (Theta,K)-grid for
C   the combined-mode with (PLeg,Exp)-DVRs.
C
C***********************************************************************

      subroutine pleghin(psi,exphin,psik,gdim1,gdim2,kdim,dim)

      integer gdim1,gdim2,kdim,dim,e,ierr
      complex*16 psi(gdim1,gdim2,dim),psik(gdim1,kdim,dim),
     +           exphin(gdim2,kdim)

      do e=1,dim
c         call mmxxzz(psi(1,1,e),exphin,psik(1,1,e),gdim1,gdim2,kdim)
      enddo

      call schmidtortho(psik,gdim1*kdim,dim,ierr)

      return
      end



C***********************************************************************
C
C              MCTDH    Subroutine    VEIGEN
C
C    Calculates eigenvectors and eigenvalues of a 1D Hamiltonian
C
C***********************************************************************

      subroutine veigen(psi,gdim,dim,ort,pop,hops,s,f,m,ipot,lab,
     +                  eval,evec,work,worki,eeigen,e0eigen,
     +                  nspfdof,spfdof,lwrite,ldvrf)

      implicit none

c#include "ALL_INC.inc"
#include "global.inc"
#include "maxdim.inc"
#include "maxkoe.inc"
#include "maxsta.inc"
#include "operdef.inc"
cc#include "channels.inc"

      integer f,f1,m,e,s,ss,g,gdim,dim,pop,np,nm,nn,n,ier,ipot
      integer kk,k,h,hh,ll,ilbl,worki(gdim),nspfdof,spfdof(maxdim)
      integer sym1,sym2,nsym
      real*8  eval(gdim),evec(gdim,gdim),ort(gdim),hops(hopsdim),
     +        eeigen(dim),e0eigen,work(gdim,gdim),sum
      complex*16   psi(gdim,dim)
      character*16 lab
      logical      lwrite, ldvrf

C-----------------------------------------------------------------------
C  Initialisation and checks. 
C-----------------------------------------------------------------------
      ilbl=index(hamlab(ipot),' ')-1
      routine = 'VEIGEN (genphi1.f) (keyword "eigenf")'
      ll = index(lab,' ')-1
      if (.not.ldvrf) then
         message='The primitive basis for '//lab(1:ll)//
     +           ' must be a DVR (e.g. not FFT).'  
cc         call errormsg
      endif
      if(nodiag(ipot)) then
         message='The Hamiltonian '//hamlab(ipot)(1:ilbl)//
     +           ' must be of "usediag" type for "eigenf".'  
c         call errormsg
      endif

C-----------------------------------------------------------------------

      ss = s
      if (khzahl(m,ss,ipot) .eq. 0) ss = 1     ! in case S1&1 is set
      if (khzahl(m,ss,ipot) .eq. 0) then
         write(message,'(5a)') 'Hamiltonian ',hamlab(ipot)(1:ilbl),
     +        ' used for "eigenf" does not contain an uncorrelated',
     +        ' operator       for DOF : ',lab(1:ll)
c         call errormsg
      elseif (khzahl(m,ss,ipot) .gt. nspfdof) then
         write(message,'(5a)') 'Hamiltonian ',hamlab(ipot)(1:ilbl),
     +        ' contains more uncorrelated operators than needed.', 
     +        ' DOF : ',lab(1:ll)
c         call warnmsg1
      endif

c.....For combined modes search for the appropriate single sum. 
      h = 0
      do kk = zhun(m,ss,ipot), zhun(m,ss,ipot)+khzahl(m,ss,ipot)-1
         hh = hamilton(f,kk)
         if( htmtype(hh) .eq. 1 ) goto 100 ! skip unit operators
         if(  htmdof(hh) .ne. f ) goto 100 ! skip mode operators

C........Skip operators which are non-unit at other DOFs than f.
         do n=1,nspfdof
            f1=spfdof(n)
            if(f1.ne.f .and. htmtype(hamilton(f1,kk)).ne.1) goto 100
         enddo
         if( h .eq. 0 ) then    ! Success, found an appropriate operator
            k = kk
            h = hh
         else                   ! Hmm, found a second operator
            write(message,'(5a)') 'Hamiltonian ',
     +           hamlab(ipot)(1:ilbl),' used for "eigenf" contains',
     +           ' more than one non-unit 1D-operator for DOF : ',
     +           lab(1:ll)
            write(6,'(/,a,2i4,a,2i7,a,2i4)') 'Veigen: k, kk =', 
     +           k, kk,' h, hh =',h, hh,' , f, m =',f,m
            write(6,*) ' htmtypes (h,hh) :',htmtype(h),htmtype(hh)
c            call errormsg
         endif
 100     continue
      enddo

      if( h.eq.0) then
         write(message,'(5a)') 'Hamiltonian ',
     +        hamlab(ipot)(1:ilbl),' used for "eigenf" does not ',
     +        'contain a non-unit 1D-operator for DOF : ',lab(1:ll)
c         call errormsg
      endif

      if (htmdof(h) .ne. f) then
         write(6,*) 's,f,h,k,htmdof(h),htmmode(h),htmmuld(h): ',
     +        s,f,h,k,htmdof(h),htmmode(h),htmmuld(h)
         message=
     +   'The uncorrelated operator found is not a DOF operator.'
c         call errormsg
      endif
      if (htmtype(h) .ne. 8) then
         write(6,*) 'Error : htmtype should be 8 (real matrix).' 
         write(6,*) 's,f,k,h : ',s,f,k,h,' ,  htmtype(h) =',htmtype(h)
         message='1D-Operator for eigenf must be a '//
     +        'real symmetric matrix (e.g. not FFT)'
c         call errormsg
      endif
      
C-----------------------------------------------------------------------
C copy Hamiltonian operator to matrix evec.
C eigenvectors "evec" of the real symmetric matrix "evec" 
C-----------------------------------------------------------------------
       write(ilog,'(/,a,i5,1x,4a,i3)') 
     +     '**** VEIGEN ****  operator: h =',h,
     +     hamlab(ipot)(1:ilbl),' , DOF = ',lab(1:12),
     +     'state =', s
      call cpmxd(hops(zhop(h)),evec,gdim,gdim)

C-----------------------------------------------------------------------
C get all eigenvalues "eval" and the corresponding normalized
C eigenvectors "evec" of the real symmetric matrix "evec" 
C-----------------------------------------------------------------------

      g=gdim*gdim
      call dsyev('V','U',gdim,evec,gdim,eval,work,g,ier)

      if (ier.ne.0) then
        write (ilog,*) 'Error in subroutine VEIGEN'
        write (ilog,*) 'Error on diagonalising the Hamiltonian matrix'
        write (ilog,*) 'Error code : ',ier
        stop
      endif
      if (abs(koeff(k)-1.d0) .gt. 1.d-9 ) then
         do g = 1, gdim
            eval(g) = eval(g)*dble(koeff(k))
         enddo  
      endif

C-----------------------------------------------------------------------
C Change the phases of the eigenvectors.
C The phase convention is: Sum_i [i*evec(i)] > 0.
C-----------------------------------------------------------------------
      do e=1,min(dim+pop,gdim)
         sum=0.d0
         do g=1,gdim
            sum = sum + g*evec(g,e)
         enddo  
         if(sum.lt.0.d0) then
            do g=1,gdim
               evec(g,e) = -evec(g,e)
            enddo
         endif
      enddo
 
C-----------------------------------------------------------------------
C Zero the integer work array.
C Save eigenvalue corresponding to eigenstate 'pop'.
C-----------------------------------------------------------------------
      call zerovxi(worki,gdim)
      e0eigen=eval(pop)

C-----------------------------------------------------------------------
C When pop is negative take every second eigenfunction. 
C-----------------------------------------------------------------------
      if (pop.lt.0) then
         pop = -pop
         sym1=2
         sym2=-mod(pop,2)
      else
         sym1=1
         sym2=0
      endif

C-----------------------------------------------------------------------
C Transfer real evec data to complex psi array. The order of the spf's is:
C pop, pop-1, pop-2, pop+1, pop-3, pop-4, pop+2, ...
C-----------------------------------------------------------------------      
      nm = pop 
      np = pop - 1
      do e=1,min(dim,gdim)
         if( pop .le. 1 )  then
            nn = e
         else
            if(((mod(e,3).eq.1).or.(nm .eq. 1)) .and. np.lt.gdim) then 
               np = np + 1
               nn = np
            else
               nm = nm - 1
               nn = nm
            end if
         end if
         nsym=sym1*nn + sym2
         nsym=mod(nsym+gdim-1,gdim) + 1
         worki(nsym)=e 
         eeigen(e)=eval(nsym)
         do g=1,gdim
            psi(g,e)=evec(g,nsym)
         enddo
      enddo

C-----------------------------------------------------------------------
C Write eigenvalues in [eV] to log file
C-----------------------------------------------------------------------
      write(ilog,*)
      do e=1,min(max(sym1*dim+1,np,(gdim+4)/3),gdim)
         if(worki(e) .gt. 0 ) then
            write(ilog,'(i4,a,f16.9,a,f18.12,2a,i4,a)') 
     +      e,'. Eigenvalue',eval(e)*ev,'[eV]',eval(e),'[au]',
     +      ' ,',worki(e),'. spf'  
         else
            write(ilog,'(i4,a,f16.9,a,f18.12,a)') 
     +           e,'. Eigenvalue',eval(e)*ev,'[eV]',eval(e),'[au]'
         endif
      enddo
      write (ilog,*)
c      call myflush(ilog)

C-----------------------------------------------------------------------
C if requested, write eigenvalues and vectors to veigen file
C-----------------------------------------------------------------------
      if (lwrite) then
         write(iveig,'(a)') '**** Eigenvalue in [eV] and [au] ****'
         do e=1,gdim
            write(iveig,'(i5,2f20.12)') e,eval(e)*ev,eval(e)
         enddo

         do e=1,gdim
            write(iveig,*)
            do g=1,gdim
               write(iveig,'(2f20.12)') ort(g),evec(g,e)
            enddo
         enddo
      endif

      return
      end

C#######################################################################
C
C            MCTDH subroutine GENSPHFBR
C
C
C sets up single particle functions for FBR-spherical harmonics basis
C set.
C
C variables passed for one mode:
C  psi: single particle functions (states)
C  dim: number of spfs (states)
C  gdim: nb of basis functions
C  jsph,msph  see subroutine initfbr in init1.f
C  qnl0,qnm   see main.inc
C
C 08/97 MCH
C
C***********************************************************************

      subroutine gensphfbr(psi,gdim,dim,jsph,msph,
     +                     isbaspar,ilog,routine,message)

      implicit none

      character*80 routine, message
      integer gdim,dim,ilog
      integer istart,e,i,linf,isbaspar(*)
      integer jsph(gdim),msph(gdim),qnl0,qnm0
      complex*16 psi(gdim,dim)

C-----------------------------------------------------------------------
C get spf basis parameters
C-----------------------------------------------------------------------
      qnl0 = isbaspar(1)
      qnm0 = isbaspar(2)

c    identification of the index istart of the basis function
c    corresponding to the initial wavepacket.

      i=0
 10   i=i+1
      if (i.gt.gdim) then
        routine = 'genphi1.F  GENSPHFBR'
        message = 'Wrong quantum number for spherical harmonics spf'
c        call errormsg
      endif

      if ((jsph(i).eq.qnl0).and.(msph(i).eq.qnm0)) then
            istart=i
            goto 20
      else 
            goto 10
      endif

 20   continue 



c    initial unoccupied spf are distributed around the initial 
c    occupied spf, as symetrically as possible
c    linf = first unoccupied spf.

      linf=max(istart-dim/2,1) 

      if (linf+dim-1.gt.gdim) linf=gdim-dim+1
 
      if (linf.le.0) then
        routine = 'GENSPHFBR'
        message = 'Wrong number of spherical harmonics spf'
c        call errormsg
      endif

c output to log
            write(ilog,'(/,a,/,a,i3,",",i3,/)')
     . '*** Initial wave packet  ***** ',
     . 'FBR basis : spherical harmonics (j,mj): ',qnl0,qnm0

c  set psi(g,e) to zero.
      call  zeromxz(psi,gdim,dim)

      psi(istart,1)=(1.0d0,0.0d0)
     
      do e=2,dim
       if ( linf.eq.istart ) linf=linf+1
       psi(linf,e)=(1.0d0,0.0d0)
       linf=linf+1
      enddo
         
      return 
      end



C#######################################################################
C GENHIGHER
C
C called by gen1, gen2, rstpsi
C calls to overlap, subcomp, normphi
C
C given an initial set of functions (1 or more), a higher set (starting
C  with the 'next' function) of Schmidt orthogonalised functions are
C  generated. 
C
C psi: single particle functions
C gdim: number of grid points
C dim: number of single particle functions
C next: first undefined spf
C ort: grid positions
C
C V4.0
C#######################################################################

      subroutine genhigher (psi,gdim,dim,next,ort,basis,ipbaspar,
     +                      power,gxeq,shift)

      implicit none

      integer gdim,dim,next,e,g,e1,basis,ibleg,ipbaspar(*),power
      complex*16 psi(gdim,dim),s
      real*8 ort(gdim),norm,x,gxeq,shift

      if (basis .eq. 2  .or. basis .eq. 12) then
         ibleg=ipbaspar(2)
         elseif(basis .eq. 7 .or. basis .eq. 9) then
         ibleg=ipbaspar(1)
      endif

C-----------------------------------------------------------------------
C basis(f) : basis type for dof f. (dof=degree of freedom)
C(phifbr), el, HO, Leg, sin, FFT, exp, sphfbr, kleg, k, pleg, rHO, Leg/r
C   -1      0   1    2    3    4    5      6    7    8    9    11   12  
C-----------------------------------------------------------------------
C Schmidt-orthogonalized coordinate power multiplied Gaussians
C  1. Preceeding function is multiplied by co-ordinate (or function 
C     thereof)
C  2. Components of all preceeding functions removed
C  3. Function normalised
C  NB: shift is for periodic grids.
C-----------------------------------------------------------------------       
      do e=next,dim
         if((basis.eq.2 .or. basis.eq.7 .or. basis.eq.9
     +        .or. basis.eq.12)      .and. ibleg.eq.0 ) then
            do g=1,gdim
               psi(g,e)=cos(ort(g))*psi(g,e-1)
            enddo
         else if((basis.eq.2 .or. basis.eq.7 .or. basis.eq.9
     +            .or. basis.eq.12)      .and. ibleg.gt.0 ) then
            do g=1,gdim
               psi(g,e)=(cos(ort(g))**2)*psi(g,e-1)
            enddo
         else 
            do g=1,gdim
               x=ort(g)-gxeq
               if(shift .ne. 0.d0 ) then
                  if(abs(x) .gt. abs(x-shift)) x = x-shift
                  if(abs(x) .gt. abs(x+shift)) x = x+shift
                  if(abs(x) .gt. abs(x-shift)) x = x-shift
               endif
               x=x**power
               psi(g,e)=x*psi(g,e-1)
            enddo
         end if
         do e1=1,e-1
            call vvaxzz(psi(1,e1),psi(1,e),s,gdim)
            call xvxxzzs(s,psi(1,e1),psi(1,e),gdim)
         enddo
         call normvxz(psi(1,e),norm,gdim)
         if( norm .lt. 1.d-12 ) norm = 1.d20
         call xvixdzo(norm,psi(1,e),gdim)
      enddo

      return
      end

C#######################################################################

      subroutine excite2(f,ex,b)

#include "maxdim.inc"
c#include "ALL_INC.inc"

      integer f, ex, b(f), i, c(maxdim), fh

C     initialize c-vector      
      c(1) = ex
      do i=2,(f-1)
         c(i) = c(i-1) - b(i-1) + 1
      enddo
      c(f) = b(f) - 1
      fh = f

C     determine next configuration
 10   continue
      c(fh) = c(fh) + 1
      if (c(fh) .gt. c(fh-1)) then

C     all configurations for ex are already calculated, set
C     ex = ex+1 and start with first higher configuration
         if (fh .eq. 2) then
            ex = ex + 1
            b(1) = ex + 1
            do i=2,f
               b(i) = 1
            enddo
            return
         endif

         c(fh) = 0
         fh = fh - 1
         goto 10
      endif

C     calculate next configuration for ex
      do i=1,(f-1)
         b(i) = c(i) - c(i+1) + 1
      enddo
      b(f) = c(f) + 1

      return 
      end


C######################################################################

      integer function excite(indx,n)

      integer i,n,indx(n)

      excite=0
      do i=1,n
         excite=excite+(indx(i)-1)
      enddo

      end

C#######################################################################

      subroutine multspf(spfs,psi,indx,vgdim,gdim,ngdim,nspf,dim)

      implicit none

      integer vgdim,gdim,ngdim,n,v,g,nspf,dim,indx
      complex*16 psi(vgdim,gdim,ngdim,dim),spfs(gdim,dim)

      do n=1,ngdim
         do g=1,gdim
            do v=1,vgdim
               psi(v,g,n,nspf)=psi(v,g,n,nspf)*spfs(g,indx)
            enddo
         enddo
      enddo

      return
      end

C#######################################################################
C 
C if initally populated spf is not the first in the list of generated
C functions, swap it with the first spf.
C remaining functions are then listed: n, n-1, n+1, n-2, n+2 ....
C
C#######################################################################

      subroutine swapspf(psi,workc,inspf,gdim,dim)

      implicit none

      integer gdim,dim,inspf,e,e1,n
      complex*16 psi(gdim,dim),workc(gdim,dim)

      call cpvxz(psi(1,inspf),workc(1,1),gdim)
      e=inspf
      e1=inspf
      n=1
   10 continue
      e=e-1
      e1=e1+1
      if (e .ge. 1) then
         n=n+1
         call cpvxz(psi(1,e),workc(1,n),gdim)
      endif
      if (e1 .le. dim) then
         n=n+1
         call cpvxz(psi(1,e1),workc(1,n),gdim)
      endif
      if (n .lt. dim) go to 10

      call cpmxz(workc,psi,gdim,dim)
      
      return
      end

