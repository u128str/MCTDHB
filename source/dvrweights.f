C#######################################################################
C               SUBROUTINE  DVRWEIGHTS
C This routine computes the DVR weights. It is called from the 
C plot-routines showsys, showd1d, showspf and showrst as well as 
C from the geninwf/genphi1.F routines gengauss and genaleg.
C HDM 03/00, New Version: 07/02
C#######################################################################
 
      subroutine dvrweights(trafo,ort,gdim1,ggdim,weight,basis,
     +                      rpbaspar,ipbaspar,ipbaspar1,lsq,lnw)

      implicit none

      integer i,blz,ibleg,g,gdim1,ggdim,basis,ipbaspar(*),ipbaspar1(*)
      real*8  trafo(ggdim,gdim1),ort(gdim1),hoxeq,hofreq,homass,
     +        x0,breite,xx,fac,weight(gdim1),rpbaspar(*)
      logical lsq,lnw

C***********************************************************************
C basis(f) : basis type for dof f. (dof=degree of freedom)
C(phifbr), el,  HO, Leg, sin, FFT, exp, sphfbr, kleg,  k, pleg, ---
C   -1      0    1    2    3    4    5      6    7     8    9    10  
C-----------------------------------------------------------------------
C  rHO, Leg/r, --- , --- , Lagu, --- ,  --- ,  --- ,  --- ,  --- 
C   11     12    13    14    15    16     17     18     19    20    
C***********************************************************************

C-----------------------------------------------------------------------
C Special cases: no-weights; el, sphfbr and k; sin, FFT and exp. 
C-----------------------------------------------------------------------
      if( lnw .or. basis.le.0 .or. basis.eq.6 .or. basis.eq.8 ) then
         do g=1,gdim1
            weight(g) = 1.d0
         enddo
         return
      end if

      if( basis.ge.3 .and. basis.le.5 ) then
         do g=1,gdim1
            weight(g) = sqrt(ort(2)-ort(1))
         enddo
         goto 100
      end if

C-----------------------------------------------------------------------
C Determine DVR parameters and compute SQRT( 2^(m+1) * m! / (2*m+1)!! ).
C-----------------------------------------------------------------------

      if (basis .eq. 1 .or. basis .eq. 11) then
         hoxeq  = rpbaspar(1)
         hofreq = rpbaspar(2)
         homass = rpbaspar(3)
      else if (basis.eq.2 .or. basis.eq.12) then
         blz    = ipbaspar(1)
         ibleg  = ipbaspar(2)
      else if (basis.eq.7 ) then
         blz    = ipbaspar1(3)
         ibleg  = ipbaspar(1)
      else if (basis.eq.9 ) then
         blz    = ipbaspar(5)
         ibleg  = ipbaspar(1)
      else if (basis.ge.15 .and. basis.le.18) then
         x0     = rpbaspar(1)
         breite = rpbaspar(2)
      endif

      if ( basis.eq.2 .or. basis.eq.7 .or.
     +     basis.eq.9 .or. basis.eq.12 ) then
         fac = 1.d0/dble(blz+1)
         do i = 1, blz+1
            fac = fac*2*i/dble(2*i-1)
         enddo
         fac = sqrt(fac)
      endif

C-----------------------------------------------------------------------
C Determine DVR weights.
C-----------------------------------------------------------------------

      do g=1,gdim1
 
        if (basis.eq.1) then
            weight(g) = trafo(1,g)*(homass*hofreq)**(-0.25d0)
     +           *exp(0.5d0*homass*hofreq*(ort(g)-hoxeq)**2)
               
c            weight(g) = trafo(1,g)*(homass*hofreq)**(-0.25d0)
c     +           *exp(0.5d0*homass*hofreq*(ort(g)-0.0d0*hoxeq)**2)
c         write(6,*)" STR dvrweights",hoxeq,ort(g),weight(g)

         elseif (basis.eq.11) then
            weight(g) = trafo(1,g)*(homass*hofreq)**(-0.75d0)
     +           *exp(0.5d0*homass*hofreq*(ort(g)-hoxeq)**2)
     +            /(ort(g)-hoxeq)
 
         elseif (basis.eq.2 .or. basis.eq.7 .or.
     +           basis.eq.9 .or. basis.eq.12) then
            if (blz.eq.0) then
               weight(g) = trafo(1,g)*fac
            elseif (blz.gt.0) then
               weight(g) = trafo(1,g)*fac/sin(ort(g))**blz
            else
               write(6,*) ' DVRweights, wrong blz: ', blz
               stop
            endif
            if((ibleg.gt.0) .and. (mod(ibleg+blz,2).eq.1))
     +         weight(g) = weight(g)/(sqrt(dble(2*blz+3))*cos(ort(g)))

         elseif (basis.eq.15 ) then
            xx = (ort(g)-x0)/breite
            weight(g) = trafo(1,g)*exp(0.5d0*xx)*sqrt(breite/xx)

         elseif (basis.eq.16 ) then
            xx = (ort(g)-x0)/breite
            weight(g) = trafo(1,g)*exp(0.5d0*xx)*sqrt(2.d0*breite)/xx

         elseif (basis.eq.17 ) then
            xx = (ort(g)-x0)/breite
            weight(g)=trafo(1,g)*exp(0.5d0*xx)*sqrt(6.d0*breite/xx)/xx

         elseif (basis.eq.18 ) then
            xx = (ort(g)-x0)/breite
            weight(g) = trafo(1,g)*exp(0.5d0*xx)*sqrt(24.d0*breite)
     +                              /(xx*xx)


         else
            write(6,*) ' DVRweights: Cannot do basis no:', basis
            stop
         endif

      enddo

 100  if(lsq) then
C.....If the weights and not the square roots of them are needed.
         do g = 1, gdim1
            weight(g) = weight(g)**2
         enddo
       endif

      return
      end

C#######################################################################



