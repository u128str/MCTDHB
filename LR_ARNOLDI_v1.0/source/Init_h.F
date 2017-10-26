      subroutine init_h(xi,xf,igdim,hdvr,weight,
     . ort,ort_p,dif1mat,DVRMETHOD)
      USE   W_INTERPARTICLE
      Implicit NONE
      integer igdim,ndim,j,kdim,i,kmin
C============For Getgauss
      integer hotyp,ierr
      real*8, ALLOCATABLE :: trafo(:,:),dif2mat(:,:)
      real*8 hofm,hoxeq,x,dx
      real*8  weight(igdim),ort(igdim),ort_p(igdim)
      real*8  rpbaspar(3), rsbaspar(5),PI1,GGauss
      real*8, ALLOCATABLE ::  workr(:)
      COMPLEX*16, ALLOCATABLE ::  psi(:,:)
      COMPLEX*16, ALLOCATABLE ::  exphin(:,:),exprueck(:,:)
c==============For SCF=============
      real*8 hdvr(igdim,igdim),dif1mat(igdim,igdim)

      integer DVRMETHOD,indimhalf
 
      real*8  DDOT,DNRM2,DZNRM2


      character*20 message
      character*16 lname 
      character*26 lnamenew,lnamerest  
c============================================
      real*8 xmassinp,xfreqinp,xeqlinp,xi,xf,xmass,l
      real*8 xall,wscl,xeql,xmass1,xfreq1,xfreq
      parameter(xmassinp=1.0d0,xfreqinp=1.0d0,xeqlinp=0.0d0)
c      parameter(xmassinp=0.5d0,xfreqinp=1.0d0,xeqlinp=0.0d0)
      parameter        (PI1 = 3.141592653589793238462643d0)
      external initho,gengauss,schmidtortho
c==============================================================
        IF(DVRMETHOD/=4) THEN
        allocate(trafo(igdim,igdim),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in Init_h: trafo"
        allocate(dif2mat(igdim,igdim),stat=ierr)
!        if(ierr /= 0)write(*,*)"allocation error in Init_h: dif2mat"
!        allocate(dif1mat(igdim,igdim),stat=ierr)
         ENDIF
c=========DVRMETHOD:  1         3       4     5
c=========meaning  :  HO DVR    SinDVR  FFT  ExpDVR
c        xi=-3.0d0*PI
c        xf= 3.0d0*PI
c        xi=-10.0d0*0.75d0
c        xf= 10.0d0*0.75d0
c        xi=-20.0d0
c        xf= 20.0d0
c        xf=6 .1092682926829269d0
c        xi=-7.5d0
c        xf= 7.5d0
c        xi=-15.0d0
c        xf= 15.0d0
c          l=1.0d0
c        xi=-3.0d0*PI*l
c        xf= 3.0d0*PI*l
c        xi=Time_xint
c        xf=Time_xfnl
c        xi=-1.5d0*PI/10
c        xf= 1.5d0*PI/10
c===============================================================
       indimhalf=(igdim+1)/2
C=====================Construction of DVR Matrix================C
        ndim=igdim
        IF(DVRMETHOD.EQ.1) THEN
        allocate(psi(igdim,igdim),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in HO Init_h"
c=======HO DVR parametrizations ======================================
       wscl=1.0D0
       xmass=xmassinp/wscl
       xfreq=1.0D0*xfreqinp/wscl
       xeql=xeqlinp

c        write(6,*)"DVR Method Harmonic Oscillator"
         hofm=xmass*xfreq
         hoxeq=xeql
         dx=0.0D0
         dx=xf-xi
         hotyp=1
         call initho(trafo,ort,dif2mat,dif1mat,ndim
     +                   ,hofm,hoxeq,dx,hotyp,message)

c                      rpbaspar(1)=xeql
c                      rpbaspar(2)=xfreq
c                      rpbaspar(3)=xmass
c                      rsbaspar(1)=xeql
c                      rsbaspar(2)=0.0D0
c                      rsbaspar(3)=dsqrt(2.0D0*xfreq*xmass)
c                      rsbaspar(4)=0.0D0
c                      rsbaspar(5)=xmass
c======================Slightly modified
                      rpbaspar(1)=hoxeq
                      rpbaspar(2)=xfreq
                      rpbaspar(3)=hofm
                      rsbaspar(1)=hoxeq
                      rsbaspar(2)=0.0d0
                      rsbaspar(3)=dsqrt(2.0D0*xfreq*xmass)
                      rsbaspar(4)=0.0D0
                      rsbaspar(5)=xmass

                  call gengauss(psi,ndim,ndim,weight,
     +                    trafo,ort,1,
     +                    rpbaspar,0,0,
     +                    rsbaspar,0,ndim)
               do i=1,ndim
c               do j=1,ndim
c               hofunk(i,j)=DREAL(psi(i,j))
c               enddo
c               write(6,*)i," ",weight(i),"   ",weight(ndim-i+1) 
               weight(i)=weight(i)*DSQRT(DSQRT(PI)) 
               enddo
C=================Renormalizations ORT,Weights
       if(xeql.eq.0.0d0) then
c       write(6,*)" Symmetrization of the ORT and Weights Matrixes is ON"
       do i=1,indimhalf
       ort(i)=-1.0D0*ort(ndim-i+1) 
       weight(i)=weight(ndim-i+1) 
       enddo
       endif 
      write(6,'(a50,f20.16,a3,f20.16,a8,i10)')
     ."DVR Method Harmonic Oscillator on interval  [",
     .ort(1),":",ort(ndim),"] DIM ",ndim
      ELSEIF(DVRMETHOD.EQ.5) THEN
c        write(6,*)"DVR Method Exponential DVR"
        allocate(exphin(igdim,igdim),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in Exp Init_h"
        allocate(exprueck(igdim,igdim),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in Exp Init_h"
c========== Interval adjastment for periodic boundaries==========
      xall=xf-xi
      xi=xi
      xf=xf-(xall)/(ndim)
        kmin=1
        kdim=1
      CALL initexp(trafo,ort,dif2mat,dif1mat,exphin,exprueck,ndim,
     +                   xi,xf,kmin,kdim)
c        write(6,*)"DVR Method Exponential DVR from",xi," to ",xf
                      rpbaspar(1)=0.d0
                      rpbaspar(2)=0.d0
                      rpbaspar(3)=0.d0
                      rsbaspar(1)=0.d0
                      rsbaspar(2)=0.D0
                      rsbaspar(3)=0.d0
                      rsbaspar(4)=0.D0
                      rsbaspar(5)=0.d0
      CALL dvrweights(trafo,ort,ndim,ndim,weight,5,
     +                rpbaspar,0,0,.false.,.false.)
      write(6,'(a50,f20.16,a3,f20.16,a8,i10)')
     ."DVR Method Exponential DVR     on interval  [",
     .ort(1),":",ort(ndim),") DIM ",ndim
      ELSEIF(DVRMETHOD.EQ.3) THEN
        allocate(workr(igdim*igdim),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in SIN Init_h"
c========== Interval adjastment for periodic boundaries==Particle-in-a-box ========
c       xi=-0.5d0+1.d0/(ndim+1)
c       xf=+0.5d0-1.d0/(ndim+1)
       xall=xf-xi
      xi=xi+(xall)/(ndim+1)
      xf=xf-(xall)/(ndim+1)
      CALL initsin(trafo,ort,dif2mat,dif1mat,ndim,xi,xf,workr)
c        write(6,*)"DVR Method SIN  DVR from",xi," to ",xf
                      rpbaspar(1)=0.d0
                      rpbaspar(2)=0.d0
                      rpbaspar(3)=0.d0
                      rsbaspar(1)=0.d0
                      rsbaspar(2)=0.d0
                      rsbaspar(3)=0.d0
                      rsbaspar(4)=0.d0
                      rsbaspar(5)=0.d0
      CALL dvrweights(trafo,ort,ndim,ndim,weight,3,
     +                rpbaspar,0,0,.false.,.false.)
      write(6,'(a50,f20.16,a3,f20.16,a8,i10)')
     ."DVR Method          SIN DVR    on interval  [",
     .ort(1),":",ort(ndim),"] DIM ",ndim
      ELSEIF(DVRMETHOD.EQ.4) THEN
      call getpos(xi,xf,igdim,ort)
      call getmom(xi,xf,igdim,ort_p)
C=================Renormalizations ORT,Weights
         do i=1,ndim
         weight(i) = sqrt(ort(2)-ort(1))
         enddo

      write(6,'(a50,f20.16,a3,f20.16,a8,i10)')
     ."DVR Method          FFT DVR    on interval  [",
     .ort(1),":",ort(ndim),"] DIM ",ndim
      return
c      stop
      ENDIF

C=====================Construction of Hamiltonian Matrix================C
         xmass1=xmassinp
         xfreq1=xfreqinp
         xmass1=Time_mass
         Do i=1,ndim
         Do j=1,ndim
         hdvr(i,j)=(-1.0D0/(2.0D0*xmass1))*dif2mat(i,j)
         enddo
         enddo 
C============================= Symmetrization =========================
c         Do i=1,ndim
c         Do j=i,ndim
c         hdvr(j,i)=hdvr(i,j)
c         enddo
c         enddo
c         Do i=1,ndim
c       hamORG(i,i)=hamORG(i,i)+ pt(ort(i))
c         Do i=1,igdim
c         write(6,*)" HDVR",(hdvr(i,J),j=1,igdim)
c         enddo
c         write(6,*)" ORT",ort
c         write(6,*)" HDVR",weight
c===========================================
c         stop
!          write(6,*)"D/DX",dif1mat(1,:)
!           pause 
        deallocate(trafo,dif2mat)
        IF(DVRMETHOD.EQ.1) deallocate(psi)
        IF(DVRMETHOD.EQ.3) deallocate(workr)
        IF(DVRMETHOD.EQ.5) deallocate(exphin,exprueck)
         return
         end subroutine init_h
