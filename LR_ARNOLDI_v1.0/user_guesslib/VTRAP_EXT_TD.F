      SUBROUTINE  VTRAP_EXT_TD(time)
      USE   SHARED_DIMS
      USE   W_INTERPARTICLE
      USE   DVR_ALL
      USE   rR_hW
      USE   USR_PAR
      use interpreter !external
      IMPLICIT NONE
      logical                    ::      file_exists1,file_exists2
      integer                    ::      i,ind,J,K,heaviside
      real*8                     ::      x,x0,xs,a,s0,s2,s3 !for a spline
      real*8, intent(in)         ::      time
      real*8 t,y,z,D,GR,w,x_0,delta,b,c 
      real*8 sigma, dx, Vmax, t2, r, phase, ramp
       complex*16 compt,job_kind

           t=time

!      if(stringV.ne.'')  call init (stringV, variables, statusflag)
       statusflag='ndef'
        if(stringV.ne.'Using defaults from VTRAP_EXT_TD.F') 
     .       call init (stringV, variables, statusflag)
!        if(statusflag .ne. 'ok') then
!        STOP ( "IN V trap external")
!        end if
!      endif
c!=========================================================================
      Dim_1_2_3D: SELECT CASE (DIM_MCTDHB)

       CASE (1) ! 1D  ! 1D  ! 1D  ! 1D  ! 1D  ! 1D  ! 1D  ! 1D  ! 1D
       do i=1,NDX

          
         x         =       ort_X(i)
     
!        1D harmonic trap 
           VTRAP_EXT(I)=x*x

!        1D anharmonic trap 
!           VTRAP_EXT(I)=0.5*x*x+0.01*x*x*x*x

!        1D lattice 
!           VTRAP_EXT(I)=5.0*(cos(pi/(1.0)*x))**2
!           IF(ort_x(I).lt.(-1.5).or.ort_x(I).gt.1.5)
!     .      VTRAP_EXT(I)=500.0

!        TRAP LENZ EXP
!           x_0=1.5; delta=0.1; b=0.04; c=0.0565044
!           VTRAP_EXT(I)=c+0.05/2.0*(x*x-sqrt(b*b+(delta-2*x*x_0)*
!     .                   (delta-2*x*x_0)))



       if(statusflag.eq.'ok') then
           variablesvalues(1)=x
           variablesvalues(10)=t
           VTRAP_EXT(I)= evaluate (variablesvalues) 
       endif

       EndDO

       CASE (2) ! 2D  2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D

!        inquire(file="POTENTIAL.dat",exist=file_exists1)
!        inquire(file="x_y_potential",exist=file_exists2)

!        open(unit=56,file="POTENTIAL.dat",status='replace')
!        open(unit=57,file="x_y_potential.dat")


c         write(*,*)"This potential was used in 1D"
c         stop
         ind=1
         Do J=1,NDY
         Do I=1,NDX
         x=ort_X(I)
         y=ort_Y(J)
!       sigma=2.6d0
!       dx=0.0d0
!       Vmax=4.0d0
!       t2=1.0d0
!         a=(x+dx)*(x+dx)/(2*sigma**2)
!     .          +Vmax*t2*Exp(-(x+dx)*(x+dx)/(2*sigma**2))             
!         dx=0.1d0
!         a=30*Exp(-(y+dx)*(y+dx)/(0.02))             
!         dx=0.5d0
!         a=Exp(-(y+dx)*(y+dx)*0.3846153846153846153846153846d0)
!         VTRAP_EXT(ind)=4*y*y + 2*x*x
!       VTRAP_EXT(ind)=0.5*y*y + 0.5*x*x



        r=Dsqrt(x**2+y**2)
        Rb=4.089
        Rc=3.0

        job_kind=(0.0,-1.0)
         

      if(Real(JOB_prefac).eq.0.0) then 

   
        if(Rb.le.3.0) then

! Relaxation IN
          if(r.eq.Rb) then
              VTRAP_EXT(ind)=1.0  !change to 1
          else
             if(r.eq.Rc) then
                 VTRAP_EXT(ind)=200.0
             else
 
             VTRAP_EXT(ind)=heaviside(Rb-r)*(0.0*exp(-2.0*(r-Rb)**4)  !change to 1
     .     +200.0*exp(-0.5*(r-Rc)**4))                                !change to 1
     .+heaviside(r-Rb)*heaviside(Rc-r)*max(1.0,200.0*exp(-0.5*(r-Rc)**4))
     .)   
     .+heaviside(r-Rc)*200.0
             endif
          endif
        else

! Relaxation OUT 
          if(r.eq.Rb.or.r.eq.Rc) then
            if(r.eq.Rb) VTRAP_EXT(ind)=1.0  !change to 1.0
            if(r.eq.Rc) VTRAP_EXT(ind)=200.0  !change to 1.0
          else
 
 
             VTRAP_EXT(ind)=heaviside(Rb-r)*1.0  !change to 1.0
     .+heaviside(r-Rb)*heaviside(Rc-r)*(1.0*exp(-2.0*(r-Rb)**4) !change to 1.0
     .+200.0*exp(-0.5*(r-Rc)**4))+
     .heaviside(r-Rc)*200.0
          endif
        end if
     
!      end if    !end RELAXATION part





      else

!    Symmetric Trap forward prop
         VTRAP_EXT(ind)=heaviside(Rc-r)*(0.0*exp(-2.0*(r-Rb)**4)   !change to 1
     .+200.0*exp(-0.5*(r-Rc)**4))+
     .heaviside(r-Rc)*200.0

       if(r.eq.Rc) VTRAP_EXT(ind)=200.0

      end if


      VTRAP_EXT(ind)=0.5*(x**2+y**2)


!        2D lattice 
!           VTRAP_EXT(ind)=25.0-25.0*(cos(pi/(1.0)*(x-0.5)))**2
!     .                 *(cos(pi/(1.0)*(y-0.5)))**2

         ind=ind+1
         Enddo
         Enddo



!       close(56)
!       close(57)

c======================================================================================

       CASE (3) ! 3D ! 3D ! 3D ! 3D ! 3D ! 3D ! 3D ! 3D ! 3D ! 3D ! 3D

!         write(*,*)"This potential was used in 1D"
!         stop
!         write(6,*)" 3D  Trap is OK"
         ind=1
         Do K=1,NDZ
         Do J=1,NDY
         Do I=1,NDX
c         r=DSQRT(ort_X(i1)**2+ort_Y(j1)**2)
         gr=4.0d0
         D=1.0d0
         x=ort_X(i)
         y=ort_y(j)
         z=ort_z(k)
         VTRAP_EXT(ind)=0.5*y*y + 0.5*x*x + 0.5*z*z
       if(statusflag.eq.'ok') then
           variablesvalues(1)=x
           variablesvalues(2)=y
           variablesvalues(3)=z
           variablesvalues(10)=t
           VTRAP_EXT(ind)= evaluate (variablesvalues) 
       endif
          
           ind=ind+1
         Enddo
         Enddo
         Enddo

       END SELECT Dim_1_2_3D

       if(stringV.ne.'Using defaults from VTRAP_EXT_TD.F')
     .         call destroyfunc()
       END SUBROUTINE VTRAP_EXT_TD



c============================================
c======  User-defined functions =============
 
       function heaviside (x)
 
        integer heaviside
        real*8 x
 
        heaviside = 0
        if (x > 0) then
         heaviside = 1
        else
         heaviside = 0
        end if
 
       end function heaviside



       function ramp (t,r)
 
         real*8 t,ramp,beta,r

         beta=10.0
 
         if(r.le.9.0.and.r.ge.8.0) then

!            if(t.le.1) ramp=0
 
!            if(t.le.181) ramp=beta*t
 
            if(t.le.20) ramp=180.0
 
            if(t.ge.20.and.t.le.100) ramp=180.0-beta*(t-20)

            if(t.ge.100) ramp=0

          else
               ramp=0
          end if
 
       end function ramp


       
       function phase (x,y)

        real*8 x,y,phase,beta,r,stepx,stepy
        real*8 :: pi=3.141593

        stepx=24.0/128.0
        stepy=stepx 

        beta=1.0

       if(x.ne.0.0) then
         phase=atan(y/x)+pi/2.0
         if(x<0) then
            phase=phase+pi
         end if
       else
         if(y.ne.0.0) phase=atan(y/(x+stepx))+pi/2.0
         if(y.eq.0.0) phase=atan(y+stepy/(x+stepx))+pi/2.0
       end if
 

       end function phase
