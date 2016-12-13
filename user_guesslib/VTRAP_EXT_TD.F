      SUBROUTINE  VTRAP_EXT_TD(time)
      USE   SHARED_DIMS
      USE   W_INTERPARTICLE
      USE   DVR_ALL
      USE   rR_hW
      use interpreter !external
      USE USR_PAR
      IMPLICIT NONE
      integer                    ::      i,ind,J,K
      real*8                     ::      x,x0,xs,a,s0,s2,s3 !for a spline
      real*8, intent(in)         ::      time
      real*8 t,y,z,D,GR 
      real*8 sigma, dx, Vmax, t2
       complex*16 compt

           t=time

      IF (LPAR1) THEN
      WRITE(6,*) "VTARP L",LPAR10
      WRITE(6,*) "VTARP D",DPAR20
      WRITE(6,*) "VTARP S",SPAR10
      WRITE(6,*) "VTARP I",IPAR10
      END IF 
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
      
           VTRAP_EXT(I)=(x*x/2)
       if(statusflag.eq.'ok') then
           variablesvalues(1)=x
           variablesvalues(10)=t
           VTRAP_EXT(I)= evaluate (variablesvalues) 
       endif

       EndDO

       CASE (2) ! 2D  2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D

c         write(*,*)"This potential was used in 1D"
c         stop
         ind=1
         Do J=1,NDY
         Do I=1,NDX
         x=ort_X(I)
         y=ort_Y(J)
       sigma=2.6d0
       dx=0.0d0
       Vmax=4.0d0
       t2=1.0d0
         a=(x+dx)*(x+dx)/(2*sigma**2)
     .          +Vmax*t2*Exp(-(x+dx)*(x+dx)/(2*sigma**2))             
         dx=0.1d0
         a=30*Exp(-(y+dx)*(y+dx)/(0.02))             
         dx=0.5d0
         a=Exp(-(y+dx)*(y+dx)*0.3846153846153846153846153846d0)
         VTRAP_EXT(ind)=4*y*y + 2*x*x
       VTRAP_EXT(ind)=0.5*y*y + 0.5*x*x
       if(statusflag.eq.'ok') then
           variablesvalues(1)=x
           variablesvalues(2)=y
           variablesvalues(10)=t
           VTRAP_EXT(ind)= evaluate (variablesvalues) 
       endif
         ind=ind+1
         Enddo
         Enddo

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

