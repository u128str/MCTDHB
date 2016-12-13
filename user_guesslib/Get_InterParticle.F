!  Subroutine to provide All information on interparticle interaction to module W_INTERPARTICLE
!   W(r-r')= SUM_i (i=1,NtVec) F_i F_ij
          SUBROUTINE Get_InterPart(WSLTime)
          USE SHARED_DIMS
          USE DVR_ALL
          USE W_INTERPARTICLE
          use interpreter !external

          IMPLICIT NONE
          REAL*8 :: VINTERP,r,D,fpot,xerr
          REAL*8 :: x,y,z,dx,dy,dz,xall,yall,zall,Sigma
          REAL*8 :: x1,y1,z1,x2,y2,z2
          INTEGER i1,j1,i2,j2,ind,ndim,m
          INTEGER :: I,J,K
          EXTERNAL DCOPY,VINTERP
          REAL*8  :: WSLTime,t  ! For time dependent inter.particle interactions
c==============================================================================
c===============!!!!!!!!!!! IN 1 D and W(r-r') it FASTER to use Wxx_TYPE=2 instead of Wxx_TYPE=1
c=============== Wxx_TYPE=2 works only with equidistant grids!!!!!!!!!!


c=============== Select CASE of  W(r-r') in use:
C From namelist input.inp         Wxx_TYPE=0
!           write(6,*)" Wxx_TYPE from Namelist=",Wxx_TYPE
c=============== Wxx_TYPE=0:  W(r,r')=Delta(x-x')
c=============== Wxx_TYPE=1:  W(r,r')=W(x-x')*W(y-y')*W(z-z')
c=============== Wxx_TYPE=2:  W(r,r')=W(r-r') ONLY FOR Equidistant grids!!!!!!
c=============== Wxx_TYPE=3:  W(r,r')=W(r,r') one need to store matrix W(i,j) Dim(W)=NDX*NDY*NDZ * NDX*NDY*NDZ
          statusflag='ndef'
        if(stringW.ne.'Using Defaults from Get_InterParticle.F') 
     .               call init (stringW, variables, statusflag)

         WxxTYPE: SELECT CASE (Wxx_TYPE)
c============================== Delta Function Case
         CASE (0) !Delta(x-x')
!         write(6,*)" Case of Delta-Function Interaction W=Delta(R-R')*lambda(t)"
           t=WSLTime
           xlambda0=xlambda_0 !*(1.0d0-0.5*sin(t*1))
       if(statusflag.eq.'ok') then
           variablesvalues(10)=t
            xlambda0=xlambda_0*evaluate (variablesvalues) 
       endif
c============================== Case Separable interparticel potential WOp_X*WOp_Y*WOp_Z
         CASE (1) !W(x-x')*W(y-y')*W(z-z')
            t=WSLTime
!          write(6,*)"    Case of SEPARABLE W=W(R-R')=WOp_X*WOp_Y*WOp_Z"
!          write(6,*)" Three WOP_X operators are needed to evaluate W_SL"
         D=0.1d0
         D=0.05d0
         Do I=1,NDX
         Do J=1,NDX
         x1=ort_X(J)
         x2=ort_X(I)
         r=(ort_X(J)-ort_X(I))
         WOp_X(I,J)=Exp(-r**2/(2*D*D))/(DSQRT(2*PI)*D)
         WOp_X(I,J)=(x1-x2)**2 
!     & *(1.0d0-0.5*sin(t*1))
       if(statusflag.eq.'ok') then
           variablesvalues(11)=r
           variablesvalues(10)=t
           WOp_X(I,J)= evaluate (variablesvalues) 
       endif
         Enddo
         Enddo
         IF(DIM_MCTDHB.ge.2) THEN
        if(stringW.ne.'Using Defaults from Get_InterParticle.F') 
     .     STOP 'Change W in Get_InterparticleF'
         Do I=1,NDY
         Do J=1,NDY
          y1=ort_Y(J)
          y2=ort_Y(I)
         r=(ort_Y(I)-ort_Y(J))
         fpot=(Exp(-200*r**2))*(200/PI)
         fpot=Exp(-r**2/(2*D*D))/(DSQRT(2*PI)*D)
         WOp_Y(I,J)=Exp(-r**2/(2*D*D))/(DSQRT(2*PI)*D)
         WOp_Y(I,J)=(y1-y2)**2 
         Enddo
         Enddo
         ENDIF
         IF(DIM_MCTDHB.eq.3) THEN
        if(stringW.ne.'Using Defaults from Get_InterParticle.F') 
     .   STOP "Change W in Get_Interparticle.F"
         Do I=1,NDZ
         Do J=1,NDZ
          z1=ort_Z(J)
          z2=ort_Z(I)
         r=(ort_Z(I)-ort_Z(J))
         WOp_Z(I,J)=Exp(-r**2/(2*D*D))/(DSQRT(2*PI)*D)
         WOp_Z(I,J)=(z1-z2)**2 
         Enddo
         Enddo
         ENDIF
c============================== Case UnSeparable interparticel potential W(r-r')
         CASE (2) !W(r,r')=W(r-r')
          write(6,*)DIM_MCTDHB,
     ."D- Case of W=W(R-R') only one complex vector   "
        write(6,*)"W2xx(2*NDX*NDY*NDZ-1) one needs to evaluate W_SL  "
        write(6,*)" THIS TRUE ONLY FOR EQUIDISTANT DVR-GRID!!!"
       xerr=ort_X(NDX)-ort_X(NDX-1)-ort_X(NDX-1)+ort_X(NDX-2)
      IF(ABS(xerr).ge.0.0000001d0) then
      write(6,*)" Your X-Grid is NON equidistant!!!- Abort"
      stop
      endif
        IF(DIM_MCTDHB.eq.2) then 
        xerr=ort_Y(NDY)-ort_Y(NDY-1)-ort_Y(NDY-1)+ort_Y(NDY-2)
        IF(ABS(xerr).ge.0.0000001d0) then
        write(6,*)" Your Y-Grid is NON equidistant!!!- Abort"
        stop
       endif
       endif

         ind=1
         D=1d0/20
         D=0.05d0
         Do I=1,NDZ
         Do J=1,NDY
         Do K=1,NDX
         IF(DIM_MCTDHB.eq.1) r=ort_X(K)-ort_X(1)
         IF(DIM_MCTDHB.eq.2) 
     .   r=DSQRT((ort_X(K)-ort_X(1))**2+(ort_Y(J)-ort_Y(1))**2)
         IF(DIM_MCTDHB.eq.3) r=DSQRT((ort_X(K)-ort_X(1))**2+
     .   (ort_Y(J)-ort_Y(1))**2+(ort_Z(I)-ort_Z(1))**2)
         fpot=Exp(-r**2/(2*D*D))/DSQRT(2*PI)/D
         fpot=(Exp(-200*r**2))*(200/PI)
         fpot=Exp(-r**2/(2*D*D))/(DSQRT(2*PI)*D)**DIM_MCTDHB
         fpot= r**2
         W2xx(NDX*NDY*NDZ-1+ind)=fpot*ZONER
       if(statusflag.eq.'ok') then
           variablesvalues(11)=r
           variablesvalues(10)=t
           W2xx(NDX*NDY*NDZ-1+ind)= evaluate (variablesvalues) 
       endif
         if(ind.gt.1)  then
         W2xx(NDX*NDY*NDZ-ind+1)=fpot*ZONER
       if(statusflag.eq.'ok') then
           variablesvalues(11)=r
           variablesvalues(10)=t
           W2xx(NDX*NDY*NDZ-ind+1)= evaluate (variablesvalues) 
       endif
           endif
         ind=ind+1
         Enddo
         Enddo
         Enddo
         ndim=NDX*NDY*NDZ
      write(6,'(a3,20F10.5)')"W2xx",(DREAL(W2xx(I)),I=ndim-5,ndim+5)
c============================== General Case UnSeparable interparticel potential W(r,r')
         CASE (3)
         Dim_1_2_3D: SELECT CASE (DIM_MCTDHB)
         CASE (1)
         write(6,'(a30)')" W(x-x') in 1D"
         Do I=1,NDX
         Do J=1,NDX
         D=0.1d0
         D=0.05d0
         r=DSQRT((ort_X(I)-ort_X(J))**2)
         fpot=Exp(-r**2/(2*D*D))/DSQRT(2*PI)/D
         fpot=(Exp(-200*r**2))*DSQRT(200/PI)
         fpot=Exp(-r**2/(2*D*D))/(DSQRT(2*PI)*D)**DIM_MCTDHB
         fpot=r**2
         Wxx(I,J)=ZONER*fpot
       if(statusflag.eq.'ok') then
           variablesvalues(11)=r
           variablesvalues(10)=t
           Wxx(I,J)= evaluate (variablesvalues) 
       endif
         Enddo
         Enddo
         DO J=NDX/2-5,NDX/2+5
         write(6,'(a3,20F10.5)')" ",(DREAL(Wxx(J,I)),I=NDX/2-5,NDX/2+5)
         ENDDO
          write(6,*) "Integral V(x-x')=",SUM(Wxx)
         CASE (2)
c         write(6,'(a40)')" they sould be NDX*NDY:NDX*NDY !!!"
c         stop
         Do I=1,NDX*NDY
         IF(I-INT(I/NDX)*NDX.eq.0) then
         i1=NDX
         j1=INT(I/NDX)
         else
         i1=I-INT(I/NDX)*NDX
         j1=INT(I/NDX)+1
         endif
         Do J=1,NDX*NDY
         IF(J-INT(J/NDX)*NDX.eq.0) then
         i2=NDX
         j2=INT(J/NDX)
         else
         i2=J-INT(J/NDX)*NDX
         j2=INT(J/NDX)+1
         endif
         D=0.1d0
         D=0.05d0
         r=sqrt((ort_X(i1)-ort_X(i2))**2+(ort_Y(j1)-ort_Y(j2))**2)
         Wxx(I,J)=ZONER*Exp(-r/(2*D*D))/(2*PI)/D/D
         Wxx(I,J)=r**2
       if(statusflag.eq.'ok') then
           variablesvalues(11)=r
           variablesvalues(10)=t
           Wxx(I,J)= evaluate (variablesvalues) 
       endif
         Enddo
         Enddo
c         pause
          DO J=NDX*NDY/2-5,NDX*NDY/2+5
       write(6,'(a3,20F10.5)')" ",
     . (DREAL(Wxx(J,I)),I=NDX*NDY/2-5,NDX*NDY/2+5)
          ENDDO
          write(6,*) "Integral V(x-x')=",SUM(Wxx)
         CASE (3)
         write(6,'(a40)')" W(r-r') in 3D - is not implemented"
         stop
         END SELECT Dim_1_2_3D
CCCCCCCCCCCCCCCCCCCCCCC FOR IMEST
      CASE (4) !W(r,r')=W(r-r') FOR FFT_WSL
!          write(6,*) "In W t=",wsltime
!          write(6,*)DIM_MCTDHB,
!     ."D- Case of W=W(R-R') only one complex vector, IMEST ALGORITHM"
!        write(6,*)"W3xx and W3xxFFT(NDX*NDY*NDZ) "
      IF(DIM_MCTDHB.ge.1) then
      xerr=ort_X(NDX)-ort_X(NDX-1)-ort_X(NDX-1)+ort_X(NDX-2)
      IF(ABS(xerr).ge.0.0000001d0) then
      write(6,*)" Your X-Grid is NON equidistant!!!- Abort"
      stop
      endif
      endif
      IF(DIM_MCTDHB.ge.2) then
      xerr=ort_Y(NDY)-ort_Y(NDY-1)-ort_Y(NDY-1)+ort_Y(NDY-2)
      IF(ABS(xerr).ge.0.0000001d0) then
      write(6,*)" Your Y-Grid is NON equidistant!!!- Abort"
      stop
      endif
      endif
      IF(DIM_MCTDHB.ge.3) then
      xerr=ort_Z(NDZ)-ort_Z(NDZ-1)-ort_Z(NDZ-1)+ort_Z(NDZ-2)
      IF(ABS(xerr).ge.0.0000001d0) then
      write(6,*)" Your Z-Grid is NON equidistant!!!- Abort"
      stop
      endif
      endif

       dx=0d0
       dy=0d0
       dz=0d0

                          dx=ort_x(2)-ort_x(1)
      IF(DIM_MCTDHB.ge.2) dy=ort_y(2)-ort_y(1)
      IF(DIM_MCTDHB.ge.3) dz=ort_z(2)-ort_z(1)

CCCCCCCCCCCCCCCCCC
CCCCCCC definition of spatial grid
      xall=ort_x(NDX)-ort_x(1)+dx
      yall=ort_y(NDY)-ort_y(1)+dy
      zall=ort_z(NDZ)-ort_z(1)+dz
      do m=1,NDX*NDY*NDZ
         call get_ijk_from_m(m,ndx,ndy,i,j,k)
         x = (i-1)*dx-(xall/2.d0)
         y = (j-1)*dy-(yall/2.d0)
         z = (k-1)*dz-(zall/2.d0)
         sigma=3.0d0
         w3xx(m) = xlambda0/((DSQRT(2.d0*pi*sigma**2))**DIM_MCTDHB)*
     &      DEXP(-1.d0*(x**2+y**2+z**2)/(2.d0*sigma**2))
               r=sqrt(x**2+y**2+z**2)
           w3xx(m) = xlambda0*exp(-((r/D)**6)/2)
            t=WSLTime
           w3xx(m) = xlambda0*r**2
!       *(1-0.1d0*0.2d0*sin(t)*sin(t))
       if(statusflag.eq.'ok') then
           variablesvalues(11)=r
           variablesvalues(10)=t
           w3xx(m)= xlambda0*evaluate (variablesvalues) 
       endif
      enddo
      call computevtilde(w3xxfft,w3xx)
         END SELECT WxxTYPE
        if(stringW.ne.'Using Defaults from Get_InterParticle.F') 
     .          call destroyfunc()
!           write(6,*)" INteraction time =",WSLTime
         END  SUBROUTINE Get_InterPart
