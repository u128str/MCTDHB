          SUBROUTINE GUESS_PSI(PSI)
          USE   SHARED_DIMS
          USE   W_INTERPARTICLE
          USE   DVR_ALL
          USE   rR_hW
       use interpreter
          IMPLICIT NONE
c========================================================
c        COMPLEX*16  ::  PSI(NDX*NDY*NDZ,Morb)
c        COMPLEX*16, dimension(:,:) ::   PSI
c        COMPLEX*16 ::  PSI(:,:)
       COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb)  :: PSI
c       COMPLEX*16, DIMENSION(NDX*NDY*NDZ)  :: temp1,temp2,temp3
c       REAL*8,     DIMENSION(NDX*NDY*NDZ) :: Vexternal
        COMPLEX*16, ALLOCATABLE ::  temp1(:)
        COMPLEX*16, ALLOCATABLE ::  temp2(:)
        COMPLEX*16, ALLOCATABLE ::  temp3(:)
        REAL*8,     ALLOCATABLE ::  Vexternal(:)
        REAL*8,     ALLOCATABLE ::  phase(:)

       LOGICAL :: NewStart
       REAL*8 :: time, tau,Dx,Dy
       REAL*8 x0in, gr,r,D,x,y,z,x1,y1,gr1,a1,a2,a3,a4,a5,a6,a7,a8,a9,x2
       REAL*8 b1,b2,b3,b4,b5,b6,b7,b8,b9
       COMPLEX*16 :: Z1,ZDOTC,zz1,zz2,zk
       REAL*8 :: sigma,Vmax,t2,compt,theta
       REAL*8, dimension(100) :: aa
       REAL*8 rx,ry,rz
c===========================================================
          INTEGER :: I,J,K,DM,l,ll,p
          INTEGER :: i1,j1,i2,j2,ind,ierr
          real*4   start,finish, exec_time
          REAL*8 ::   DZNRM2,p3,p4
!          call read_parse_V_W_Psi
c================================================
         allocate(phase(NDX*NDY*NDZ),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in Guess_PSI"
         allocate(temp1(NDX*NDY*NDZ),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in Guess_PSI"
         allocate(temp2(NDX*NDY*NDZ),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in Guess_PSI"
         allocate(temp3(NDX*NDY*NDZ),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in Guess_PSI"
         allocate(Vexternal(NDX*NDY*NDZ),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in Guess_PSI"
c================================================
           write(6,*) " MEM for PSI  is ok!", SIZE(PSI),Morb*NDX*NDY*NDZ
            PSI =ZERO
            Vexternal=0d0
c======================== G U E S S and External T R A P ========================
c======================== G U E S S and External T R A P ========================
c======================== G U E S S and External T R A P ========================
c======================== G U E S S and External T R A P ========================
c======================== G U E S S and External T R A P ========================
       statusflag='ndef'
         Dim_1_2_3D: SELECT CASE (DIM_MCTDHB)
         CASE (1) ! 1D  ! 1D  ! 1D  ! 1D  ! 1D  ! 1D  ! 1D  ! 1D  ! 1D 
         D=1.00d0
         DO I=1,NDX
          D=1/Sqrt(2.0d0)
          D=1.0d0
         gr  =0.04d0*99/2/2
         x0in=4d0
         r=ort_X(I)
!         r=ort_X(I)+x0in
!         aa(2)=aa(1)+weight_X(I)*DSQRT(gr/2)*2/(Exp(gr*r)+Exp(-gr*r))
!         aa(1)=aa(1)-weight_X(I)*DSQRT(gr/2)*2/(Exp(gr*r)+Exp(-gr*r))

!        1D standard gaussian
          aa(1)=Exp(-(r)**2/(2*D*D))/DSQRT(2*PI)/D
          aa(2)=r*Exp(-(r)**2/(2*D*D))/DSQRT(2*PI)/D
          aa(3)=r*r*Exp(-(r)**2/(2*D*D))/DSQRT(2*PI)/D

!        1D lattice
          aa(1)=1.0*sqrt((cos(2.0*pi/(2.0)*(r+.5)))**2)
          aa(2)=1.0*sqrt((cos(2.0*pi/(2.0)*(r+.25)))**2)
          aa(3)=1.0*sqrt((cos(2.0*pi/(2.0)*(r+.25)))**2)
          aa(4)=1.0*sqrt((cos(2.0*pi/(2.0)*(r+1.0)))**2)
          aa(5)=1.0*sqrt((cos(2.0*pi/(2.0)*(r+1.0)))**2)
          aa(6)=1.0*sqrt((cos(2.0*pi/(2.0)*(r+.3)))**2)
          aa(7)=1.0*sqrt((cos(2.0*pi/(2.0)*(r+.3)))**2)
          aa(8)=1.0*sqrt((cos(2.0*pi/(2.0)*(r+.1)))**2)
          aa(9)=1.0*sqrt((cos(2.0*pi/(2.0)*(r-.1)))**2)


!          aa(1)=Exp(-(r)**2/(2*D*D))/DSQRT(2*PI)/D
!          aa(1)=sqrt(sin(pi*100/(2.0*50.0)*(r+50.0))**2)
!          aa(2)=sin(pi*100/(2.0*50.0)*(r+50.0))**1


!        x0in=4.0d0
!        r=ort_X(I)
!        gr  =0.002d0*1999/2/2
!        gr  =1.0d0
!        a1=weight_X(I)*DSQRT(gr/2)*2/(Exp(gr*r)+Exp(-gr*r))
!        r=ort_X(I)+4.0d0
!        a2=weight_X(I)*DSQRT(gr/2)*2/(Exp(gr*r)+Exp(-gr*r))
!        r=ort_X(I)-4.0d0
!        a3=weight_X(I)*DSQRT(gr/2)*2/(Exp(gr*r)+Exp(-gr*r))
!        aa(1)=a2+a3
!        aa(2)=a2-a3

                 r=ort_X(I)
              Do j=1,Morb
        statusflag='ndef'
        if(stringPSI(j).ne.'Using Defaults from Guess_PSI.F') then
           call init (stringPSI(j), variables, statusflag)
           variablesvalues(1)=ort_X(I)
           variablesvalues(11)=ort_X(I) ! r==x
           aa(j)= evaluate (variablesvalues) 
           call destroyfunc()
        endif
              enddo


          IF(MORB.ge.1) PSI(I,1)=aa(1)
          IF(MORB.ge.2) PSI(I,2)=aa(2)
          IF(MORB.ge.3) PSI(I,3)=aa(3)
            r=ort_X(I)
            Do k=3,MORB
          aa(k)=aa(k-1)*r
           IF(MORB.ge.k) PSI(I,k)=aa(k)
            EndDo
            EndDO 
c               EndDo
c=====================================
            CALL VTRAP_EXT_TD(0d0)
!> Initial guess taken as the ground-state of a given potential Vexternal(I)
            DO I=1,NDX
           x= ort_x(i)-0
       Vexternal(I)=(x-1)**2*0.05
       Vexternal(I)=0.5d0*x*x 
c       Vexternal(I)=VTRAP_EXT(I)
            EndDO 
        
        write(6,*)"ORB_DIAG ORB_Rest ORB_Binr",
     . ORB_DIAG, ORB_Rest, ORB_Binr
       IF(ORB_DIAG) CALL Guess_Diag(PSI,Vexternal)
       IF(ORB_Rest) CALL Guess_Read_ORB(PSI)
!=============== Imprinting MOmentum =======================
! Momentum is applied        
            DO I=1,NDX
           x= ort_x(i)
           zk= 0*(0d0,2.0)*x
        if(stringMOM.ne.'Imprint defaults K from Guess_PSI.F') then
           call init (stringMOM, variables, statusflag)
           variablesvalues(1)=ort_X(I)
           zk=(0d0,1.0d0)*evaluate (variablesvalues) 
           call destroyfunc()
         endif
            DO J=1,Morb
            PSI(i,J)=PSI(i,J)*Exp(zk)
            EndDO 
            EndDO 

        

         CASE (2) ! 2D  2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D
         write(6,'(a35)')"I'm constructing 2D guess"
            CALL VTRAP_EXT_TD(0d0)
!> Initial guess is constructed by hands

 
!         open(unit=33,file="phase.dat",status='replace')

         ind=1
         Do J=1,NDY
            Do I=1,NDX
               
                phase(ind)=atan2(ort_y(j),ort_x(i))!+pi/2.0

!              if(ort_x(i).ne.0.0) then
!                if(ort_x(i)<0) then
!                 phase(ind)=phase(ind)+pi
!                end if
!              else
!                if(ort_y(j)>0.0) phase(ind)=pi
!                if(ort_y(j)<0.0) phase(ind)=0.0
!                if(ort_y(j).eq.0.0) phase(ind)=0.0
!              end if
!              phase(ind)=phase(ind)*1.0 
!              if (phase(ind).ge.2*pi) then
!                 phase(ind)=mod(phase(ind),2*pi)
!              end if
!              write(33,*) ort_x(i),ort_y(j),phase(ind)

               D=1.d0
               x=Exp(-(ort_X(i))**2/(2*D*D))/DSQRT(2*PI)/D
               x2=Exp(-(ort_X(i)+1.0)**2/(2*D*D))/DSQRT(2*PI)/D
               y=Exp(-(ort_Y(j))**2/(2*D*D))/DSQRT(2*PI)/D
               rx=ort_X(I)
               ry=ort_Y(J)
!               rx=ort_X(I)+0.25d0     ! if a big number of orbitals is used it is better to diplace
!               ry=ort_Y(J)+0.25d0     ! the orbitals in order to populate the offdiagonal 1body
                                       ! density matrix elements which are inverted
!               PSI(ind,1)=x*y

!            if(sqrt(rx**2+ry**2).ge.4.42.and
!     ..sqrt(rx**2+ry**2).le.6.6) phase(ind)=0.0
!            if(sqrt(rx**2+ry**2).ge.4.42) phase(ind)=0

!              write(33,*) ort_x(i),ort_y(j),phase(ind)
 

       p=0.0

!     Inner Part      
               PSI(ind,1)=exp(-0.2*(sqrt(rx**2+ry**2))**2
     .*(2*PI*0.25*0.25)**(-1))
     .*sqrt(rx**2+ry**2)**p*exp(ZONEI*p*phase(ind))              



               IF(Morb.ge.2) PSI(ind,2)=sqrt(rx**2+ry**2)**p
     .*exp(-(sqrt(rx**2+ry**2))**2)*exp(ZONEI*p*phase(ind))



       p3=0.0
       p4=0.0  

               IF(Morb.ge.3) PSI(ind,3)=sqrt(rx**2+ry**2)**p
     .*exp(-(sqrt(rx**2+ry**2)-0.3)**2)*exp(ZONEI*p*phase(ind))

!               IF(Morb.ge.4) PSI(ind,4)=sqrt(rx**2+ry**2)**p
!     .*exp(-(sqrt(rx**2+ry**2)-0.4)**2)*exp(ZONEI*p*phase(ind))

!               IF(Morb.ge.5) PSI(ind,5)=sqrt(rx**2+ry**2)**p
!     .*exp(-(sqrt(rx**2+ry**2)-0.45)**2)*exp(ZONEI*p*phase(ind))

!    spot state
!               PSI(ind,1)=exp(-0.2*(abs(rx-0.5))**2
!     .*(2*PI*0.25*0.25)**(-1))*exp(ZONEI*p*phase(ind))
!     .*exp(-0.2*(ry**2)
!     .*(2*PI*0.25*0.25)**(-1))
!     .*sqrt(rx**2+ry**2)**p              

!               IF(Morb.ge.2) PSI(ind,2)=exp(-0.2*(abs(rx-0.4))**2
!     .*(2*PI*0.25*0.25)**(-1))*exp(ZONEI*p*phase(ind))
!     .*exp(-0.2*(ry**2)
!     .*(2*PI*0.25*0.25)**(-1))
!     .*sqrt(rx**2+ry**2)**p              



!     Outer Part
!               IF(Morb.ge.2) PSI(ind,2)=exp(-(sqrt(rx**2+ry**2)-6.7)**2)
!     .*exp(ZONEI*p*phase(ind))
!     .*sqrt(rx**2+ry**2)**p              

!               PSI(ind,1)= exp(-2*(sqrt(rx**2+ry**2)-6)
!     .**2*(2*PI*0.25*0.25)**(-1))*exp(ZONEI*p*phase(ind))
!     .*sqrt(rx**2+ry**2)**p              

!       p3=2.0
!       p4=0.0 

!               PSI(ind,3)= exp(-2*(sqrt(rx**2+ry**2)-6.5)
!     .**2*(2*PI*0.25*0.25)**(-1))*exp(ZONEI*p3*phase(ind))
!     .*sqrt(rx**2+ry**2)**p              

!               PSI(ind,4)= exp(-2*(sqrt(rx**2+ry**2)-6.8)
!     .**2*(2*PI*0.25*0.25)**(-1))*exp(ZONEI*p4*phase(ind))
!     .*sqrt(rx**2+ry**2)**p              

!              PSI(ind,5)= exp(-2*(sqrt(rx**2+ry**2)-7.9)
!     .**2*(2*PI*0.25*0.25)**(-1))*exp(ZONEI*p*phase(ind))
!     .*sqrt(rx**2+ry**2)**p              


               PSI(ind,1)=x*y
               IF(Morb.ge.2)  PSI(ind,2)=x*y*rx**1*ry**0
               IF(Morb.ge.3)  PSI(ind,3)=x*y*rx**0*ry**1
!               IF(Morb.ge.3)  PSI(ind,3)=x*y*rx**0*ry**1
               IF(Morb.ge.4)  PSI(ind,4)=x*y*rx**1*ry**1
               IF(Morb.ge.5)  PSI(ind,5)=x*y*rx**2*ry**1
               IF(Morb.ge.6)  PSI(ind,6)=x*y*rx**1*ry**2
               IF(Morb.ge.7)  PSI(ind,7)=x*y*rx**1*ry**3
               IF(Morb.ge.8)  PSI(ind,8)=x*y*rx**3*ry**1
               IF(Morb.ge.9)  PSI(ind,9)=x*y*rx**3*ry**2
               IF(Morb.ge.10)  PSI(ind,10)=x*y*rx**2*ry**3
               IF(Morb.ge.11)  PSI(ind,11)=x*y*rx**4*ry**1
               IF(Morb.ge.12) then
                 do k=12,MORB
                       PSI(ind,k)=x*y*rx**1*ry**4
                 end do
               end if
              Do k=1,Morb
        statusflag='ndef'
        if(stringPSI(k).ne.'Using Defaults from Guess_PSI.F') then
           call init (stringPSI(k), variables, statusflag)
           variablesvalues(1)=ort_X(I)
           variablesvalues(2)=ort_Y(J)
           variablesvalues(11)=sqrt(ort_X(I)**2+ort_Y(J)**2)
           PSI(ind,k)= evaluate (variablesvalues) 
           call destroyfunc()
        endif
              enddo
               ind=ind+1
            Enddo
         Enddo
!      close(33)

!> Initial guess taken as the ground-state of a given potential Vexternal(I)
c         DO I=1,NDX*NDY
c            Vexternal(I)=DREAL(VTRAP_EXT(I))
c         EndDO 
            Vexternal=DREAL(VTRAP_EXT)
       IF(ORB_DIAG) CALL Guess_Diag(PSI,Vexternal)
       IF(ORB_Rest) CALL Guess_Read_ORB(PSI)

         ind=1
         Do J=1,NDY
            y= ort_y(J)
            Do I=1,NDX
               x= ort_x(I)
           zk= 0*(0d0,2.0)*x+0*(0d0,1.0d0)*y
        if(stringMOM.ne.'Imprint defaults K from Guess_PSI.F') then
           call init (stringMOM, variables, statusflag)
           variablesvalues(1)=ort_X(I)
           variablesvalues(2)=ort_Y(J)
           zk=(0d0,1.0d0)*evaluate (variablesvalues) 
           call destroyfunc()
         endif
                 do k=1,MORB
        PSI(ind,k)= PSI(ind,k)*Exp(zk)
                 end do
               ind=ind+1
            Enddo
         Enddo
c======================================================================================
         CASE (3) ! 3D ! 3D ! 3D ! 3D ! 3D ! 3D ! 3D ! 3D ! 3D ! 3D ! 3D
!         write(6,'(a35)')" W(r-r') in 3D - is implemented"
         write(6,'(a35)')"I'm constructing 3D guess"
         ind=1
         Do K=1,NDZ
         Do J=1,NDY
         Do I=1,NDX
               rx=ort_X(I)
               ry=ort_Y(J)
               rz=ort_Z(K)
c         r=DSQRT(ort_X(i1)**2+ort_Y(j1)**2)
         gr=4.0d0
         D=1.0d0
         x=weight_X(i)*Exp(-(ort_X(i))**2/(2*D*D))/DSQRT(2*PI)/D
         y=weight_Y(j)*Exp(-(ort_Y(j))**2/(2*D*D))/DSQRT(2*PI)/D
         z=weight_Z(k)*Exp(-(ort_Z(k))**2/(2*D*D))/DSQRT(2*PI)/D
                              PSI(Ind,1)=x*y*z
               IF(Morb.ge.2)  PSI(ind,2)=x*y*z*rx**1
               IF(Morb.ge.3)  PSI(ind,3)=x*y*z*ry**1
               IF(Morb.ge.4)  PSI(ind,4)=x*y*z*rz**1
               IF(Morb.ge.5)  PSI(ind,5)=x*y*z*rx*ry*rz
         do l=6,morb
            PSI(IND,l)=PSI(IND,1)**(l*1.d0)
         enddo
              Do l=1,Morb
        statusflag='ndef'
        if(stringPSI(l).ne.'Using Defaults from Guess_PSI.F') then
         call init (stringPSI(l), variables, statusflag)
         variablesvalues(1)=ort_X(I)
         variablesvalues(2)=ort_Y(J)
         variablesvalues(3)=ort_Z(K)
         variablesvalues(11)=sqrt(ort_X(I)**2+ort_Y(J)**2+ort_Z(K)**2)
         PSI(ind,l)= evaluate (variablesvalues) 
         call destroyfunc()
        endif
              enddo
           ind=ind+1
         Enddo
         Enddo
         Enddo
            Vexternal=DREAL(VTRAP_EXT)
       IF(ORB_DIAG) STOP 'Diagonalization 3D matrix, seriously?' !CALL Guess_Diag(PSI,Vexternal)
       IF(ORB_Rest) CALL Guess_Read_ORB(PSI)

         ind=1
       Do K=1,NDZ
        z= ort_z(K)
         Do J=1,NDY
            y= ort_y(J)
            Do I=1,NDX
               x= ort_x(I)
           zk= 0*(0d0,2.0)*x+0*(0d0,1.0d0)*y+0+0*(0d0,1.0)*z
        if(stringMOM.ne.'Imprint defaults K from Guess_PSI.F') then
           call init (stringMOM, variables, statusflag)
           variablesvalues(1)=ort_X(I)
           variablesvalues(2)=ort_Y(J)
           variablesvalues(3)=ort_Z(K)
           zk=(0d0,1.0d0)*evaluate (variablesvalues) 
           call destroyfunc()
         endif
!         print *, stringMOM
                 do l=1,MORB
       PSI(ind,l)= PSI(ind,l)*Exp(zk)
!     .                       *Exp(+1*(0d0,+3.0)*y)
!     .                       *Exp(+1*(0d0,-1.0)*z)
                 end do
               ind=ind+1
            Enddo
         Enddo
        Enddo
c======================================================================================
         END SELECT Dim_1_2_3D
c==============================================================================================
           write(6,*)" Guess is OK for N.orb=",MORB
         write(6,*) "IN PSI_Guess V",SUM(VTRAP_EXT)
         Do I=1,Morb
         Z1=ZONER/DSQRT(ABS(ZDOTC(NDX*NDY*NDZ,PSI(:,I),1,PSI(:,I),1)))
         write(6,*)I,"IN norm before renormalization =",Z1
         CALL ZSCAL(NDX*NDY*NDZ,Z1,PSI(:,I),1)
         Z1=ZONER/DSQRT(ABS(ZDOTC(NDX*NDY*NDZ,PSI(:,I),1,PSI(:,I),1)))
         write(6,*)I,"IN norm AFTER renormatization =",Z1
         ENDDO
        DM=NDX*NDY*NDZ
       IF(Morb.ne.1) call  schmidtortho(PSI,DM,Morb,ierr)
         deallocate(temp1)
         deallocate(temp2)
         deallocate(temp3)
         deallocate(Vexternal)
         deallocate(phase)
         END SUBROUTINE GUESS_PSI 
