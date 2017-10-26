
       subroutine User_OP(time,PSI,exp_val_U_op)
       USE   PASS_ARG
       USE   SHARED_DIMS
       USE   rR_hW
       USE   CI_All
!       USE   PROP_MB
       USE   W_INTERPARTICLE
       USE   DVR_ALL
       USE   Parallel_Orb
        IMPLICIT NONE
       integer ::  i,j,k,Iorb,Jorb,ind,icnt,Korb,Lorb,P
       complex*16 :: UOp_ij(Morb,Morb,10)
       complex*16 :: W,E,XSQ,YSQ,ZSQ
       real*8 :: exp_val_U_op(20),x,y,z,t,r,time
       COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb)  :: PSI
       COMPLEX*16, DIMENSION(NDX*NDY*NDZ,10)  :: UOp_x_y_z 
       integer :: cK,cJ,cL,cI,L2B
       complex*16 :: rho_jk(Morb,Morb)
       complex*16 :: rho_ijkl(Morb,Morb,Morb,Morb)

!          srtingOP(1)='x'
!          stringOP(2)='x**2'
!        statusflag='ndef'
!        if(stringPSI(l).ne.'Using Defaults from Guess_PSI.F') then
!         call init (stringPSI(l), variables, statusflag)
!         variablesvalues(1)=ort_X(I)
!         variablesvalues(2)=ort_Y(J)
!         variablesvalues(3)=ort_Z(K)
!         variablesvalues(11)=sqrt(ort_X(I)**2+ort_Y(J)**2+ort_Z(K)**2)
!          op1(ind)= evaluate (variablesvalues) 
!         call destroyfunc()
!         endif


         t=time
         ind=1
         Do K=1,NDZ
         Do J=1,NDY
         Do I=1,NDX
               x=ort_X(I)
               y=ort_Y(J)
               z=ort_Z(K)
         UOp_x_y_z(ind,1)=x
         UOp_x_y_z(ind,2)=y
         UOp_x_y_z(ind,3)=z
         UOp_x_y_z(ind,4)=x*x
         UOp_x_y_z(ind,5)=y*y
         UOp_x_y_z(ind,6)=z*z
         ind=ind+1
         Enddo
         Enddo
         Enddo

        DO Iorb=1,Morb
           DO Jorb=Iorb,Morb
              UOp_ij(Iorb,Jorb,1)=
     .  SUM(PSI(:,Jorb)*UOp_x_y_z(:,1)*dconjg(PSI(:,Iorb)))
              UOp_ij(Iorb,Jorb,2)=
     .  SUM(PSI(:,Jorb)*UOp_x_y_z(:,2)*dconjg(PSI(:,Iorb)))
              UOp_ij(Iorb,Jorb,3)=
     .  SUM(PSI(:,Jorb)*UOp_x_y_z(:,3)*dconjg(PSI(:,Iorb)))
              UOp_ij(Iorb,Jorb,4)=
     .  SUM(PSI(:,Jorb)*UOp_x_y_z(:,4)*dconjg(PSI(:,Iorb)))
              UOp_ij(Iorb,Jorb,5)=
     .  SUM(PSI(:,Jorb)*UOp_x_y_z(:,5)*dconjg(PSI(:,Iorb)))
              UOp_ij(Iorb,Jorb,6)=
     .  SUM(PSI(:,Jorb)*UOp_x_y_z(:,6)*dconjg(PSI(:,Iorb)))
              UOp_ij(Iorb,Jorb,7)=
     .  SUM(PSI(:,Jorb)*VTRAP_EXT(:)*dconjg(PSI(:,Iorb)))
              UOp_ij(Jorb,Iorb,1)=Conjg(UOp_ij(Iorb,Jorb,1))    
              UOp_ij(Jorb,Iorb,2)=Conjg(UOp_ij(Iorb,Jorb,2))    
              UOp_ij(Jorb,Iorb,3)=Conjg(UOp_ij(Iorb,Jorb,3))    
              UOp_ij(Jorb,Iorb,4)=Conjg(UOp_ij(Iorb,Jorb,4))    
              UOp_ij(Jorb,Iorb,5)=Conjg(UOp_ij(Iorb,Jorb,5))    
              UOp_ij(Jorb,Iorb,6)=Conjg(UOp_ij(Iorb,Jorb,6))    
              UOp_ij(Jorb,Iorb,7)=Conjg(UOp_ij(Iorb,Jorb,7))    
           EndDO   
        EndDO 

        Rho_ijkl=0.d0
        CALL Get_RhoIJKL(Rho_IJKL)


        XSQ=0.d0
        YSQ=0.d0
        ZSQ=0.d0
        DO CI=1,Morb
          DO CJ=1,Morb
            DO CK=1,Morb
              DO CL=1,Morb
           XSQ=XSQ+Rho_ijkl(CI,CJ,CK,CL)*UOp_ij(CI,CK,1)
     .             *UOp_ij(CJ,CL,1)
           YSQ=YSQ+Rho_ijkl(CI,CJ,CK,CL)*UOp_ij(CI,CK,2)
     .             *UOp_ij(CJ,CL,2)
           ZSQ=ZSQ+Rho_ijkl(CI,CJ,CK,CL)*UOp_ij(CI,CK,3)
     .             *UOp_ij(CJ,CL,3)
!           WRITE(9999,*) cI,cJ,cK,cL,Rho_ijkl(CI,CJ,CK,CL)
              END DO
            END DO
          END DO
       END DO
 
!       WRITE(9998,*) XSQ,Npar

       W=Zero
       E=SUM(AllHIJ*AllZRIJ)
       W=SUM(ZRIJKL*WIJKL*(DREAL(RESCALE_2B)))
     . +SUM(Conjg(ZRIJKL)*Conjg(WIJKL)*(DIMAG(RESCALE_2B)))

        exp_val_U_op=0.d0
        exp_val_U_op(1)=Real(SUM(UOp_ij(:,:,1)*AllZRIJ)) !x
        exp_val_U_op(2)=Real(SUM(UOp_ij(:,:,2)*AllZRIJ)) !y
        exp_val_U_op(3)=Real(SUM(UOp_ij(:,:,3)*AllZRIJ)) !z
        exp_val_U_op(4)=Real(XSQ) !(x*x)2BODY
        exp_val_U_op(5)=Real(SUM(UOp_ij(:,:,4)*AllZRIJ)) !(x*x)1BODY
        exp_val_U_op(6)=(exp_val_U_op(4)+exp_val_U_op(5)
     .                           -exp_val_U_op(1)**2)/dble(Npar) ! DeltaX
        exp_val_U_op(7)=Real(YSQ) !y*y
        exp_val_U_op(8)=Real(SUM(UOp_ij(:,:,5)*AllZRIJ)) !(y*y)1BODY
        exp_val_U_op(9)=(exp_val_U_op(7)+exp_val_U_op(8)
     .                           -exp_val_U_op(2)**2)/dble(Npar) ! DeltaY
        exp_val_U_op(10)=Real(ZSQ)  !(z*z)2BODY
        exp_val_U_op(11)=Real(SUM(UOp_ij(:,:,6)*AllZRIJ))  !(z*z)1BODY
        exp_val_U_op(12)=(exp_val_U_op(10)+exp_val_U_op(11)
     .                           -exp_val_U_op(3)**2)/dble(Npar) ! DeltaZ
        exp_val_U_op(13)=Real(SUM(UOp_ij(:,:,7)*AllZRIJ)) !Trap == Potential energy
        exp_val_U_op(14)=Real(E)-exp_val_U_op(13)   !Kinetic energy
        exp_val_U_op(15)=Real(W)                   !Interaction energy
        exp_val_U_op(16)=Real(E+W)                 !Total energy

        end subroutine User_OP
