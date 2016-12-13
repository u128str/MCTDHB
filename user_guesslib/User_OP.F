 
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
       integer ::  i,j,k,Iorb,Jorb,ind,icnt,Korb,Lorb,P,ierr
C       complex*16 :: UOp_ij(Morb,Morb,10)
       complex*16 :: W,E,XSQ,YSQ,ZSQ,PXSQ,PYSQ,PZSQ
       real*8 :: exp_val_U_op(30),t,time
       COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb)  :: PSI
C       COMPLEX*16, DIMENSION(NDX*NDY*NDZ,10)  :: UOp_x_y_z 
       integer :: cK,cJ,cL,cI,L2B
       complex*16 :: rho_jk(Morb,Morb)
       complex*16 :: rho_ijkl(Morb,Morb,Morb,Morb)
C      SK 08/04/2015
C      MATRIX REPRESENTATION OF THE OPERATORS IN THE WORKING ORBITALS BASIS
       DOUBLE COMPLEX, DIMENSION(MORB,MORB) :: XMAT,YMAT,ZMAT
       DOUBLE COMPLEX, DIMENSION(MORB,MORB) :: X2MAT,Y2MAT,Z2MAT
       DOUBLE COMPLEX, DIMENSION(MORB,MORB) :: PXMAT,PYMAT,PZMAT
       DOUBLE COMPLEX, DIMENSION(MORB,MORB) :: PX2MAT,PY2MAT,PZ2MAT
       DOUBLE COMPLEX, DIMENSION(MORB,MORB) :: VMAT
C      ********************************************************************
       DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: X,Y,Z,X2,Y2,Z2
C      ********************************************************************
       COMPLEX*16, DIMENSION(NDX*NDY*NDZ)  :: DPSI
       COMPLEX*16, DIMENSION(NDX*NDY*NDZ)  :: DDPSI

       ALLOCATE(X(NDX*NDY*NDZ),stat=ierr) 
        if(ierr /= 0) write(*,*)"allocation error User_Op"
       ALLOCATE(Y(NDX*NDY*NDZ),stat=ierr) 
        if(ierr /= 0) write(*,*)"allocation error User_Op"
       ALLOCATE(Z(NDX*NDY*NDZ),stat=ierr) 
        if(ierr /= 0) write(*,*)"allocation error User_Op"
       ALLOCATE(X2(NDX*NDY*NDZ),stat=ierr) 
        if(ierr /= 0) write(*,*)"allocation error User_Op"
       ALLOCATE(Y2(NDX*NDY*NDZ),stat=ierr) 
        if(ierr /= 0) write(*,*)"allocation error User_Op"
       ALLOCATE(Z2(NDX*NDY*NDZ),stat=ierr) 
        if(ierr /= 0) write(*,*)"allocation error User_Op"
         ind=1
         Do K=1,NDZ
           Do J=1,NDY
             Do I=1,NDX
               X(ind)=ort_X(I)
               Y(ind)=ort_Y(J)
               Z(ind)=ort_Z(K)
               X2(ind)=ort_X(I)**2
               Y2(ind)=ort_Y(J)**2
               Z2(ind)=ort_Z(K)**2
               ind=ind+1
             Enddo
           Enddo
         Enddo




       t=time
       XMAT=0.d0
       YMAT=0.d0
       ZMAT=0.d0
       X2MAT=0.d0
       Y2MAT=0.d0
       Z2MAT=0.d0
       PXMAT=0.d0
       PYMAT=0.d0
       PZMAT=0.d0
       PX2MAT=0.d0
       PY2MAT=0.d0
       PZ2MAT=0.d0
       DO Iorb=1,Morb
          DO Jorb=Iorb,Morb
            XMAT(Iorb,Jorb)=
     .                SUM(PSI(:,Jorb)*(X(:))*dconjg(PSI(:,Iorb)))
            YMAT(Iorb,Jorb)=
     .                SUM(PSI(:,Jorb)*(Y(:))*dconjg(PSI(:,Iorb)))
            ZMAT(Iorb,Jorb)=
     .                SUM(PSI(:,Jorb)*(Z(:))*dconjg(PSI(:,Iorb)))
            X2MAT(Iorb,Jorb)=
     .                SUM(PSI(:,Jorb)*(X2(:))*dconjg(PSI(:,Iorb)))
            Y2MAT(Iorb,Jorb)=
     .                SUM(PSI(:,Jorb)*(Y2(:))*dconjg(PSI(:,Iorb)))
            Z2MAT(Iorb,Jorb)=
     .                SUM(PSI(:,Jorb)*(Z2(:))*dconjg(PSI(:,Iorb)))
            VMAT(Iorb,Jorb)=
     .                SUM(PSI(:,Jorb)*VTRAP_EXT(:)*dconjg(PSI(:,Iorb)))
               DPSI=PSI(:,Jorb)
               CALL Get_d_PSI(DPSI,1)
               DDPSI=DPSI
               CALL Get_d_PSI(DDPSI,1)
            PXMAT(Iorb,Jorb)=
     .                 DCMPLX(0.d0,-1.d0)*SUM(DPSI*dconjg(PSI(:,Iorb)))
            PX2MAT(Iorb,Jorb)=-1.d0*SUM(DDPSI*dconjg(PSI(:,Iorb)))
            IF (DIM_MCTDHB.GT.1) THEN
               DPSI=PSI(:,Jorb)
               CALL Get_d_PSI(DPSI,2)
               DDPSI=DPSI
               CALL Get_d_PSI(DDPSI,2)
            PYMAT(Iorb,Jorb)=
     .                 DCMPLX(0.d0,-1.d0)*SUM(DPSI*dconjg(PSI(:,Iorb)))
            PY2MAT(Iorb,Jorb)=-1.d0*SUM(DDPSI*dconjg(PSI(:,Iorb)))
            ENDIF
            IF (DIM_MCTDHB.GT.2) THEN
               DPSI=PSI(:,Jorb)
               CALL Get_d_PSI(DPSI,3)
               DDPSI=DPSI
               CALL Get_d_PSI(DDPSI,3)
            PZMAT(Iorb,Jorb)=
     .                 DCMPLX(0.d0,-1.d0)*SUM(DPSI*dconjg(PSI(:,Iorb)))
            PZ2MAT(Iorb,Jorb)=-1.d0*SUM(DDPSI*dconjg(PSI(:,Iorb)))
            ENDIF
            XMAT(Jorb,Iorb)=Conjg(XMAT(Iorb,Jorb))    
            YMAT(Jorb,Iorb)=Conjg(YMAT(Iorb,Jorb))    
            ZMAT(Jorb,Iorb)=Conjg(ZMAT(Iorb,Jorb))    
            X2MAT(Jorb,Iorb)=Conjg(X2MAT(Iorb,Jorb))    
            Y2MAT(Jorb,Iorb)=Conjg(Y2MAT(Iorb,Jorb))    
            Z2MAT(Jorb,Iorb)=Conjg(Z2MAT(Iorb,Jorb))    
            VMAT(Jorb,Iorb)=Conjg(VMAT(Iorb,Jorb))    
            PXMAT(Jorb,Iorb)=Conjg(PXMAT(Iorb,Jorb))    
            PYMAT(Jorb,Iorb)=Conjg(PYMAT(Iorb,Jorb))    
            PZMAT(Jorb,Iorb)=Conjg(PZMAT(Iorb,Jorb))    
            PX2MAT(Jorb,Iorb)=Conjg(PX2MAT(Iorb,Jorb))    
            PY2MAT(Jorb,Iorb)=Conjg(PY2MAT(Iorb,Jorb))    
            PZ2MAT(Jorb,Iorb)=Conjg(PZ2MAT(Iorb,Jorb))    
 
           EndDO   
        EndDO 

        Rho_ijkl=0.d0
        CALL Get_RhoIJKL(Rho_IJKL)

C       TWO-BODY CONTRIBUTION TO X2/Y2/Z2/PX2/PY2/PZ2
C       <A^2>_2Body=SUM_IJKL RHO(I,J,K,L)*A(I,K)*A(J,L)
        XSQ=0.d0
        YSQ=0.d0
        ZSQ=0.d0
        PXSQ=0.d0
        PYSQ=0.d0
        PZSQ=0.d0
        DO CI=1,Morb
          DO CJ=1,Morb
            DO CK=1,Morb
              DO CL=1,Morb
           XSQ=XSQ+Rho_ijkl(CI,CJ,CK,CL)*XMAT(CI,CK)*XMAT(CJ,CL)
           YSQ=YSQ+Rho_ijkl(CI,CJ,CK,CL)*YMAT(CI,CK)*YMAT(CJ,CL)
           ZSQ=ZSQ+Rho_ijkl(CI,CJ,CK,CL)*ZMAT(CI,CK)*ZMAT(CJ,CL)
           PXSQ=PXSQ+Rho_ijkl(CI,CJ,CK,CL)*PXMAT(CI,CK)*PXMAT(CJ,CL)
           PYSQ=PYSQ+Rho_ijkl(CI,CJ,CK,CL)*PYMAT(CI,CK)*PYMAT(CJ,CL)
           PZSQ=PZSQ+Rho_ijkl(CI,CJ,CK,CL)*PZMAT(CI,CK)*PZMAT(CJ,CL)
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
C       CALCULATING EXPECTATION VALUES PER PARTICLE
        exp_val_U_op=0.d0
C **************************************************************************
C                            BEGIN ENERGY PART
C **************************************************************************
        exp_val_U_op(1)=Real(E+W)/DBLE(Npar)                       !<E>/N
        exp_val_U_op(3)=Real(SUM(VMAT(:,:)*AllZRIJ))/DBLE(Npar)    !<V>/N
        exp_val_U_op(2)=Real(E)/DBLE(Npar)-exp_val_U_op(3)        !<T>/N
        exp_val_U_op(4)=Real(W)/DBLE(Npar)                         !<W>/N
C **************************************************************************
C                            BEGIN POSITION PART
C **************************************************************************
        exp_val_U_op(5)=Real(SUM(XMAT(:,:)*AllZRIJ))/DBLE(Npar) !<x>/N
        exp_val_U_op(6)=Real(SUM(YMAT(:,:)*AllZRIJ))/DBLE(Npar) !<y>/N
        exp_val_U_op(7)=Real(SUM(ZMAT(:,:)*AllZRIJ))/DBLE(Npar) !<z>/N
        exp_val_U_op(8)=                                        !<x^2>/N
     .       (Real(XSQ)+Real(SUM(X2MAT(:,:)*AllZRIJ)))/DBLE(Npar)
        exp_val_U_op(9)=                                        !<y^2>/N
     .       (Real(YSQ)+Real(SUM(Y2MAT(:,:)*AllZRIJ)))/DBLE(Npar)
        exp_val_U_op(10)=                                        !<z^2>/N
     .       (Real(ZSQ)+Real(SUM(Z2MAT(:,:)*AllZRIJ)))/DBLE(Npar)
        exp_val_U_op(11)=                                        !Var(x)/N
     .        exp_val_U_op(8)-exp_val_U_op(5)**2*DBLE(Npar) 
        exp_val_U_op(12)=                                        !Var(y)/N
     .        exp_val_U_op(9)-exp_val_U_op(6)**2*DBLE(Npar) 
        exp_val_U_op(13)=                                        !Var(z)/N
     .        exp_val_U_op(10)-exp_val_U_op(7)**2*DBLE(Npar) 
C **************************************************************************
C                            BEGIN MOMENTUM PART
C **************************************************************************
        exp_val_U_op(14)=Real(SUM(PXMAT(:,:)*AllZRIJ))/DBLE(Npar) !<px>/N
        exp_val_U_op(15)=Real(SUM(PYMAT(:,:)*AllZRIJ))/DBLE(Npar) !<py>/N
        exp_val_U_op(16)=Real(SUM(PZMAT(:,:)*AllZRIJ))/DBLE(Npar) !<pz>/N
        exp_val_U_op(17)=                                        !<px^2>/N
     .       (Real(PXSQ)+Real(SUM(PX2MAT(:,:)*AllZRIJ)))/DBLE(Npar)
        exp_val_U_op(18)=                                        !<py^2>/N
     .       (Real(PYSQ)+Real(SUM(PY2MAT(:,:)*AllZRIJ)))/DBLE(Npar)
        exp_val_U_op(19)=                                        !<pz^2>/N
     .       (Real(PZSQ)+Real(SUM(PZ2MAT(:,:)*AllZRIJ)))/DBLE(Npar)
        exp_val_U_op(20)=                                        !Var(px)/N
     .        exp_val_U_op(17)-exp_val_U_op(14)**2*DBLE(Npar) 
        exp_val_U_op(21)=                                        !Var(py)/N
     .        exp_val_U_op(18)-exp_val_U_op(15)**2*DBLE(Npar) 
        exp_val_U_op(22)=                                        !Var(pz)/N
     .        exp_val_U_op(19)-exp_val_U_op(16)**2*DBLE(Npar) 
C **************************************************************************
C                            BEGIN ****** PART
C **************************************************************************



        DEALLOCATE(X,Y,Z,X2,Y2,Z2)
        end subroutine User_OP
