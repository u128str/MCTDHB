! Subroutine to provide WIJ -> -part of two-body orbital part
! NtVec - number of the natural vectors in the expansion of W(r-r') over
! natural potentials  and respective occupation numbers F  Dim(F)= NtVec   
! Number of Natural potentials FIJ is of course equal to the number of
! natural potential occupation numbers DIM(FIJ) = NDX*NDY*NDZ : NtVec 
!        I checked both possibilites  with and without natpot!!!!!!!!!!!

          SUBROUTINE Get_WSL(WSL,psiS,psiL,WSLTime)
          USE SHARED_DIMS
          USE W_INTERPARTICLE
          USE DVR_ALL
 
          IMPLICIT NONE
          REAL*8 :: WSLTime
c           COMPLEX*16 ::  WSL(:),psiS(:),psiL(:)
         COMPLEX*16,DIMENSION(NDX*NDY*NDZ) :: psiS,psiL
         COMPLEX*16,DIMENSION(NDX*NDY*NDZ) :: WSL
c         COMPLEX*16,DIMENSION(NDX*NDY*NDZ) ::temp
         COMPLEX*16, allocatable :: temp(:)

          COMPLEX*16 :: ZDOTU,ZDOTC,qq
c==============================================================================
          INTEGER :: I,J,K,ierr
c==============================================================================
c============================== If (.FALSE.) W(x-x')=delta(x-x')*xlambda0 is used
c============================== or if (.TRUE.) - true interparticle potential 
c============================== W(r-r')*xlambda0 is used
c============================== defined by Wxx(I,J)=W(r,r') in DVR grid
c============================== or by vector W2xx(J)=W(ort(1)-ort(J))
c============================== or for separable case: WOP_X*WOP_Y*WOP_Z
c=========================== Delta case ========================================
c         print *, Wxx_TYPE,"WSL0", SUM(WSL*WSL)
         IF(WTD.eqv..TRUE.) CALL Get_InterPart(WSLTime)  !AIS 15JUl2014 For time-dependent interactions
          
         WxxTYPE: SELECT CASE (Wxx_TYPE)
c=========================== Delta case ========================================
         CASE (0)
         WSL=xlambda0*Conjg(psiS)*PsiL/weight/weight
         RETURN
c============================General Cases ======================================
c=================== Case Separable W WOp_X
         CASE (1)
                            WSL=Conjg(psiS)*psiL*xlambda0
                            CALL Get_Op_PSI(WSL,WOp_X,NDX,1)
         IF(DIM_MCTDHB.eq.2) CALL Get_Op_PSI(WSL,WOp_Y,NDY,2)
         IF(DIM_MCTDHB.eq.3) CALL Get_Op_PSI(WSL,WOp_Z,NDZ,3)
c=================== Case Unseparable W(x-x') - only one coloumn is needed
         CASE (2)
        allocate(temp(NDX*NDY*NDZ),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in Get_WSL"
         temp=Conjg(psiS)*psiL*xlambda0
         Do I=1,NDX*NDY*NDZ
         WSL(I)=ZDOTU(NDX*NDY*NDZ,W2xx(NDX*NDY*NDZ+1-I),1,temp,1)
         Enddo
         deallocate(temp)
         RETURN
c=================== Case Unseparable W(x,x') - all elements Wxx are needed 
         CASE (3)
        allocate(temp(NDX*NDY*NDZ),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in Get_WSL"
         temp=Conjg(psiS)*psiL*xlambda0
       CALL ZSYMV('L',NDX*NDY,ZONER,Wxx,NDX*NDY,temp,1,ZERO,WSL,1)
         deallocate(temp)
         RETURN
         case (4)
           call computewsl(wsl,psis,psil,w3xxFFT)
         RETURN
c===============================================================================
         CASE DEFAULT
         WRITE(6,*)" In WSL case Default???? Why?,Wxx_TYPE",Wxx_TYPE
         WRITE(6,*)" STOP ERRRR "
         stop 
         END SELECT WxxTYPE
         END  SUBROUTINE Get_WSL
