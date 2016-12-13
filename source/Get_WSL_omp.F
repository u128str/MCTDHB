! Subroutine to provide WIJ -> -part of two-body orbital part
! NtVec - number of the natural vectors in the expansion of W(r-r') over
! natural potentials  and respective occupation numbers F  Dim(F)= NtVec   
! Number of Natural potentials FIJ is of course equal to the number of
! natural potential occupation numbers DIM(FIJ) = NDX*NDY*NDZ : NtVec 
!        I checked both possibilites  with and without natpot!!!!!!!!!!!

          SUBROUTINE Get_WSL_omp(WSL,psiS,psiL,WSLTime)
          USE SHARED_DIMS
          USE W_INTERPARTICLE
          USE DVR_ALL
          USE omp_lib
 
          IMPLICIT NONE
         REAL*8 :: WSLTime
         COMPLEX*16,DIMENSION(NDX*NDY*NDZ) :: psiS,psiL
         COMPLEX*16,DIMENSION(NDX*NDY*NDZ) :: WSL
c         COMPLEX*16,DIMENSION(NDX*NDY*NDZ) ::temp
         COMPLEX*16, allocatable :: temp(:)

          COMPLEX*16 :: ZDOTU,ZDOTC,qq
c==============================================================================
          INTEGER :: I,J,K,ierr,ii
c==============================================================================
c============================== If (.FALSE.) W(x-x')=delta(x-x')*xlambda0 is used
c============================== or if (.TRUE.) - true interparticle potential 
c============================== W(r-r')*xlambda0 is used
c============================== defined by Wxx(I,J)=W(r,r') in DVR grid
c============================== or by vector W2xx(J)=W(ort(1)-ort(J))
c============================== or for separable case: WOP_X*WOP_Y*WOP_Z
c=========================== Delta case ========================================
!         print *, Wxx_TYPE,"WSL0", WTD,WSLTime
         IF(WTD.eqv..TRUE.) CALL Get_InterPart(WSLTime)  !AIS 15JUl2014 For time-dependent interactions
          
         WxxTYPE: SELECT CASE (Wxx_TYPE)
c=========================== Delta case ========================================
         CASE (0)
!$OMP PARALLEL DO
!$OMP& PRIVATE(ii)  
!$OMP& SHARED (WSL,PSIL,PSIS,xlambda0,weight,ndx,ndy,ndz)
                Do ii=1, ndx*ndy*ndz
        WSL(ii)=xlambda0*Conjg(psiS(ii))*PsiL(ii)/weight(ii)/weight(ii)
                Enddo
!$END OMP PARALLEL DO
c          WSL=xlambda0*Conjg(psiS)*PsiL/weight/weight
         RETURN
c============================General Cases ======================================
c=================== Case Separable W WOp_X
         CASE (1)
!$OMP PARALLEL DO
!$OMP& PRIVATE(ii)  
!$OMP& SHARED (WSL,PSIL,PSIS,xlambda0,ndx,ndy,ndz)
                Do ii=1, ndx*ndy*ndz
        WSL(ii)=xlambda0*Conjg(psiS(ii))*PsiL(ii)
                Enddo
!$END OMP PARALLEL DO
c                            WSL=Conjg(psiS)*psiL*xlambda0
                            CALL Get_Op_PSI(WSL,WOp_X,NDX,1)
         IF(DIM_MCTDHB.ge.2) CALL Get_Op_PSI(WSL,WOp_Y,NDY,2)
         IF(DIM_MCTDHB.ge.3) CALL Get_Op_PSI(WSL,WOp_Z,NDZ,3)
c=================== Case Unseparable W(x-x') - only one coloumn is needed
         CASE (2)
        allocate(temp(NDX*NDY*NDZ),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in Get_WSL"
c         temp=Conjg(psiS)*psiL*xlambda0
!$OMP PARALLEL DO
!$OMP& PRIVATE(ii)  
!$OMP& SHARED (temp,PSIL,PSIS,xlambda0,ndx,ndy,ndz)
                Do ii=1, ndx*ndy*ndz
        temp(ii)=xlambda0*Conjg(psiS(ii))*PsiL(ii)
                Enddo
!$END OMP PARALLEL DO

!$OMP PARALLEL DO
!$OMP& PRIVATE(I)  
!$OMP& SHARED (W2xx,temp,ndx,ndy,ndz)
         Do I=1,NDX*NDY*NDZ
         WSL(I)=ZDOTU(NDX*NDY*NDZ,W2xx(NDX*NDY*NDZ+1-I),1,temp,1)
         Enddo
!$END OMP PARALLEL DO

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
         END  SUBROUTINE Get_WSL_OMP
