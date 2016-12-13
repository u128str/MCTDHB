         SUBROUTINE Get_KSL(KSL,psiS,psiL)
         USE SHARED_DIMS
         USE W_INTERPARTICLE
         USE DVR_ALL
 
         IMPLICIT NONE
c           COMPLEX*16 ::  WSL(:),psiS(:),psiL(:)
         COMPLEX*16,DIMENSION(NDX*NDY*NDZ) :: psiS,psiL,WSL,temp
         COMPLEX*16,DIMENSION(NDX*NDY*NDZ,NDX*NDY*NDZ) :: KSL
         COMPLEX*16 :: ZDOTU,ZDOTC,qq
         INTEGER :: I,J,K,ierr
c==============================================================================
         Dim_1_2_3D: SELECT CASE (DIM_MCTDHB)
c==============================================================================
         CASE (1) ! 1D  ! 1D  ! 1D  ! 1D  ! 1D  ! 1D  ! 1D  ! 1D  ! 1D
            WxxTYPE: SELECT CASE (Wxx_TYPE)
            CASE (0) 
            CALL  Get_WSL_omp(WSL,psis,psil)  ! contact case - is overheaded by filling the complit matrix KSL with zeros apart from diagonal
            KSL=Zero
            Do J=1,NDX
            KSL(J,J)=WSL(J)
            Enddo
            CASE (1) !W(x-x')*W(y-y')*W(z-z')
            Do I=1,NDX
            Do J=1,NDX
            KSL(J,I)=xlambda0*Conjg(psiS(I))*PsiL(J)
     -     *WOp_X(I,J)
            Enddo
            Enddo
c=====================
!            temp=0d0
!            Do I=1,NDX
!            Do J=1,NDX
!            temp(i)=temp(i)+KSL(I,J)*PsiL(J)
!            Enddo
!            Enddo
!       OPEN(unit=117,File='testW1.out')
!       DO i=1,NDX
!             WRITE(117,'(I4,2E16.8)') i,  temp(i)
!       END DO
!       CLOSE(117)
!            stop "ok W1"
c====================================
            CASE (2) !W(r,r')=W(r-r') Case UnSeparable interparticel potential W(r-r') THIS IS TRUE ONLY FOR EQUIDISTANT DVR-GRID!!!"
            STOP "KSL case 2 not implemented yet"

            CASE (3) !======================== General Case UnSeparable interparticel potential W(r,r')
            Do I=1,NDX
            Do J=1,NDX
            KSL(J,I)=xlambda0*Conjg(psiS(I))*PsiL(J)
     .     *Wxx(I,J)
            Enddo
            Enddo
c=====================
!            temp=0d0
!            Do I=1,NDX
!            Do J=1,NDX
!            temp(i)=temp(i)+KSL(I,J)*PsiL(J)
!            Enddo
!            Enddo
!       OPEN(unit=117,File='testW3.out')
!       DO i=1,NDX
!             WRITE(117,'(I4,2E16.8)') i,  temp(i)
!       END DO
!       CLOSE(117)
!            stop "ok W3"
c=====================
            CASE (4) ! !W(r,r')=W(r-r') FOR FFT_WSL  TRUE ONLY FOR EQUIDISTANT DVR-GRID!!!"
              Do J=1,NDX
                           temp=0d0
                           temp(J)=1.0d0
           call computewsl(wsl,psis,temp,w3xxFFT)
               Do I=1,NDX
               KSL(I,J)=WSL(I)*psiL(i)
               ENDDO
               ENDDO
!            STOP "KSL case 4 is not implemented yet"
c=====================
!            temp=0d0
!            Do I=1,NDX
!            Do J=1,NDX
!            temp(i)=temp(i)+KSL(I,J)*PsiL(J)
!            Enddo
!            Enddo
!       OPEN(unit=117,File='testW4.out')
!       DO i=1,NDX
!             WRITE(117,'(I4,2E16.8)') i,  temp(i)
!       END DO
!       CLOSE(117)
!            stop "ok W4"
c=====================
            END SELECT WxxTYPE



         RETURN
         CASE (2) ! 2D  ! 2D  ! 2D  ! 2D  ! 2D  ! 2D  ! 2D  ! 2D  ! 2D
         WRITE(6,*)" STOP KSL 2D is still not  implemented"
         stop
         CASE (3) ! 3D  ! 3D  ! 3D  ! 3D  ! 3D  ! 3D  ! 3D  ! 3D  ! 3D
         WRITE(6,*)" STOP KSL 3D is still not  implemented"
         stop
c===============================================================================
         CASE DEFAULT
         WRITE(6,*)" In KSL case Default???? Why?,Wxx_TYPE",Wxx_TYPE
         WRITE(6,*)" STOP ERRRR "
         stop 
         END SELECT Dim_1_2_3D
         END  SUBROUTINE Get_KSL
