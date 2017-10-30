         SUBROUTINE Get_KSL(KSL_temp,psiS)
         USE SHARED_DIMS
         USE W_INTERPARTICLE
         USE DVR_ALL
 
         IMPLICIT NONE
c           COMPLEX*16 ::  WSL(:),psiS(:),psiL(:)
         COMPLEX*16,DIMENSION(NDX*NDY*NDZ) :: psiS,psiL,WSL,temp
         COMPLEX*16,DIMENSION(NDX*NDY*NDZ) :: KSL
         COMPLEX*16,DIMENSION(NDX*NDY*NDZ,NDX*NDY*NDZ) :: KSL_temp
         COMPLEX*16 :: ZDOTU,ZDOTC,qq
         INTEGER :: I,J,K,ierr,pos,counter,x_pos,y_pos,col,i1,i2,j1,
     .                 j2,q1,q2
c==============================================================================
!         write(*,*) "Dim",DIM_MCTDHB 
!         write(*,*) "Wxx_Type",Wxx_TYPE 
         Dim_1_2_3D: SELECT CASE (DIM_MCTDHB)
c==============================================================================

         CASE (1) ! 1D  ! 1D  ! 1D  ! 1D  ! 1D  ! 1D  ! 1D  ! 1D  ! 1D
            pos=mod(col,NDX)
            if(pos.eq.0) pos=NDX

            WxxTYPE: SELECT CASE (Wxx_TYPE)
            CASE (0)
            Do i=1, NDX 
              temp=0d0
              temp(i)=1.0d0
              CALL  Get_WSL_omp(WSL,psis,temp)  ! contact case - is overheaded by filling the complit matrix KSL with zeros apart from diagonal
              KSL_temp(:,i)=WSL(:)
            End Do

            CASE (1) !W(x-x')*W(y-y')*W(z-z')
            Do J=1,NDX
             KSL(J)=xlambda0*Conjg(psiS(pos))*PsiL(J)
     -     *WOp_X(pos,J)
             Enddo

c====================================
            CASE (2) !W(r,r')=W(r-r') Case UnSeparable interparticel potential W(r-r') THIS IS TRUE ONLY FOR EQUIDISTANT DVR-GRID!!!"
            STOP "KSL case 2 not implemented yet"

            CASE (3) !======================== General Case UnSeparable interparticel potential W(r,r')
            Do J=1,NDX
             KSL(J)=xlambda0*Conjg(psiS(pos))*PsiL(J)
     .     *Wxx(pos,J)
            Enddo

c=====================
            CASE (4) ! !W(r,r')=W(r-r') FOR FFT_WSL  TRUE ONLY FOR EQUIDISTANT DVR-GRID!!!"
               Do I=1,NDX
                           temp=0d0
                           temp(I)=1.0d0
           call computewsl(wsl,psis,temp,w3xxFFT)
c               KSL(I,J)=WSL(I)*psiL(i)
c               KSL(I)=WSL(pos)*psiL(pos)
                KSL_temp(:,I)=WSL(:)
               ENDDO
c=====================
            END SELECT WxxTYPE
         RETURN

         CASE (2) ! 2D  ! 2D  ! 2D  ! 2D  ! 2D  ! 2D  ! 2D  ! 2D  ! 2D
            pos=mod(col,NDX*NDY)
            if(pos.eq.0) pos=NDX*NDY


            WxxTYPE2: SELECT CASE (Wxx_TYPE)
            CASE (0) 
            WRITE(6,*)" contact 2D not yet implemented"
            stop
            CASE (1) !W(x-x')*W(y-y')*W(z-z')
            counter=1

             if(mod(col,NDX).eq.0) then
               y_pos=col/NDX
             else
               y_pos=(col/NDX)+1
             end if
 
             x_pos=mod(col,NDX)
             if(x_pos.eq.0) x_pos=NDX


!            Do I=1,NDX
!             do j=1, NDY

             DO j1=1,NDY
              j2=y_pos
              if(j1.eq.j2) then
               q1=1
              else
               q1=0
              endif
              Do i1=1,NDX
               i2=x_pos
               if(i1.eq.i2) then
                q2=1
               else
                q2=0
               endif
               I=NDX*(j1-1)+i1
!               Tkin_col(I)=q1*Op_X(i1,i2)+q2*Op_Y(j1,j2)
                 KSL(I)=xlambda0*Conjg(psiS(I))*PsiL(pos)
     -     *(WOp_X(i1,i2)+WOp_y(j1,j2))
                 counter=counter+1
              endDo
             ENDDO


!               KSL(counter)=xlambda0*Conjg(psiS(counter))*PsiL(pos)
!     -     *(WOp_X(pos,I)+WOp_y(pos,j))
!               counter=counter+1
!             end do 
!            Enddo
c====================================
            CASE (2) !W(r,r')=W(r-r') Case UnSeparable interparticel potential W(r-r') THIS IS TRUE ONLY FOR EQUIDISTANT DVR-GRID!!!"
            STOP "KSL case 2 not implemented yet"

            CASE (3) !======================== General Case UnSeparable interparticel potential W(r,r')
            STOP "KSL case 3 not implemented" 
            counter=1
            Do I=1,NDX
             do j=1, NDY
            KSL(counter)=xlambda0*Conjg(psiS(counter))*PsiL(pos)
     .     *Wxx(pos,counter)
               counter=counter+1
             end do 
            Enddo
c=====================
            CASE (4) ! !W(r,r')=W(r-r') FOR FFT_WSL  TRUE ONLY FOR EQUIDISTANT DVR-GRID!!!"
              counter=1
               Do I=1,NDX
                do j=1, NDY
                           temp=0d0
                           temp(counter)=1.0d0
!                           temp(pos)=1.0d0
           call computewsl(wsl,psis,temp,w3xxFFT)
c               KSL(I,J)=WSL(I)*psiL(i)
                 KSL_temp(:,counter)=WSL(:)
!                 KSL(counter)=WSL(pos)*psiL(pos)
                 counter=counter+1
                end do 
               ENDDO
c=====================
            END SELECT WxxTYPE2
         RETURN



         CASE (3) ! 3D  ! 3D  ! 3D  ! 3D  ! 3D  ! 3D  ! 3D  ! 3D  ! 3D
            WxxTYPE3: SELECT CASE (Wxx_TYPE)
            CASE (0) 
            CALL  Get_WSL_omp(WSL,psis,psil)  ! contact case - is overheaded by filling the complit matrix KSL with zeros apart from diagonal
            KSL=Zero
            KSL(pos)=WSL(pos)

            CASE (1) !W(x-x')*W(y-y')*W(z-z')
            counter=1
            Do I=1,NDX
             do j=1, NDY
              do k=1, NDZ
               KSL(counter)=xlambda0*Conjg(psiS(counter))*PsiL(pos)
     -     *(WOp_X(pos,i)+WOp_y(pos,j)+WOp_z(pos,k))
               counter=counter+1
              end do 
             end do 
            Enddo
c====================================
            CASE (2) !W(r,r')=W(r-r') Case UnSeparable interparticel potential W(r-r') THIS IS TRUE ONLY FOR EQUIDISTANT DVR-GRID!!!"
            STOP "KSL case 2 not implemented yet"

            CASE (3) !======================== General Case UnSeparable interparticel potential W(r,r')
            STOP "KSL case 3 not implemented yet"
            counter=1
            Do I=1,NDX
             do j=1, NDY
              do k=1, NDZ
            KSL(counter)=xlambda0*Conjg(psiS(counter))*PsiL(pos)
     .     *Wxx(pos,counter)
               counter=counter+1
              end do
             end do 
            Enddo
c=====================
            CASE (4) ! !W(r,r')=W(r-r') FOR FFT_WSL  TRUE ONLY FOR EQUIDISTANT DVR-GRID!!!"
              counter=1
               Do I=1,NDX
                do j=1, NDY
                  do k=1, NDZ
                           temp=0d0
                           temp(counter)=1.0d0
           call computewsl(wsl,psis,temp,w3xxFFT)
                   KSL(counter)=WSL(pos)*psiL(pos)
                   counter=counter+1
                  end do                  
                end do 
               ENDDO
c=====================
            END SELECT WxxTYPE3
         RETURN


c===============================================================================
         CASE DEFAULT
         WRITE(6,*)" In KSL case Default???? Why?,Wxx_TYPE",Wxx_TYPE
         WRITE(6,*)" STOP ERRRR "
         stop 
         END SELECT Dim_1_2_3D
         END  SUBROUTINE Get_KSL
