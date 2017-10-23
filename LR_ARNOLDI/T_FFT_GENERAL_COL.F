C=============== Explicite construction of the Kinetic energy OPERATOR COLUMN!! needed for direct diagonalization 
       subroutine T_FFT_GENERAL(h2_row,row)
       USE SHARED_DIMS
       USE DVR_ALL
       USE Parallel_Orb
       USE rR_hW
       implicit NONE
       REAL*8 :: local_time
c===================================================================
       INTEGER ::  i,j,k,ierr,counter,row
       COMPLEX*16,  DIMENSION(NDX*NDY*NDZ) :: h2_row 
       COMPLEX*16, ALLOCATABLE :: temp1(:), temp2(:)
c===================================================================
         allocate(temp1(NDX*NDY*NDZ),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in T_FFT_GENERAL"
         allocate(temp2(NDX*NDY*NDZ),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in T_FFT_GENERAL"
!            write(6,*)"T FFT ",DIM_MCTDHB, NDX*NDY*NDZ
c=============   K  I  N  E  T  I  C  S ===============================
           temp2=(0d0,0d0)
!           temp2(row)=(1.0d0,0d0)
          
 
!       counter=1
!       Do i=1,NDX
!        do j=1, NDY
!          do k=1,NDZ
 
           temp1=(0d0,0d0)
!           temp1(counter)=(1.0d0,0d0)
           temp1(row)=(1.0d0,0d0)
           temp2=temp1

         if(dim_mctdhb.eq.1) then
#ifdef SFX1D              
                   CALL Get_Op_PSI_ SFX1D (temp1,1)  !FFT 1D in X direction
#endif
         end if


         if(dim_mctdhb.eq.2) then
#ifdef SFX2D              
                   CALL Get_Op_PSI_ SFX2D (temp1,4)  !FFT 2D 
#endif
         end if
         if(dim_mctdhb.eq.3) then
#ifdef SFX3D              
                   CALL Get_Op_PSI_ SFX3D (temp1,5)  !FFT 2D 
#endif
         end if


!       IF (LZ.eqv..TRUE.) then
!        CALL Get_Op_PSI_L_FFT_MKL(temp2,3)  !FFT
cccccccccccccccccccccccccccccccccccc OMEGA_Z ccccccccccccccccccccc
!        temp1 = temp1-OMEGAZ*temp2
cccccccccccccccccccccccccccccccccccc OMEGA_Z ccccccccccccccccccccc
!       endif

!          h2_row(counter)= dot_product(temp2(:),temp1(:))
          h2_row= temp1
!          h2_row(counter)= temp1(row)
!          counter=counter+1 
!         end do
!        end do
!       EndDO


c        deallocate(tempX,tempY,tempZ)
        deallocate(temp1,temp2)
c       write(6,*)"Orbital  Assistent ",MYID," is DONE"
       return
       end subroutine T_FFT_GENERAL
