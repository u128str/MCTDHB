C=============== Explicite construction of the Kinetic energy OPERATOR needed for direct diagonalization 
       subroutine T_FFT_1D(h2)
       USE SHARED_DIMS
       USE DVR_ALL
       USE Parallel_Orb
       USE rR_hW
       implicit NONE
       REAL*8 :: local_time
c===================================================================
       INTEGER ::  i,j,ierr
       COMPLEX*16,  DIMENSION(NDX,NDX) :: h2 
       COMPLEX*16, ALLOCATABLE :: tempX(:), tempY(:)
c===================================================================
         allocate(tempX(NDX),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in T_FFT_1D"
         allocate(tempY(NDX),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in T_FFT_1D"
            write(6,*)"T FFT ",DIM_MCTDHB,Time_DVRMETHODX, NDX
c=============   K  I  N  E  T  I  C  S ===============================
        IF(DIM_MCTDHB.ne.1) stop "T_FFT_1D only in 1D"
        IF(Time_DVRMETHODX==4) THEN
       Do i=1,NDX
        tempX=(0d0,0d0)
        tempX(i)=(1.0d0,0d0)
#ifdef SFX1D              
                   CALL Get_Op_PSI_ SFX1D (tempX,1)  !FFT 1D in X direction
#endif
        Do j=1,NDX
        tempY=(0d0,0d0)
        tempY(j)=(1.0d0,0d0)
        h2(i,j)= sum(conjg(tempY(:))*tempX(:))
        End Do
       EndDO
        ENDIF
c        deallocate(tempX,tempY,tempZ)
        deallocate(tempX,tempY)
c       write(6,*)"Orbital  Assistent ",MYID," is DONE"
       return
       end subroutine T_FFT_1D
