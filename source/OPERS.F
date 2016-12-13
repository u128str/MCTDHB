C=============== OPERATORS are acting here  
       subroutine OPERS(temp,local_time)
       USE SHARED_DIMS
       USE DVR_ALL
       USE Parallel_Orb
       USE rR_hW
       implicit NONE
c=================== MPI 
       INCLUDE 'mpif.h'
       INTEGER ::  ierr,MYID,numprocs
c===================================================================
         REAL*8 :: local_time
c===================================================================
       INTEGER ::  n,Iorb,IPRC,m,nn,icntr
       COMPLEX*16,  DIMENSION(NDX*NDY*NDZ,Morb) :: temp ! STR 2013,PSI
       COMPLEX*16, ALLOCATABLE :: tempX(:), tempIN(:)
       logical :: threeD
       EXTERNAL Get_Op_PSI,ZSCAL,VTRAP_EXT_TD
c===================================================================
         allocate(tempX(NDX*NDY*NDZ),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in OPERS Asistent"
!
         allocate(tempIN(NDX*NDY*NDZ),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in OPERS Asistent"
c========================= MPI =====================================


         call MPI_COMM_RANK(MPI_COMM_WORLD,MYID,ierr)
         call MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,ierr)
c===================================================================
        IPRC=MYID+1
c        CALL MPI_ScatterV(PSI,Orb_Block,Orb_Displ,MPI_DOUBLE_COMPLEX,
c     .  temp,Orb_Block(IPRC),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)

!            temp has all orbitals
c===================== GET EXTERNAL TRAP ===========================
          CALL VTRAP_EXT_TD(local_time)
c           write(6,*)"MID",IPRC
c===================================================================
                 icntr=1
c       write(6,*)"from",Proc_Iorb_Starts(IPRC),"til",
c     . Proc_Iorb_Finish(IPRC)

      I_ORBITALS: DO Iorb=Proc_Iorb_Starts(IPRC),Proc_Iorb_Finish(IPRC)
c                write(6,*)"M/S",IPRC,"IORB",iorb
c       tempIN=temp(:,icntr)
       tempIN=temp(:,Iorb)
       threed=.TRUE.
       FFT_2D=.TRUE.
!       FFT_2D=.False.
       threed=.True.
!       threed=.False.
c=============   K  I  N  E  T  I  C  S ===============================
        tempX=tempIN
      IF((Time_DVRMETHODX==4).AND.(Time_DVRMETHODY==4).and.
     $ (threed.eqv..true.)
     $      .AND.(DIM_MCTDHB==3).AND.(Time_DVRMETHODZ==4)) THEN
#ifdef SFX3D              
         CALL Get_Op_PSI_ SFX3D (tempX,5)  !TRUE 3D FFT 3D MKL FFTW or CUDA CPP or 2D plus 1D 
#endif
          temp(:,icntr)=tempX
        Else
        IF((Time_DVRMETHODX==4).AND.(Time_DVRMETHODY==4)
     $   .and.
     $   (DIM_MCTDHB==2).AND.(FFT_2D.eqv..TRUE.)) THEN
#ifdef SFX2D              
                   CALL Get_Op_PSI_ SFX2D (tempX,4)  !FFT 2D MKL or FFTW
!             print * ,"Kernel 2D in" 
#endif
        temp(:,icntr)=tempX
        ELSE
        IF(Time_DVRMETHODX==4) THEN
#ifdef SFX1D              
                   CALL Get_Op_PSI_ SFX1D (tempX,1)  !FFT 1D in X direction
#endif
        ELSE
        CALL Get_Op_PSI(tempX,Op_X,NDX,1)
        ENDIF
        temp(:,icntr)=tempX
        IF(DIM_MCTDHB.ge.2) then
        tempX=tempIN
        IF(Time_DVRMETHODY==4) THEN
#ifdef SFX1D              
                   CALL Get_Op_PSI_ SFX1D(tempX,2)  !FFT 1D in Y direction
#endif
        ELSE
        CALL Get_Op_PSI(tempX,Op_Y,NDY,2)
        ENDIF
        temp(:,icntr)=temp(:,icntr)+tempX
        endif
        IF(DIM_MCTDHB.ge.3) then
       tempX=tempIN
      IF(Time_DVRMETHODZ==4) THEN
!        write(6,*)"STR FFT-Z direction"
#ifdef SFX1D              
                   CALL Get_Op_PSI_ SFX1D (tempX,3)  !FFT 1D in Z direction
#endif
!        write(6,*)"STR FFT-Z direction done"
      ELSE
      CALL Get_Op_PSI(tempX,Op_Z,NDZ,3)
      ENDIF
        temp(:,icntr)=temp(:,icntr)+tempX
        endif
        endif
        ENDIF
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCC ANGULAR MOMENTUM CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

!       IF (LZ.eq..TRUE.) then
!        tempX=tempIN
!        CALL Get_Op_PSI_L_FFT_MKL(tempX,3)  !FFT
cccccccccccccccccccccccccccccccccccc OMEGA_Z ccccccccccccccccccccc
!        temp(:,icntr) = temp(:,icntr)-OMEGAZ*tempX
cccccccccccccccccccccccccccccccccccc OMEGA_Z ccccccccccccccccccccc
!       endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCC ANGULAR MOMENTUM CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c=============    T  R  A  P     P O T E N T I A L ===================
c        temp(:,icntr)=PSI1(:,Iorb) + VTRAP_EXT*PSI(:,Iorb) ! KInetic + EXTERNAL TRAP
c        temp(:,icntr)=temp(:,icntr) + VTRAP_EXT*PSI(:,Iorb) ! KInetic + EXTERNAL TRAP
        temp(:,icntr)=temp(:,icntr) + VTRAP_EXT*tempIN ! KInetic + EXTERNAL TRAP
                 icntr=icntr+1

                      END DO I_ORBITALS

!                temp  has the results 
c================ Gathering togeather MPI_REDUCE ===================
c       CALL  MPI_GatherV(temp,Orb_Block(IPRC),MPI_DOUBLE_COMPLEX,PSI1,
c     . Orb_Block,Orb_Displ,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)

c        deallocate(tempX,tempY,tempZ)
!       write(6,*)"Orbital  Assistent ",sum(tempX),sum(tempIN)
!       write(6,*)"Orbital  Assistent ",size(tempX),size(tempIN)
        IF(ALLOCATED(tempX)) deallocate(tempX)
!          write(6,*) "OK"
        IF(ALLOCATED(tempIN)) deallocate(tempIN) 
!       write(6,*)"Orbital  Assistent ",MYID," is DONE"
       return
       end subroutine OPERS
