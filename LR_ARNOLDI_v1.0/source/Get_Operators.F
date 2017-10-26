! Subroutine to provide Op_x, Op_Y,Op,_Z for module DVR_ALL
         SUBROUTINE Get_Operators
         USE SHARED_DIMS
         USE DVR_ALL
         USE W_INTERPARTICLE
         IMPLICIT NONE
c         INCLUDE 'mpif.h'
c         INTEGER ::  ierr,MYID,numprocs
c         REAL*8, DIMENSION(NDX,NDY,NDZ) :: weight3D
c         REAL*8, allocatable :: weight3D(:,:,:)
c         REAL*8, DIMENSION(NDX*NDX) :: XX
c         REAL*8, DIMENSION(NDY*NDY) :: YY
c         REAL*8, DIMENSION(NDZ*NDZ) :: ZZ

c==============================================================================
          INTEGER :: I,J,K,DVRMETHOD,ierr,ind
          REAL*8 :: xi,xf
          EXTERNAL DCOPY
c==============================================================================
c            allocate(weight3D(NDX,NDY,NDZ),stat=ierr)
c         if(ierr /= 0)write(*,*)"allocation error in Get_Operators"

          Op_X=0d0
          Op_Y=0d0
          Op_Z=0d0
c=========DVRMETHOD:      1         3       5     4
c=========Corrsponds:  HO DVR    SinDVR  ExpDVR  FFT
c==============================================================================
          DVRMETHOD=Time_DVRMETHODX
c          DVRMETHOD=5
          xi=Time_xint
          xf=Time_xfnl
          CALL Init_h(xi,xf,NDX,Op_X,weight_X,ort_X,Ort_Kx,
     .  dif1matX,DVRMETHOD) 
          ort_Y=0.0d0
          ort_Z=0.0d0
          weight_Y=1.0d0
          weight_Z=1.0d0
          DVRMETHOD=Time_DVRMETHODY
          xi=Time_yint
          xf=Time_yfnl
       IF(NDY.NE.1) 
     & CALL Init_h(xi,xf,NDY,Op_Y,weight_Y,ort_Y,Ort_Ky,
     .  dif1matY,DVRMETHOD)
          DVRMETHOD=Time_DVRMETHODZ
          xi=Time_zint
          xf=Time_zfnl
       IF(NDZ.NE.1) 
     . CALL Init_h(xi,xf,NDZ,Op_Z,weight_Z,ort_Z,Ort_Kz,
     . dif1matZ,DVRMETHOD)
          write(6,*)"        Problem is ",DIM_MCTDHB," Dimensional"
          write(6,*)" total dimension of the DVR grids is ",NDX*NDY*NDZ
c=================== Formation of the 3D  -weight arrays
            ind=1
           Do K=1,NDZ
           Do J=1,NDY
           Do I=1,NDX
c          weight3D(I,J,K)=weight_X(I)*weight_Y(J)*weight_Z(K)
           weight(ind)=weight_X(I)*weight_Y(J)*weight_Z(K)
            ind=ind+1
           EndDo
           EndDo
           EndDo
!          Write(6,*)" Weight _X ",weight_X(1)
!          Write(6,*)" Weight _Y ",weight_Y(1)
!          Write(6,*)" Weight _Z ",weight_Z(1)
!          Write(6,*)" Weight(1) ",weight(1)
!           pause
c          Write(6,*)" Weight(NDX*NDY*NDZ) ",weight(NDX*NDY*NDZ)
c          Write(6,*)" Ort_X ",Ort_X(1)
c          Write(6,*)" Ort_Y ",Ort_Y(1)
c          Write(6,*)" Ort_Z ",Ort_Z(1)
c          Write(6,*)" Weight  ",weight
cc          Write(6,*) SIZE(weight3D),SIZE(weight)
c          pause
c          CALL DCOPY(NDX*NDY*NDZ,weight3D,1,weight,1)
c============================== MPI PART ===============================          
c           write(6,*) NDX*NDY*NDZ
c           pause
c========================= MPI =====================================
C         call MPI_COMM_RANK(MPI_COMM_WORLD,MYID,ierr)
C         call MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,ierr)
C          Write(6,*)MYID," MPI",NUMPROCS
Cc===================== GET Op_X Op_Y Op_Z  =========================
C      CALL DCOPY(NDX*NDX,Op_X,1,XX,1)
C      CALL MPI_BCAST(XX,NDX*NDX,MPI_DOUBLE_PRECISION,0,
C     . MPI_COMM_WORLD,ierr)
C      write(6,*)"MASTER ",MYID," SEND X OP matrix",ierr
C      CALL DCOPY(NDY*NDY,Op_Y,1,YY,1)
C      CALL MPI_BCAST(YY,NDY*NDY,MPI_DOUBLE_PRECISION,0,
C     . MPI_COMM_WORLD,ierr)
C      write(6,*)"MASTER ",MYID," SEND Y OP matrix",ierr
C      CALL DCOPY(NDZ*NDZ,Op_Z,1,ZZ,1)
C      CALL MPI_BCAST(ZZ,NDZ*NDZ,MPI_DOUBLE_PRECISION,0,
C     . MPI_COMM_WORLD,ierr)
C      write(6,*)"MASTER ",MYID," SEND ALL OP matrices"
Cc      CALL MPI_Barrier(MPI_COMM_WORLD)
c          deallocate(weight3D)
         END  SUBROUTINE Get_Operators
