! PROGRAM TO TEST THE OPTIONAL ARGUMENTS TO RESHAPE
! Subroutine to provide Op Psi -> psi Op - is operator Matrix in the
! respective DVR presentation DIM(Op)=(NDX,NDX)
! X_i index pointing to the T_x,T_y or T_z will be calculated - only for
! 2D and 3D problems
! NdimMax = 1D(NDX)
! NdimMax = 2D(NDX*NDY)
! NdimMax = 3D(NDX*NDY*NDZ)
          SUBROUTINE Get_Op_PSI_Gnrl(psi,Op,NDop,X_i) 
          USE SHARED_DIMS
          USE DVR_ALL
          IMPLICIT NONE
          INTEGER                     :: X_i,Ndop
          COMPLEX*16,   DIMENSION(NDX*NDY*NDZ) :: psi
          REAL*8,       DIMENSION(NDop,NDop) :: Op
c          COMPLEX*16,   DIMENSION(NDop,NDop) :: COp
c          COMPLEX*16   Op(*)
          
          COMPLEX*16,   DIMENSION(NDX) :: ZX
          COMPLEX*16,   DIMENSION(NDY) :: ZY
          COMPLEX*16,   DIMENSION(NDZ) :: ZZ
c==================== Local Arrays for 3D Transformations 
          INTEGER :: I,J,K,I1,incr,j1
          EXTERNAL ZSYMV,ZCOPY
c           Cop=ZONER*Op
c===============================================================================
            Xi1: SELECT CASE (X_i)
         CASE (1)
          ZX=ZERO
c           write(6,*) "X dim ==",NDX
            Do i=1,NDY*NDZ
            i1=1+(i-1)*NDX
c       CALL ZSYMV('L',NDX,ZONER,Cop,NDop,PSI(i1),1,ZERO,PSI(i1),1)
c       CALL ZSYMV('L',NDX,ZONER,Cop,NDop,PSI(i1),1,ZERO,ZX,1)
!       CALL ZSYMV('L',NDX,ZONER,ZONER*op,NDop,PSI(i1),1,ZERO,ZX,1)
!           CAL  ZSYMV( UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )
c       call  ZGEMV('N', NDX, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)
       call  ZGEMV('N', NDX,NDX,ZONER,ZONER*op,NDop,PSI(i1),1,ZERO,ZX,1)
       CALL ZCOPY(NDX,ZX,1,PSI(i1),1)
            Enddo
c                write(6,*) "X dim ==",PSI
         CASE (2)
             ZY=Zero
c                write(6,*) "Y dim =="
            Do I=1,NDX*NDY*NDZ,NDX*NDY
            Do J=1,NDX
            i1=I-1+J
c      CALL  ZTRMV('U','N','N',NDY,COp,NDop,psi(i1),NDX)
c       CALL ZCOPY(NDY,PSI(i1),NDX,ZY,1)
c                write(6,*) "X dim ==",ZY
c       CALL ZSYMV('L',NDY,ZONER,Cop,NDop,PSI(i1),NDX,ZERO,ZY,1)
!     CALL      ZSYMV('L',NDY,ZONER,ZONER*op,NDop,PSI(i1),NDX,ZERO,ZY,1)
      Call ZGEMV('N', NDY,NDY,ZONER,ZONER*op,NDop,PSI(i1),NDX,ZERO,ZY,1)
       CALL ZCOPY(NDY,ZY,1,PSI(i1),NDX)
c            pause
            Enddo
            Enddo
c                write(6,*) "Y dim ==",PSI
         CASE (3)
             ZZ=Zero
c               write(6,*) "Z dim =="
            Do i=1,NDX*NDY
            i1=i
c      CALL  ZTRMV('U','N','N',NDZ,COp,NDop,psi(i1),NDX*NDY)
c       CALL ZSYMV('L',NDZ,ZONER,Cop,NDop,PSI(i1),NDX*NDY,ZERO,ZZ,1)
!          CALL ZSYMV('L',NDZ,ZONER,ZONER*op,NDop,PSI(i1),NDX*NDY,ZERO,ZZ,1)
      Call ZGEMV('N', NDZ,NDZ,ZONER,ZONER*op,NDop,
     . PSI(i1),NDX*NDY,ZERO,ZZ,1)
      CALL ZCOPY(NDZ,ZZ,1,PSI(i1),NDX*NDY)
            Enddo
         END SELECT Xi1
c================================================================================
          END  SUBROUTINE Get_Op_PSI_Gnrl
