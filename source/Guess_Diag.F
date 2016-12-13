          SUBROUTINE Guess_Diag(PSI,Vext)
          USE   SHARED_DIMS
          USE   DVR_ALL
       USE Parallel_Orb
!       USE rR_hW
          IMPLICIT NONE
       INTEGER :: IL,IU,M,INFO,I,J,i1,i2,j1,j2,ind,q1,q2
       REAL*8 :: VL,VU,ABSTOL
c       INTEGER, PARAMETER :: ndim=NDX*NDY*NDZ
c       INTEGER, PARAMETER :: LWORK=8*NDX*NDY*NDZ
       REAL*8, DIMENSION(8*NDX*NDY*NDZ)  :: WORK
       INTEGER :: ndim ! =NDX*NDY*NDZ
       INTEGER :: LWORK! LWORK=8*NDX*NDY*NDZ
c========================================================
       COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb)  :: PSI
       REAL*8, DIMENSION(NDX*NDY*NDZ)  :: Vext
       REAL*8, DIMENSION(NDX*NDY*NDZ,Morb+4)  :: vec
Corg       REAL*8, DIMENSION(LWORK)  :: WORK
       INTEGER, DIMENSION (5*NDX*NDY*NDZ) :: IWORK
       INTEGER, DIMENSION   (NDX*NDY*NDZ) :: IFAIL
       REAL*8, DIMENSION(NDX*NDY*NDZ,NDX*NDY*NDZ) :: ham
       COMPLEX*16, DIMENSION(NDX,NDX) :: hamZ
       REAL*8, DIMENSION(NDX*NDY*NDZ) :: eigval

       ndim=NDX*NDY*NDZ
       LWORK=8*NDX*NDY*NDZ




c===============================================================
         PSI =ZERO
         VEC=0d0
         eigval=0d0
         WORK=0d0
         IWORK=0
         IFAIL=0
c===============================================================
      Dim_1_2_3D: SELECT CASE (DIM_MCTDHB)

       CASE (1) ! 1D  ! 1D  ! 1D  ! 1D  ! 1D  ! 1D  ! 1D  ! 1D  ! 1D
              write(6,*)" DIAG GUESS",Time_DVRMETHODX
          IF(Time_DVRMETHODX==4) then !AIS 22JUN2014 T FFT is added
           call T_FFT_1D(hamZ)
            ham=REAL(hamZ)
          else
                 ham=Op_X
          endif
!         ham=Op_X
         Do I=1,ndim
         ham(i,i)=ham(i,i)+Vext(i)
         EnD DO
       CASE (2) ! 2D  2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D
         write(6,*)"2D case in Guess_Diag"
c         ham=Op_X+Op_Y+Vext

         Do j1=1,NDY
         Do j2=1,NDY
         IF(j1.eq.j2) then
         q1=1
         else
         q1=0
         endIF
            Do i1=1,NDX
            Do i2=1,NDX
            IF(i1.eq.i2) then
            q2=1
            else
            q2=0
            endIF
            I=NDX*(j1-1)+i1
            J=NDX*(j2-1)+i2
            ham(I,J)=q1*Op_X(i1,i2)+q2*Op_Y(j1,j2)
        IF(I.eq.J) ham(I,I)=ham(I,I)+Vext(I)
c      write(6,'(i2,i2,a10,4i2)') I,J,"i1i2_j1j2",q1*i1,q1*i2,q2*j1,q2*j2
            endDo
            endDo
         endDo
         endDo
       CASE (3) ! 3D ! 3D ! 3D ! 3D ! 3D ! 3D ! 3D ! 3D ! 3D ! 3D ! 3D
       write(*,*)"3D case in Guess_Diag"
       return 
       END SELECT Dim_1_2_3D
c===============================================================
        VL=0.0d0
        VU=1.0d0
        IL=1
        IU=MAX(Morb+4,IL+1)
        M=IU
        ABSTOL=0.0D0
         write(6,*)"Egenvalues from external potential!!!!"
c          pause
c         write(6,*) size(ham),size(Op_X)
c         write(6,*) Morb,ndim,VL,VU,IL,IU,ABSTOL
c         write(6,*) ham(1,1)
c         write(6,*) eigval(1)
c         write(6,*) vec(1,1)
c         write(6,*) work(1)
c         write(6,*) iwork(1)
c===============================================================
C       SUBROUTINE DSYEVX( JOBZ, RANGE, UPLO,  N,  A,  LDA,  VL,  VU,  IL,  IU,
C                          ABSTOL,  M,  W,  Z,  LDZ, WORK, LWORK, IWORK, IFAIL,
C                          INFO )
c      SUBROUTINE DSYEVR( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU,
c     $                   ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK,
c     $                   IWORK, LIWORK, INFO )
        CALL DSYEVX('N','I','U',ndim,ham,ndim,VL,VU,IL,IU,ABSTOL,
     .  M,eigval,vec, ndim, WORK, LWORK , IWORK, IFAIL, INFO)
         write(6,'(30(F10.5,a3))')(eigval(i),"   ",i=1,M)
          open(unit=10,file='Init_Trap_Eigen.dat',form='formatted')
          Do i=1,M
          write(10,'((a2,i2,a2,F20.16))')"E(",i,")=",eigval(i)
          Enddo
          close(10)
c            stop
c===============================================================
      Dim_1_2_3D_copy: SELECT CASE (DIM_MCTDHB)

       CASE (1) ! 1D  ! 1D  ! 1D  ! 1D  ! 1D  ! 1D  ! 1D  ! 1D  ! 1D
          IF(Time_DVRMETHODX==4) then !AIS 22JUN2014 T FFT is added
         call T_FFT_1D(hamZ)
           ham=REAL(hamZ)
          else
                 ham=Op_X
          endif
!         ham=Op_X
         Do I=1,ndim
         ham(i,i)=ham(i,i)+Vext(i)
         EnD DO
       CASE (2) ! 2D  2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D 2D
         write(6,*)"2D case in Guess_Diag"
c         ham=Op_X+Op_Y+Vext

         Do j1=1,NDY
         Do j2=1,NDY
         IF(j1.eq.j2) then
         q1=1
         else
         q1=0
         endIF
            Do i1=1,NDX
            Do i2=1,NDX
            IF(i1.eq.i2) then
            q2=1
            else
            q2=0
            endIF
            I=NDX*(j1-1)+i1
            J=NDX*(j2-1)+i2
            ham(I,J)=q1*Op_X(i1,i2)+q2*Op_Y(j1,j2)
        IF(I.eq.J) ham(I,I)=ham(I,I)+Vext(I)
c      write(6,'(i2,i2,a10,4i2)') I,J,"i1i2_j1j2",q1*i1,q1*i2,q2*j1,q2*j2
            endDo
            endDo
         endDo
         endDo
       CASE (3) ! 3D ! 3D ! 3D ! 3D ! 3D ! 3D ! 3D ! 3D ! 3D ! 3D ! 3D
       write(*,*)"3D case in Guess_Diag"
       return 
       END SELECT Dim_1_2_3D_copy
c===============================================================

        CALL DSYEVX('V','I','U',ndim,ham,ndim,VL,VU,IL,IU,ABSTOL,
     .  M,eigval,vec, ndim, WORK, LWORK , IWORK, IFAIL, INFO)
c            pause

c        CALL DSYEV( 'N','U',ndim, ham, ndim, eigval,WORK,LWORK,INFO )
cc           DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
c         CALL DCOPY (ndim,weigval, INCX, weigval1, INCY)
C==================== writing the initiall guesses for VEC, VEC1
         write(6,*)M,"Initial egenvalues are:",SIZE(eigval)
         write(6,'(30(F10.5,a3))')(eigval(i),"   ",i=1,Morb+4)
c         pause

         Do I=1,Morb
         PSI(:,I)=vec(:,I)*(1.0d0,0d0)
c         PSI(:,I)=ham(:,I)*Zoner
         END DO
c         call  schmidtortho(PSI,NDX*NDY*NDZ,2,ierr)
c         write(6,*) "SCHMIDT ", ierr
c          pause

          END SUBROUTINE Guess_Diag 
