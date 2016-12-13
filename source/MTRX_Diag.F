          SUBROUTINE MTRX_Diag
          USE   SHARED_DIMS
          USE   DVR_ALL
          USE   rR_hW
          IMPLICIT NONE
       INTEGER :: ndim ! =NDX*NDY*NDZ
       INTEGER :: IL,IU,M,INFO,I,NDM,indx
       REAL*8 :: VL,VU,ABSTOL
c       INTEGER, PARAMETER :: ndim=NDX*NDY*NDZ
c       INTEGER, PARAMETER :: LWORK=8*NDX*NDY*NDZ
       INTEGER :: LWORK! LWORK=8*NDX*NDY*NDZ
       INTEGER :: cL,cK,cJ,cI,P,J
c========================================================
       REAL*8, DIMENSION(Rdim,Rdim)  :: vec
       REAL*8, DIMENSION(8*Rdim)  :: WORK
       INTEGER, DIMENSION (5*Rdim) :: IWORK
       INTEGER, DIMENSION   (Rdim) :: IFAIL
       REAL*8, DIMENSION(Rdim,Rdim) :: MTRX
       REAL*8, DIMENSION(Rdim) :: eigval

       LWORK=8*Rdim
       ndim=Rdim
c===============================================================
         VEC=0d0
         MTRX=0d0
         eigval=0d0
         WORK=0d0
         IWORK=0
         IFAIL=0
c===============================================================
         indx=0
         Do I=1,ndim
         Do J=I,ndim
         indx=indx+1
         MTRX(I,J)=WIJKL(indx)
         P=TERM_INDEX_2B(indx)
c================ Unpack cI cJ cK cL from P
         cL= INT(P/1000000)
         cK= INT((P-cL*1000000)/10000)
         cJ= INT((P-cL*1000000-cK*10000)/100)
         cI= P-cL*1000000-cK*10000-cJ*100 
         IF(cL.eq.cJ.and.cK.eq.cI)  THEN
c         write(6,'(4i3,a4,F26.16,a4,F26.16,a4,F26.16)') 
c     .  cI,cJ,cK,CL
c     ,"W",Real(WIJKL(I)), ."Re:",DReal(ZRIJKL(I)),"Im:",DImag(ZRIJKL(I))
c        MTRX(I,J)=MTRX(I,J)+Real(AllHIJ(cL,cL))+Real(AllHIJ(cK,cK))
        MTRX(I,J)=MTRX(I,J)/2
            ENDIF
        MTRX(I,J)=MTRX(I,J)+Real(AllHIJ(cL,cJ))+Real(AllHIJ(cK,cI))
c         write(101,'(4i3,a4,F26.16,a4,F26.16,a4,F26.16)') 
c     .  cI,cJ,cK,CL,"W",Real(WIJKL(I)), ."Re:",DReal(ZRIJKL(I)),"Im:",DImag(ZRIJKL(I))
         EnD DO
         EnD DO
c===============================================================
        VL=0.0d0
        VU=1.0d0
        IL=1
        IU=MAX(Morb+4,IL)
        IU=10
        M=IU
        M=10
        ABSTOL=0.0D0
c         write(6,*)"Full CI Analitic!!!!"
c         write(6,*) size(ham),size(Op_X)
c         write(6,*) Morb,ndim,VL,VU,IL,IU,ABSTOL
          do i=1,4
         write(6,'(4F10.5)') (Real(AllHIJ(i,j)),j=1,4)
           enddo
c          do i=1,10
c         write(6,'(10F10.5)') (MTRX(i,J),j=1,10)
c           enddo
c         write(6,*) eigval(1)
c         write(6,*) vec(1,1)
c         write(6,*) work(1)
c         write(6,*) iwork(1)
c              MTRX=0d0
c          do i=1,4
c          do j=1,4
c           MTRX(i,J)= Real(AllHIJ(i,j))
c           enddo
c           enddo
c           IU=4
c===============================================================
C       SUBROUTINE DSYEVX( JOBZ, RANGE, UPLO,  N,  A,  LDA,  VL,  VU,  IL,  IU,
C                          ABSTOL,  M,  W,  Z,  LDZ, WORK, LWORK, IWORK, IFAIL,
C                          INFO )
c      SUBROUTINE DSYEVR( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU,
c     $                   ABSTOL, M, W, Z, LDZ, ISUPPZ, WORK, LWORK,
c     $                   IWORK, LIWORK, INFO )
c===============================================================
        CALL DSYEVX('N','I','U',ndim,MTRX,ndim,VL,VU,IL,IU,ABSTOL,
     .  M,eigval,vec, ndim, WORK, LWORK , IWORK, IFAIL, INFO)
c            pause

c        CALL DSYEV( 'N','U',ndim, ham, ndim, eigval,WORK,LWORK,INFO )
cc           DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
c         CALL DCOPY (ndim,weigval, INCX, weigval1, INCY)
C==================== writing the initiall guesses for VEC, VEC1
         write(6,*)M,"FCI egenvalues are:",SIZE(eigval)
         write(6,'(20(F10.5,a3))')(eigval(i),"   ",i=1,10)
c         pause

c         Do I=1,Morb
c         PSI(:,I)=vec(:,I)*Zoner
c         PSI(:,I)=ham(:,I)*Zoner
c         END DO
c         call  schmidtortho(PSI,NDX*NDY*NDZ,2,ierr)
c         write(6,*) "SCHMIDT ", ierr
c          pause

          END SUBROUTINE MTRX_Diag 
