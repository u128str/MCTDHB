          SUBROUTINE Guess_Read_ORB(PSI)
          USE   SHARED_DIMS
          USE   DVR_ALL
          IMPLICIT NONE
c========================================================
       COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb)  :: PSI
       REAL*8, DIMENSION(NDX*NDY*NDZ,100)  :: vec
       REAL*8, DIMENSION(25*NDX*NDY*NDZ)  :: WORK
       INTEGER, DIMENSION (5*NDX*NDY*NDZ) :: IWORK
       INTEGER, DIMENSION   (NDX*NDY*NDZ) :: IFAIL
       REAL*8, DIMENSION(NDX*NDY*NDZ)  :: Vexternal
       REAL*8, DIMENSION(1000)  :: V
       REAL*8, DIMENSION(NDX*NDY*NDZ,NDX*NDY*NDZ) :: ham
       REAL*8, DIMENSION(NDX*NDY*NDZ) :: eigval
       INTEGER :: ndim,IL,IU,M,INFO,I
       REAL*8 :: VL,VU,ABSTOL
       character*26 lnamerest
       character*18 lname
       REAL*8 :: xlambda0_rest,efull1,xy,vec_a,vec_b
       REAL*8 :: x1,x2,x3,x4,x5,x6,x7,x8,x9
       INTEGER :: ndim_rest,ierr,jj,kk
       INTEGER :: J,K,L,ind,Nrestorb,Norb
c===============================================================
         PSI =ZERO
         ndim=NDX*NDY*NDZ
c===============================================================
             goto 2
         lnamerest='50.0000000psi1_psi2_LC.da'
         lnamerest='99.8746133_____psi2_LC.da'
       open(unit=09,file=lnamerest,Err=101,form='formatted')
       read(09,*,err=101,end=101) xlambda0_rest,ndim_rest,efull1
              if(ndim_rest.ne.ndim) then
              write(6,*)"ERRRRRRRRRR 1",ndim_rest, ndim
c              stop
              endif
       Do i=1, ndim
       read(09,*,err=101,end=101) xy,vec_a,vec_b
       PSI(i,1)=vec_a
       PSI(i,2)=vec_b
c       if(numberX(i).eq.ndimx) read(09,'(a2)') lnl
       enddo
       close(09)
c         write(6,*) "SCHMIDT ", ierr
c          pause
c===============================================================
1        continue
         Norb=Morb
c===============================================================
       lname='0.00000000time.dat'
       lname='external__file.dat'
       write(6,*)" Read Guess from the file",lname
       open(unit=11,file=lname,Err=101,form='formatted')
c       read(11,*)
       Do i=1, ndim
       read(11,*,err=101,end=101) x1,x2,x3,x4,x5,x6
!,x7,x8,x9
!       write(6,*) x1,x6
       PSI(i,1)=(0.d0,0d0)
       PSI(i,2)=(0.d0,0d0)
       PSI(i,1)=x2+(0.d0,1.0d0)*x3
       PSI(i,2)=x4+(0.d0,1.0d0)*x5
       enddo
       close(11)
       goto 9999 
c===============================================================
2        continue
         Norb=Morb
c         Norb=3
c===============================================================
      lname='IN.0000000time.daT'
c       lname='140.000000time.daT'
         lname=Time_Res_Orb_File_Name
c       write(6,*)" Read Guess from the file ",lname
      write(6,*)"Orbitals Guess is taken from file: ",lname
       open(unit=11,file=lname,Err=101,form='formatted')
       read(11,*,err=101,end=101)
       read(11,*,err=101,end=101)
       ind=1
       do K=1,NDZ
       do J=1,NDY
       do I=1,NDX
      read(11,*,end=101,err=101) (V(jj),jj=1,10+5*Norb)
c      write(6,*)" Read",(V(jj),jj=1,10+5*Norb)
c      pause
      Do kk=1,Norb
      PSI(ind,kk)=ZONER*V(10+2*(kk-1))+ZONEI*V(10+2*(kk-1)+1) !Working orbitals
c      PSI(ind,kk)=ZONER*V(10+2*Norb+2*(kk-1))
c     &           +ZONEI*V(10+2*Norb+2*(kk-1)+1) !Natural orbitals
      EndDo
c       ort_X(I),"  ",ort_Y(J),"  ",ort_Z(K),"  ",ww,"  "
c     $,DREAL(VTRAP_EXT(ind)),"  "
c     $,DREAL(rhoPR/xN),"   ",DIMAG(rhoPR/xN),"   "
c     $,DREAL(yy/ww/ww/xN),"   ",DIMAG(yy/ww/ww/xN),"   "
c     $,(DREAL(PSI(ind,jj)),"  ",DIMAG(PSI(ind,jj)),"  ",jj=1,Morb)
c     $,(DREAL(v(jj)),"  ",DIMAG(v(jj)),"   ",jj=1,Morb)
c     $,(Nocc(jj),"   ",jj=1,Morb)
c     $,t_total
       ind=ind+1
       enddo
       enddo
       enddo
       close(11)
c===============================================================
9999      continue
         write(6,*)"PSI-file has been successfully read"
c         write(6,*)"PSI",SUM(Conjg(PSI(:,1))*PSI(:,1))
c         write(6,*)"PSI",SUM(Conjg(PSI(:,1))*PSI(:,2))
c         write(6,*)"PSI",SUM(Conjg(PSI(:,2))*PSI(:,2))
       call  schmidtortho(PSI,NDX*NDY*NDZ,Norb,ierr)
         write(6,*)"PSI-file has been successfully Schmidted"
         return
101    write(6,*)"The input PSI DATA file is either empty or corrupted"
       write(6,*)"Change Guess"
           stop 
          END SUBROUTINE Guess_Read_Orb
c===============================================================
c===============================================================
c===============================================================
       SUBROUTINE Guess_Read_CI(VIN)
       USE   SHARED_DIMS
       USE   rR_hW
       USE   CI_All
       USE   W_INTERPARTICLE
       USE   DVR_ALL

      IMPLICIT NONE
c      implicit real*8 (a-h,o-z)
c      implicit integer (i-n)
      COMPLEX*16, DIMENSION(Nconf) :: VIN
      REAL*8, DIMENSION(1000)  :: v
      real*8 time,t_total,ww,CNK
      integer :: i,jj,k,Ncnf,Norb,ind,n,m,Nread,Mread,ind_new
      integer :: CI_READ
      character*18 lname
      integer, DIMENSION(100)   :: Ivec
      integer, DIMENSION(100)   :: Nvec
      integer :: GetInd
c==================================================
      VIN=Zero
c==================================================
      lname='0.50000000coef.dat'
      lname='18.0000000coef.dat'
      lname='2.00000000coef.dat'
      lname='1.30000000coef.daT'
      lname='5.00000000coef.daT'
      lname='40.0000000coef.daT'
      lname='1.00000000coef.daT'
      lname='50.0000000coef.daT'
      lname='7.00000000coef.daT'
      lname='50.GS00000coef.daT'
      lname='0.00000000coef.daT'
      lname='GS.III0000coef.daT'
      lname='GS.IV00000coef.daT'
c      lname='GS.VIII000coef.dat'
c      lname='7.50000000coef.daT'
c      lname='GS.V000000coef.daT'
c      lname='GS.IIL0000coef.daT'
c      lname='GS.II00000coef.daT'
c      lname='GS.IIHOT00coef.daT'
c      lname='GS.0000000coef.dat'
c      lname='GS.IITRN00coef.daT'
c      lname='GS.IIHOT00coef.daT'
c       lname='GS.0000000coef.daT'
c      lname='GS.III0000coef.daT'
c      lname='1.50000000coef.daT'
c      lname='6.00000000coef.daT'
      lname='44.0000000coef.daT'
      lname='140.000000coef.daT'
      lname='IN.0000000coef.daT'
      lname=Time_Res_CIc_File_Name
      write(6,*)"CI coeff Guess is taken from file: ",lname
c==================================================
       CI_READ=0 ! Default  Reading N=Npar M=Morb
c       CI_READ=2 ! Not normal reading N<=Npar M<=Morb
c==================================================
        CIRD  : SELECT CASE (CI_READ)
                  CASE (0)
c===================DEFAULT CASE N=Npar M=Morb ===============================
              N=Npar
              M=Morb
         Ncnf=NINT(CNK(N+M-1,N)) 
      write(6,*)" Data taken from similar system N=",N," Morb=",M
      open(unit=11,file=lname,Err=101,form='formatted')
      read(11,*,err=101,end=101)
      read(11,*,err=101,end=101)
      do I=1,Ncnf
      read(11,*,end=10,err=101) (V(jj),jj=1,3)
      VIN(NINT(V(1)))=ZONER*V(2)+ZONEI*V(3)
      EndDo
10      close(11)
                  CASE (2)
c===================Not normal reading from CASE N ne Npar or  M ne Morb ===============================
              N=Npar
              M=Morb-1
c              M=2
      IF(N.GT.Npar) STOP "READ CI: Npar  wrong"
      IF(M.GT.Morb) STOP "READ CI: M orb wrong"
      write(6,*)" Data taken from DIFFERENT system N=",N," M=",M
      open(unit=11,file=lname,Err=101,form='formatted')
      read(11,*,end=101,err=101)
      read(11,*,end=101,err=101)
         Ncnf=NINT(CNK(Npar+Morb-1,Npar)) 
      do I=1,Ncnf
      read(11,*,end=11,err=101) (V(jj),jj=1,3)
      ind=NINT(V(1))
      Ivec=0
      Nvec=0
      call GET_Ind_vs_ii(Ind,N,M,Ivec,Nvec)
      Ind_new=GetInd(Npar,Morb,Nvec)
c      write(6,*)"Reading:", Ind_new,V(1),Nvec
      VIN(Ind_new)=ZONER*V(2)+ZONEI*V(3)
      EndDo
11      close(11)
                END SELECT CIRD
      write(6,*)" CI guess Read is done",N," Morb=",M
          return
101    write(6,*)"The input CI DATA file is either empty or corrupted"
       write(6,*)"Change Guess"
        stop
      END SUBROUTINE Guess_Read_CI
