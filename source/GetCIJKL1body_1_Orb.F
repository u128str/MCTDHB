C=============== Input: C^IC^JC_KC_L as a packed INTEGER number P=I+100*J+10000*K+1000000*L
C===============        VECIN - incoming complex*16 vector
C=============== Output: VECOUT - outcoming complex*16 vector: VECOUT=C^IC^JC_KC_L*VecIN
C===============        < VECIN^| C^IC^JC_KC_L*VecIN>== rho_ijkl
C=============== Lexicografical Ordering and numberging -> inverse combinadic (only zero's are adressed) 
       subroutine GetCIJKL1body_1(MYID,VIN)
c     . GetCIJKL1body_I(N,M,Nconf,VIN,ZHJK,ZRJK,TRM_REQ,Nterms)

c       subroutine GetCIJKL1body(VIN,HJK,Escale,ZRHO,P,N,M,Nconf)
       USE SHARED_DIMS
       USE rR_hW
       USE CI_All
       implicit NONE
c=========================================================
       INTEGER ::  MYID,IPRC
c=========================================================
       integer::  P,N,M 
       integer::  i,j,k,l
       integer::  i1,i2,i3,i4,i5,i6,i7,i8,i9,i10
       integer::  j1,j2,j3,j4,j5,j6,j7,j8,j9,j10
       integer::  n1,n2,n3,n4,n5,n6,n7,n8,n9,n10
       integer::  m1,m2,m3,m4,m5,m6,m7,m8,m9,m10
       integer::  Sh_m1,Sh_m2,Sh_m3,Sh_m4,Sh_m5
       integer::  Sh_m6,Sh_m7,Sh_m8,Sh_m9,Sh_m10
       integer::  Sh_p1,Sh_p2,Sh_p3,Sh_p4,Sh_p5
       integer::  Sh_p6,Sh_p7,Sh_p8,Sh_p9,Sh_p10
       integer::  I_current_term,Nterms
       integer::  Sh_k,Sh_l

       integer::  cI,cJ,cK,cL
c=========================================================
c        COMPLEX*16, DIMENSION(Nconf) :: VIN,VOUT
        COMPLEX*16 :: VIN(:)

        COMPLEX*16 :: zrho,Escale,zdir,zinv,zrhocnjg
        INTEGER, DIMENSION(7) :: nvecin,nvecout,Sh_m,Sh_p
        real*8 :: xbprefac
c==============================================
       real*4::    start,finish,exec_time,cnk_time       
c=====================================================================
c=====================================================================
         Nterms=MaxTrm1b
c=====================================================================
c=====================================================================
                IPRC=MYID+1           
        ZRIJ(1)=ZONER*Npar
        VIN(1)=ZONER*Npar*HIJ(1)
c        VIN=VOUT
c           write(6,*) " 1B done "
c        pause
      end subroutine GetCIJKL1body_1
