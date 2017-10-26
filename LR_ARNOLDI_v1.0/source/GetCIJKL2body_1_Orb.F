C=============== Input: C^IC^JC_KC_L as a packed INTEGER number P=I+100*J+10000*K+1000000*L
C===============        VECIN - incoming complex*16 vector
C=============== Output: VECOUT - outcoming complex*16 vector: VECOUT=C^IC^JC_KC_L*VecIN
C===============        < VECIN^| C^IC^JC_KC_L*VecIN>== rho_ijkl
C=============== Lexicografical Ordering and numberging -> inverse combinadic (only zero's are adressed) 
cORG       subroutine GetCIJKL2body_I(VIN,WIJKL,Escale,ZRHO,P,N,M,Nconf)
c       subroutine  .GetCIJKL2body_Par(N,M,Nconf,VIN,ZWIJKL,ZRIJKL,TRM_REQ,FromN,TillN)
       subroutine  GetCIJKL2body_1(MYID,VIN)
       USE SHARED_DIMS
       USE rR_hW
       USE CI_All
       USE Parallel_CI
       implicit NONE
c=========================================================
       INTEGER ::  MYID,IPRC
c=========================================================
c=================== MPI
c       INCLUDE 'mpif.h'
c       INTEGER ::  ierr,MYID,numprocs
c=========================================================
       integer :: P,N,M 
       integer :: i,j,jcnjg,k,l
       integer :: i1,i2,i3,i4,i5,i6,i7,i8,i9,i10
       integer :: j1,j2,j3,j4,j5,j6,j7,j8,j9,j10
       integer :: n1,n2,n3,n4,n5,n6,n7,n8,n9,n10
       integer :: m1,m2,m3,m4,m5,m6,m7,m8,m9,m10
       integer :: k1,k2
       integer :: Sh_m1,Sh_m2,Sh_m3,Sh_m4,Sh_m5
       integer :: Sh_m6,Sh_m7,Sh_m8,Sh_m9,Sh_m10
       integer :: Sh_p1,Sh_p2,Sh_p3,Sh_p4,Sh_p5
       integer :: Sh_p6,Sh_p7,Sh_p8,Sh_p9,Sh_p10
       integer :: Sh_k,Sh_l,jst

       integer :: cI,cJ,cK,cL,Nadr,Wadress,NadrW
       INTEGER ::  FromN,TillN

c=========================================================
c        COMPLEX*16, DIMENSION(Nconf) :: VIN,VOUT
        COMPLEX*16 :: VIN(:)
c        COMPLEX*16, ALLOCATABLE  :: VOUT(:)
c        COMPLEX*16, DIMENSION(Nconf) :: VOUT

        COMPLEX*16 zrho,zrhocnjg
        COMPLEX*16 zdir,zinv
        INTEGER   , DIMENSION(7) :: nvecin,nvecout,Sh_m,Sh_p
      integer :: I_current_term,Nterms
      real*8 ::  xbprefac        
      real*4 ::  start,finish, exec_time , cnk_time       
c==========================================================================================================================
c==========================================================================================================================
              IPRC=MYID+1
c==========================================================================================================================
              FromN=CI_Proc_From(IPRC)
              TillN=CI_Proc_Till(IPRC)
c==========================================================================================================================
c         write(6,*)" CIJKL_2B  ",FromN,TillN
c        call cpu_time(start)
c==========================================================================================================================
c==========================================================================================================================
c======= MEMORY Allocation For Vout =====================================
c        ALLOCATE( VOUT(SIZE(VIN)) )
c        IF(ALLOCATED(VOUT).eqv..FALSE.) 
c     . write(6,*)" MEM for VOUT is NOT ok!",SIZE(VOUT)
c         VOUT=ZERO
c=====================================================================
c======================== Here all the required arrays and constants are READY
c=============================================
        ZRIJKL(1)=ZONER*Npar*(Npar-1)
        VIN(1)=ZONER*Npar*(Npar-1)*WIJKL(1)*DREAL(RESCALE_2B(1))
c        VIN=VOUT
c         write(6,*)" CIJKL_2B  DONE"
      end subroutine GetCIJKL2body_1
