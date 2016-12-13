       subroutine GetCIJKL2body_56_OMP(MYID,VIN)
c       subroutine GetCIJKL1body_12_OMP(MYID,VIN)
c     . GetCIJKL1body_I(N,M,Nconf,VIN,ZHJK,ZRJK,TRM_REQ,Nterms)

c       subroutine GetCIJKL1body(VIN,HJK,Escale,ZRHO,P,N,M,Nconf)
       USE SHARED_DIMS
       USE rR_hW
       USE CI_All
       USE Parallel_CI
       USE omp_lib
       USE CI_prod
       implicit NONE
c=========================================================
       INTEGER ::  MYID,IPRC
c=========================================================
       integer::  P,N,M
       integer::  i,j,k,l,ii
       integer::  cI,cJ,cK,cL
c=========================================================
       integer::  I_current_term,Nterms
c=========================================================
        COMPLEX*16 :: VIN(:)
        COMPLEX*16, ALLOCATABLE  :: VOUT(:,:)
        integer :: TID,this_thread
        COMPLEX*16, DIMENSION(100) :: sumvout1=Zoner
        COMPLEX*16 :: zrho,Escale,zdir,zinv,zrhocnjg
c        INTEGER, DIMENSION(12) :: nvecin,nvecout,Sh_m,Sh_p
        real*8 :: xbprefac
c==============================================
       real*4::    start,finish,exec_time,cnk_time
c=====================================================================
       INTEGER ::  FromN,TillN

       INTEGER, DIMENSION(56) :: nvecin,nvecout,Sh_m,Sh_p
       COMPLEX*16, DIMENSION(1274406) :: RIJKL
       integer :: i1,j1
       integer :: i2,j2
       integer :: i3,j3
       integer :: i4,j4
       integer :: i5,j5
       integer :: i6,j6
       integer :: i7,j7
       integer :: i8,j8
       integer :: i9,j9
       integer :: i10,j10
       integer :: i11,j11
       integer :: i12,j12
       integer :: i13,j13
       integer :: i14,j14
       integer :: i15,j15
       integer :: i16,j16
       integer :: i17,j17
       integer :: i18,j18
       integer :: i19,j19
       integer :: i20,j20
       integer :: i21,j21
       integer :: i22,j22
       integer :: i23,j23
       integer :: i24,j24
       integer :: i25,j25
       integer :: i26,j26
       integer :: i27,j27
       integer :: i28,j28
       integer :: i29,j29
       integer :: i30,j30
       integer :: i31,j31
       integer :: i32,j32
       integer :: i33,j33
       integer :: i34,j34
       integer :: i35,j35
       integer :: i36,j36
       integer :: i37,j37
       integer :: i38,j38
       integer :: i39,j39
       integer :: i40,j40
       integer :: i41,j41
       integer :: i42,j42
       integer :: i43,j43
       integer :: i44,j44
       integer :: i45,j45
       integer :: i46,j46
       integer :: i47,j47
       integer :: i48,j48
       integer :: i49,j49
       integer :: i50,j50
       integer :: i51,j51
       integer :: i52,j52
       integer :: i53,j53
       integer :: i54,j54
       integer :: i55,j55
       integer :: i56,j56
        call cpu_time(start)
c==========================================================================================================================
              IPRC=MYID+1
c==========================================================================================================================
              FromN=CI_Proc_From(IPRC)
              TillN=CI_Proc_Till(IPRC)
c==========================================================================================================================
c        call cpu_time(start)
c======= MEMORY Allocation For Vout =====================================
        TID = OMP_GET_MAX_THREADS()
        ALLOCATE( VOUT(SIZE(VIN),TID) )
        IF(ALLOCATED(VOUT).eqv..FALSE.)
     . write(6,*)" MEM for VOUT is NOT ok!",SIZE(VOUT)
         VOUT=ZERO
         RIJKL=ZERO
c=====================================================================
c======================== Here all the required arrays and constants are READY
c=============================================
c=============== DO LOOP over the desired TRM_REQ
c============== Nterms is the number of the required terms
           L=0
!$OMP PARALLEL  DEFAULT(PRIVATE)
!$OMP& shared(VIN,VOUT,RIJKL,FromN,TillN,TERM_INDEX_2B,
!$OMP& TERM_REQ_2B,MCNK,RESCALE_2B,WIJKL,TID,
!$OMP& Ind_CI_2b,Prefactors_2b,CI_Prd,Npar,Morb)
        this_thread= OMP_GET_THREAD_NUM()+1
!$OMP DO  REDUCTION(+:RIJKL)
       Iloop: DO I=FromN,TillN
            I_current_term= TERM_REQ_2B(I)
        RIJKL(I_current_term)=ZERO
            P=TERM_INDEX_2B(I_current_term)
c================ Unpack cI cJ cK cL from P
        zrho=ZERO
        zrhocnjg=ZERO
c        ZRIJKL(I_current_term)=ZERO
        cL= INT(P/1000000)
        cK= INT((P-cL*1000000)/10000)
        cJ= INT((P-cL*1000000-cK*10000)/100)
        cI= P-cL*1000000-cK*10000-cJ*100
c==================== Formation of the Annihilation and Creation patterns which are shifts: Sh_m(minus) Sh_p(plus)
c============== From
             Sh_m=0
            Sh_m(cL)=Sh_m(cL)+1
            Sh_m(cK)=Sh_m(cK)+1
c============== To
            Sh_p=0
            Sh_p(cI)=Sh_p(cI)+1
            Sh_p(cJ)=Sh_p(cJ)+1

        Do i1=Npar+Morb-1-Sh_m(1),Morb-1,-1
        nvecin(1)=Npar+Morb-1-i1
        Do i2=i1-1-Sh_m(2)   ,Morb-2,-1
        nvecin(2)=i1-i2-1
        Do i3=i2-1-Sh_m(3)   ,Morb-3,-1
        nvecin(3)=i2-i3-1
        Do i4=i3-1-Sh_m(4)   ,Morb-4,-1
        nvecin(4)=i3-i4-1
        Do i5=i4-1-Sh_m(5)   ,Morb-5,-1
        nvecin(5)=i4-i5-1
        Do i6=i5-1-Sh_m(6)   ,Morb-6,-1
        nvecin(6)=i5-i6-1
        Do i7=i6-1-Sh_m(7)   ,Morb-7,-1
        nvecin(7)=i6-i7-1
        Do i8=i7-1-Sh_m(8)   ,Morb-8,-1
        nvecin(8)=i7-i8-1
        Do i9=i8-1-Sh_m(9)   ,Morb-9,-1
        nvecin(9)=i8-i9-1
        Do i10=i9-1-Sh_m(10)   ,Morb-10,-1
        nvecin(10)=i9-i10-1
        Do i11=i10-1-Sh_m(11)   ,Morb-11,-1
        nvecin(11)=i10-i11-1
        Do i12=i11-1-Sh_m(12)   ,Morb-12,-1
        nvecin(12)=i11-i12-1
        Do i13=i12-1-Sh_m(13)   ,Morb-13,-1
        nvecin(13)=i12-i13-1
        Do i14=i13-1-Sh_m(14)   ,Morb-14,-1
        nvecin(14)=i13-i14-1
        Do i15=i14-1-Sh_m(15)   ,Morb-15,-1
        nvecin(15)=i14-i15-1
        Do i16=i15-1-Sh_m(16)   ,Morb-16,-1
        nvecin(16)=i15-i16-1
        Do i17=i16-1-Sh_m(17)   ,Morb-17,-1
        nvecin(17)=i16-i17-1
        Do i18=i17-1-Sh_m(18)   ,Morb-18,-1
        nvecin(18)=i17-i18-1
        Do i19=i18-1-Sh_m(19)   ,Morb-19,-1
        nvecin(19)=i18-i19-1
        Do i20=i19-1-Sh_m(20)   ,Morb-20,-1
        nvecin(20)=i19-i20-1
        Do i21=i20-1-Sh_m(21)   ,Morb-21,-1
        nvecin(21)=i20-i21-1
        Do i22=i21-1-Sh_m(22)   ,Morb-22,-1
        nvecin(22)=i21-i22-1
        Do i23=i22-1-Sh_m(23)   ,Morb-23,-1
        nvecin(23)=i22-i23-1
        Do i24=i23-1-Sh_m(24)   ,Morb-24,-1
        nvecin(24)=i23-i24-1
        Do i25=i24-1-Sh_m(25)   ,Morb-25,-1
        nvecin(25)=i24-i25-1
        Do i26=i25-1-Sh_m(26)   ,Morb-26,-1
        nvecin(26)=i25-i26-1
        Do i27=i26-1-Sh_m(27)   ,Morb-27,-1
        nvecin(27)=i26-i27-1
        Do i28=i27-1-Sh_m(28)   ,Morb-28,-1
        nvecin(28)=i27-i28-1
        Do i29=i28-1-Sh_m(29)   ,Morb-29,-1
        nvecin(29)=i28-i29-1
        Do i30=i29-1-Sh_m(30)   ,Morb-30,-1
        nvecin(30)=i29-i30-1
        Do i31=i30-1-Sh_m(31)   ,Morb-31,-1
        nvecin(31)=i30-i31-1
        Do i32=i31-1-Sh_m(32)   ,Morb-32,-1
        nvecin(32)=i31-i32-1
        Do i33=i32-1-Sh_m(33)   ,Morb-33,-1
        nvecin(33)=i32-i33-1
        Do i34=i33-1-Sh_m(34)   ,Morb-34,-1
        nvecin(34)=i33-i34-1
        Do i35=i34-1-Sh_m(35)   ,Morb-35,-1
        nvecin(35)=i34-i35-1
        Do i36=i35-1-Sh_m(36)   ,Morb-36,-1
        nvecin(36)=i35-i36-1
        Do i37=i36-1-Sh_m(37)   ,Morb-37,-1
        nvecin(37)=i36-i37-1
        Do i38=i37-1-Sh_m(38)   ,Morb-38,-1
        nvecin(38)=i37-i38-1
        Do i39=i38-1-Sh_m(39)   ,Morb-39,-1
        nvecin(39)=i38-i39-1
        Do i40=i39-1-Sh_m(40)   ,Morb-40,-1
        nvecin(40)=i39-i40-1
        Do i41=i40-1-Sh_m(41)   ,Morb-41,-1
        nvecin(41)=i40-i41-1
        Do i42=i41-1-Sh_m(42)   ,Morb-42,-1
        nvecin(42)=i41-i42-1
        Do i43=i42-1-Sh_m(43)   ,Morb-43,-1
        nvecin(43)=i42-i43-1
        Do i44=i43-1-Sh_m(44)   ,Morb-44,-1
        nvecin(44)=i43-i44-1
        Do i45=i44-1-Sh_m(45)   ,Morb-45,-1
        nvecin(45)=i44-i45-1
        Do i46=i45-1-Sh_m(46)   ,Morb-46,-1
        nvecin(46)=i45-i46-1
        Do i47=i46-1-Sh_m(47)   ,Morb-47,-1
        nvecin(47)=i46-i47-1
        Do i48=i47-1-Sh_m(48)   ,Morb-48,-1
        nvecin(48)=i47-i48-1
        Do i49=i48-1-Sh_m(49)   ,Morb-49,-1
        nvecin(49)=i48-i49-1
        Do i50=i49-1-Sh_m(50)   ,Morb-50,-1
        nvecin(50)=i49-i50-1
        Do i51=i50-1-Sh_m(51)   ,Morb-51,-1
        nvecin(51)=i50-i51-1
        Do i52=i51-1-Sh_m(52)   ,Morb-52,-1
        nvecin(52)=i51-i52-1
        Do i53=i52-1-Sh_m(53)   ,Morb-53,-1
        nvecin(53)=i52-i53-1
        Do i54=i53-1-Sh_m(54)   ,Morb-54,-1
        nvecin(54)=i53-i54-1
        Do i55=i54-1-Sh_m(55) ,Morb-55+Sh_m(56),-1
        nvecin(55)=i54-i55-1
        nvecin(56)=i55-1
        K=1+MCNK(i1-1,Morb-1) 
     $  +MCNK(i2-1,Morb-2)   
     $  +MCNK(i3-1,Morb-3)   
     $  +MCNK(i4-1,Morb-4)   
     $  +MCNK(i5-1,Morb-5)   
     $  +MCNK(i6-1,Morb-6)   
     $  +MCNK(i7-1,Morb-7)   
     $  +MCNK(i8-1,Morb-8)   
     $  +MCNK(i9-1,Morb-9)   
     $  +MCNK(i10-1,Morb-10)   
     $  +MCNK(i11-1,Morb-11)   
     $  +MCNK(i12-1,Morb-12)   
     $  +MCNK(i13-1,Morb-13)   
     $  +MCNK(i14-1,Morb-14)   
     $  +MCNK(i15-1,Morb-15)   
     $  +MCNK(i16-1,Morb-16)   
     $  +MCNK(i17-1,Morb-17)   
     $  +MCNK(i18-1,Morb-18)   
     $  +MCNK(i19-1,Morb-19)   
     $  +MCNK(i20-1,Morb-20)   
     $  +MCNK(i21-1,Morb-21)   
     $  +MCNK(i22-1,Morb-22)   
     $  +MCNK(i23-1,Morb-23)   
     $  +MCNK(i24-1,Morb-24)   
     $  +MCNK(i25-1,Morb-25)   
     $  +MCNK(i26-1,Morb-26)   
     $  +MCNK(i27-1,Morb-27)   
     $  +MCNK(i28-1,Morb-28)   
     $  +MCNK(i29-1,Morb-29)   
     $  +MCNK(i30-1,Morb-30)   
     $  +MCNK(i31-1,Morb-31)   
     $  +MCNK(i32-1,Morb-32)   
     $  +MCNK(i33-1,Morb-33)   
     $  +MCNK(i34-1,Morb-34)   
     $  +MCNK(i35-1,Morb-35)   
     $  +MCNK(i36-1,Morb-36)   
     $  +MCNK(i37-1,Morb-37)   
     $  +MCNK(i38-1,Morb-38)   
     $  +MCNK(i39-1,Morb-39)   
     $  +MCNK(i40-1,Morb-40)   
     $  +MCNK(i41-1,Morb-41)   
     $  +MCNK(i42-1,Morb-42)   
     $  +MCNK(i43-1,Morb-43)   
     $  +MCNK(i44-1,Morb-44)   
     $  +MCNK(i45-1,Morb-45)   
     $  +MCNK(i46-1,Morb-46)   
     $  +MCNK(i47-1,Morb-47)   
     $  +MCNK(i48-1,Morb-48)   
     $  +MCNK(i49-1,Morb-49)   
     $  +MCNK(i50-1,Morb-50)   
     $  +MCNK(i51-1,Morb-51)   
     $  +MCNK(i52-1,Morb-52)   
     $  +MCNK(i53-1,Morb-53)   
     $  +MCNK(i54-1,Morb-54)   
     $  +MCNK(i55-1,Morb-55)   
        j1=Npar+Morb-1 - (nvecin(1) - Sh_m(1) + Sh_p(1))
         j2= j1-1-(nvecin(2) - Sh_m(2) + Sh_p(2))
         j3= j2-1-(nvecin(3) - Sh_m(3) + Sh_p(3))
         j4= j3-1-(nvecin(4) - Sh_m(4) + Sh_p(4))
         j5= j4-1-(nvecin(5) - Sh_m(5) + Sh_p(5))
         j6= j5-1-(nvecin(6) - Sh_m(6) + Sh_p(6))
         j7= j6-1-(nvecin(7) - Sh_m(7) + Sh_p(7))
         j8= j7-1-(nvecin(8) - Sh_m(8) + Sh_p(8))
         j9= j8-1-(nvecin(9) - Sh_m(9) + Sh_p(9))
         j10= j9-1-(nvecin(10) - Sh_m(10) + Sh_p(10))
         j11= j10-1-(nvecin(11) - Sh_m(11) + Sh_p(11))
         j12= j11-1-(nvecin(12) - Sh_m(12) + Sh_p(12))
         j13= j12-1-(nvecin(13) - Sh_m(13) + Sh_p(13))
         j14= j13-1-(nvecin(14) - Sh_m(14) + Sh_p(14))
         j15= j14-1-(nvecin(15) - Sh_m(15) + Sh_p(15))
         j16= j15-1-(nvecin(16) - Sh_m(16) + Sh_p(16))
         j17= j16-1-(nvecin(17) - Sh_m(17) + Sh_p(17))
         j18= j17-1-(nvecin(18) - Sh_m(18) + Sh_p(18))
         j19= j18-1-(nvecin(19) - Sh_m(19) + Sh_p(19))
         j20= j19-1-(nvecin(20) - Sh_m(20) + Sh_p(20))
         j21= j20-1-(nvecin(21) - Sh_m(21) + Sh_p(21))
         j22= j21-1-(nvecin(22) - Sh_m(22) + Sh_p(22))
         j23= j22-1-(nvecin(23) - Sh_m(23) + Sh_p(23))
         j24= j23-1-(nvecin(24) - Sh_m(24) + Sh_p(24))
         j25= j24-1-(nvecin(25) - Sh_m(25) + Sh_p(25))
         j26= j25-1-(nvecin(26) - Sh_m(26) + Sh_p(26))
         j27= j26-1-(nvecin(27) - Sh_m(27) + Sh_p(27))
         j28= j27-1-(nvecin(28) - Sh_m(28) + Sh_p(28))
         j29= j28-1-(nvecin(29) - Sh_m(29) + Sh_p(29))
         j30= j29-1-(nvecin(30) - Sh_m(30) + Sh_p(30))
         j31= j30-1-(nvecin(31) - Sh_m(31) + Sh_p(31))
         j32= j31-1-(nvecin(32) - Sh_m(32) + Sh_p(32))
         j33= j32-1-(nvecin(33) - Sh_m(33) + Sh_p(33))
         j34= j33-1-(nvecin(34) - Sh_m(34) + Sh_p(34))
         j35= j34-1-(nvecin(35) - Sh_m(35) + Sh_p(35))
         j36= j35-1-(nvecin(36) - Sh_m(36) + Sh_p(36))
         j37= j36-1-(nvecin(37) - Sh_m(37) + Sh_p(37))
         j38= j37-1-(nvecin(38) - Sh_m(38) + Sh_p(38))
         j39= j38-1-(nvecin(39) - Sh_m(39) + Sh_p(39))
         j40= j39-1-(nvecin(40) - Sh_m(40) + Sh_p(40))
         j41= j40-1-(nvecin(41) - Sh_m(41) + Sh_p(41))
         j42= j41-1-(nvecin(42) - Sh_m(42) + Sh_p(42))
         j43= j42-1-(nvecin(43) - Sh_m(43) + Sh_p(43))
         j44= j43-1-(nvecin(44) - Sh_m(44) + Sh_p(44))
         j45= j44-1-(nvecin(45) - Sh_m(45) + Sh_p(45))
         j46= j45-1-(nvecin(46) - Sh_m(46) + Sh_p(46))
         j47= j46-1-(nvecin(47) - Sh_m(47) + Sh_p(47))
         j48= j47-1-(nvecin(48) - Sh_m(48) + Sh_p(48))
         j49= j48-1-(nvecin(49) - Sh_m(49) + Sh_p(49))
         j50= j49-1-(nvecin(50) - Sh_m(50) + Sh_p(50))
         j51= j50-1-(nvecin(51) - Sh_m(51) + Sh_p(51))
         j52= j51-1-(nvecin(52) - Sh_m(52) + Sh_p(52))
         j53= j52-1-(nvecin(53) - Sh_m(53) + Sh_p(53))
         j54= j53-1-(nvecin(54) - Sh_m(54) + Sh_p(54))
         j55= j54-1-(nvecin(55) - Sh_m(55) + Sh_p(55))
        J=1+MCNK(j1-1,Morb-1) 
     $  +MCNK(j2-1,Morb-2)   
     $  +MCNK(j3-1,Morb-3)   
     $  +MCNK(j4-1,Morb-4)   
     $  +MCNK(j5-1,Morb-5)   
     $  +MCNK(j6-1,Morb-6)   
     $  +MCNK(j7-1,Morb-7)   
     $  +MCNK(j8-1,Morb-8)   
     $  +MCNK(j9-1,Morb-9)   
     $  +MCNK(j10-1,Morb-10)   
     $  +MCNK(j11-1,Morb-11)   
     $  +MCNK(j12-1,Morb-12)   
     $  +MCNK(j13-1,Morb-13)   
     $  +MCNK(j14-1,Morb-14)   
     $  +MCNK(j15-1,Morb-15)   
     $  +MCNK(j16-1,Morb-16)   
     $  +MCNK(j17-1,Morb-17)   
     $  +MCNK(j18-1,Morb-18)   
     $  +MCNK(j19-1,Morb-19)   
     $  +MCNK(j20-1,Morb-20)   
     $  +MCNK(j21-1,Morb-21)   
     $  +MCNK(j22-1,Morb-22)   
     $  +MCNK(j23-1,Morb-23)   
     $  +MCNK(j24-1,Morb-24)   
     $  +MCNK(j25-1,Morb-25)   
     $  +MCNK(j26-1,Morb-26)   
     $  +MCNK(j27-1,Morb-27)   
     $  +MCNK(j28-1,Morb-28)   
     $  +MCNK(j29-1,Morb-29)   
     $  +MCNK(j30-1,Morb-30)   
     $  +MCNK(j31-1,Morb-31)   
     $  +MCNK(j32-1,Morb-32)   
     $  +MCNK(j33-1,Morb-33)   
     $  +MCNK(j34-1,Morb-34)   
     $  +MCNK(j35-1,Morb-35)   
     $  +MCNK(j36-1,Morb-36)   
     $  +MCNK(j37-1,Morb-37)   
     $  +MCNK(j38-1,Morb-38)   
     $  +MCNK(j39-1,Morb-39)   
     $  +MCNK(j40-1,Morb-40)   
     $  +MCNK(j41-1,Morb-41)   
     $  +MCNK(j42-1,Morb-42)   
     $  +MCNK(j43-1,Morb-43)   
     $  +MCNK(j44-1,Morb-44)   
     $  +MCNK(j45-1,Morb-45)   
     $  +MCNK(j46-1,Morb-46)   
     $  +MCNK(j47-1,Morb-47)   
     $  +MCNK(j48-1,Morb-48)   
     $  +MCNK(j49-1,Morb-49)   
     $  +MCNK(j50-1,Morb-50)   
     $  +MCNK(j51-1,Morb-51)   
     $  +MCNK(j52-1,Morb-52)   
     $  +MCNK(j53-1,Morb-53)   
     $  +MCNK(j54-1,Morb-54)   
     $  +MCNK(j55-1,Morb-55)   
c====================== Bosonic prefactors
          nvecout=nvecin
c          IF(nvecin(4).eq.24.and.nvecin(6).eq.1)
c     .    write(6,*) nvecin,": 0 0 0 24 0 1 0 ",K
c          stop
c===================== cI,cJ,cK,cL two-body
             xbprefac=1.0d0
             xbprefac=xbprefac*nvecout(cL)
             nvecout(cL)=nvecout(cL)-1
             xbprefac=xbprefac*nvecout(cK)
             nvecout(cK)=nvecout(cK)-1
             nvecout(cJ)=nvecout(cJ)+1
             xbprefac=xbprefac*nvecout(cJ)
             nvecout(cI)=nvecout(cI)+1
             xbprefac=xbprefac*nvecout(cI)
             IF(CI_Prd.eqv..True.) THEN
             Ind_CI_2b(K,I-FromN+1)=J
             Prefactors_2b(K,I-FromN+1)=xbprefac
             ENDIF
c=========================== Normal    C^_I C^_J C_K C_L  rho_IJ_KL
           zdir=VIN(K)* SQRT(xbprefac)
c           VOUT(J)=zdir
           zrho=zrho+Conjg(VIN(J))*zdir
        RIJKL(I_current_term)=RIJKL(I_current_term)+Conjg(VIN(J))*zdir
c=========================== Conjugate C^_K C^_L C_I C_J   rho_IJ_KL=rho^*_KL_IJ
           zinv=VIN(J)* SQRT(xbprefac)
c           VOUT(K)=zinv
           zrhocnjg=zrhocnjg+Conjg(VIN(K))*zinv
c============================= Explicit construction of
c================== |VOUT>=(Re(Escale)*Wijkl *C^_I C^_J C_K C_L + Im(Escale)*Wklij*C^_K C^_L C_I C_J) |PSI>
c=============================  Re(Escale), Im(Escale) can be IIII (0.5:0.5), IJIJ (2:2), IJKL (4:4), IIKL(2:2), IILL (1:1)
        VOUT(J,this_thread)=VOUT(J,this_thread)
     .+zdir*       WIJKL(I_current_term)
     .                      *DREAL(RESCALE_2B(I_current_term))
        VOUT(K,this_thread)=VOUT(K,this_thread)
     .+zinv* Conjg(WIJKL(I_current_term))
     .                      *DIMAG(RESCALE_2B(I_current_term))
c=============================================================
c            L=L+1
c==============================

c       write(6,'(i7,i7,a10,7(i6,a),a10,i7,a7,7(i6,a))') K,L,
c     $" Bosonic:|",(nvecin(i)," ", i=1,7),
c     $" |-----> ",J," | ", (nvecout(i)," ", i=1,7)


            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
            EndDo
                  EndDO Iloop
!$OMP END DO
!$OMP END PARALLEL
           ZRIJKL=RIJKL
      CALL ZGEMV('N',Nconf,TID,Zoner,VOUT,Nconf,sumvout1,1,ZERO,VIN,1)

c        VIN=VOUT
       DEALLOCATE(VOUT)
c       CI_Production_2b=.TRUE.
       CI_Production_2b=CI_Prd
c       write(6,*)MYID," CIJKL_2B  DONE",FromN,TillN
c       call cpu_time(finish)
c       exec_time=finish-start
c       write(6,*)MYID,"IN 2-Bbody: Execution DNS time:",exec_time,L
c================================================================
c          call cpu_time(finish)
c          exec_time=finish-start
c          write(6,*)IPRC,TillN-FromN,"CPU Get 2b",exec_time



        end subroutine GetCIJKL2body_56_OMP
