        subroutine GetCIJKL1body_27_OMP(MYID,VIN)
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

       INTEGER, DIMENSION(27) :: nvecin,nvecout,Sh_m,Sh_p
       COMPLEX*16, DIMENSION(378) :: RIJ
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
           call cpu_time(start)
       IPRC=MYID+1
       FromN=CI_Proc_From(IPRC)
       TillN=CI_Proc_Till(IPRC)
c=====================================================================
         Nterms=MaxTrm1b
c=====================================================================
c=====================================================================
                IPRC=MYID+1
c       write(6,*)" 1B ID",IPRC," is doing ",Nterms
c======= MEMORY Allocation =====================================
c       write(6,*)" MEM for  VIN 1B is ok!",SIZE(VIN)
        TID = OMP_GET_MAX_THREADS()
        ALLOCATE( VOUT(SIZE(VIN),TID) )
        IF(ALLOCATED(VOUT).eqv..FALSE.)
     . write(6,*)" MEM for VOUT is NOT ok!",SIZE(VOUT)
            VOUT=ZERO
            RIJ=ZERO
c=============== DO LOOP over the desired TRM_REQ
c============== Nterms is the number of the required terms
           L=0
!$OMP PARALLEL  DEFAULT(PRIVATE)
!$OMP& shared(VIN,VOUT,RIJ,Nterms,TERM_INDEX_1B,
!$OMP& TERM_REQ_1B,MCNK,RESCALE_1B,HIJ,TID,FromN,TillN,
!$OMP& Ind_CI_1b,Prefactors_1b,CI_Prd,Npar,Morb)
        this_thread= OMP_GET_THREAD_NUM()+1
!$OMP DO  REDUCTION(+:RIJ)
           DO I=FromN,TillN
           I_current_term= TERM_REQ_1B(I)
           P=TERM_INDEX_1B(I_current_term)
c================ Unpack cI cJ cK cL from P
        zrho=ZERO
        zrhocnjg=ZERO
        RIJ(I_current_term)=ZERO
        cK= INT(P/100)
        cJ= P-cK*100
c         write(6,*)"cJ-cK:",cJ,cK
c         write(6,*)"cJ-cK:",cJ,cK,ZHJK(I_current_term),RESCALE(I)
c==================== Annihilation and Creation
c============== From
            Sh_m=0
            Sh_m(cK)=Sh_m(cK)+1
c============== To
            Sh_p=0
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
        Do i26=i25-1-Sh_m(26) ,Morb-26+Sh_m(27),-1
        nvecin(26)=i25-i26-1
        nvecin(27)=i26-1
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
c====================== Bosonic prefactors
             nvecout=nvecin
c           write(6,*)K," n ",nvecin
c===================== cJ,cK one-body
             xbprefac=1.0d0
             xbprefac=xbprefac*nvecout(cK)
             nvecout(cK)=nvecout(cK)-1
             nvecout(cJ)=nvecout(cJ)+1
             xbprefac=xbprefac*nvecout(cJ)

c             VOUT(J)=VIN(K)* SQRT(xbprefac)
c             zrho=zrho+Conjg(VIN(J))*VOUT(J)
              IF(CI_Prd.eqv..True.) THEN
              Ind_CI_1b(K,I-FromN+1)=J
              Prefactors_1b(K,I-FromN+1)=xbprefac
              ENDIF
c=========================== Normal    C^_0 C^_J C_K C_0  rho_0J_K0
           zdir=VIN(K)* DSQRT(xbprefac)

c           VOUT(J)=zdir
           zrho=zrho+Conjg(VIN(J))*zdir
      RIJ(I_current_term)=RIJ(I_current_term)+Conjg(VIN(J))*zdir
c=========================== Conjugate C^_K C^_0 C_0 C_J   rho_0J_K0=rho^*_K0_0J
           zinv=VIN(J)* DSQRT(xbprefac)
c           VOUT(K)=zinv
           zrhocnjg=zrhocnjg+Conjg(VIN(K))*zinv
c============================= Explicit construction of
c================== |VOUT>=(Re(Escale)*Hjk *C^_0 C^_J C_K C_0 + Im(Escale)*Hkj*C^_K C^_0 C_0 C_J) |PSI>
c=============================  Re(Escale), Im(Escale) can be IIII (0.5:0.5), IL (1:1)
c              VOUT(J)=VOUT(J)+zdir*HJK*DREAL(Escale)
c              VOUT(K)=VOUT(K)+zinv*HJK*DIMAG(Escale)
        VOUT(J,this_thread)=VOUT(J,this_thread)
     .+zdir*   HIJ(I_current_term)
     .                      *DREAL(RESCALE_1B(I_current_term))
        VOUT(K,this_thread)=VOUT(K,this_thread)
     .+zinv*Dconjg(HIJ(I_current_term))
     .                      *DIMAG(RESCALE_1B(I_current_term))


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
                  EndDO
!$OMP END DO
!$OMP END PARALLEL
           ZRIJ=RIJ
      CALL ZGEMV('N',Nconf,TID,Zoner,VOUT,Nconf,sumvout1,1,ZERO,VIN,1)



c          call cpu_time(finish)
c          exec_time=finish-start
c          write(6,*)"IN 1Body Execution DNS time:",exec_time,L
c          write(6,*)"IN: Execution CNK time:",cnk_time
c            write(6,*)" Number of elements",L, " out of ",Nconf,Nconf-L
c        write(6,*)" From :",(Sh_m(i),i=1,7)
c        write(6,*)" -->  :",(Sh_p(i),i=1,7)
c       write(6,'(a20,F26.16,F26.16)')"Dir DNS element",zrho
c       write(6,'(a20,F26.16,F26.16)')"Cng DNS element",zrhocnjg
c================================================================
c              Do k=1,Nconf
c           IF(VOUT(k).ne.ZERO)
c     .      write(6,*)k,"IN",VIN(k),"OUT",VOUT(k)
c              Enddo
c               stop


c        VIN=VOUT
       DEALLOCATE(VOUT)
c       CI_Production_1b=.TRUE.
       CI_Production_1b=CI_Prd

c          call cpu_time(finish)
c          exec_time=finish-start
c          write(6,*)IPRC,TillN-FromN,"CPU Get 1b",exec_time

c           write(6,*)MYID, " 1B done "

        end subroutine GetCIJKL1body_27_OMP
