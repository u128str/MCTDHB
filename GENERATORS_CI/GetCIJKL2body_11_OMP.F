       subroutine GetCIJKL2body_11_OMP(MYID,VIN)
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

       INTEGER, DIMENSION(11) :: nvecin,nvecout,Sh_m,Sh_p
       COMPLEX*16, DIMENSION(2211) :: RIJKL
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
        Do i10=i9-1-Sh_m(10) ,Morb-10+Sh_m(11),-1
        nvecin(10)=i9-i10-1
        nvecin(11)=i10-1
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



        end subroutine GetCIJKL2body_11_OMP
