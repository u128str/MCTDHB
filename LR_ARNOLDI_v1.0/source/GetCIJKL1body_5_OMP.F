        subroutine GetCIJKL1body_5_OMP(MYID,VIN)
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

       INTEGER, DIMENSION(5) :: nvecin,nvecout,Sh_m,Sh_p
       COMPLEX*16, DIMENSION(15) :: RIJ
        integer :: i1,j1
        integer :: i2,j2
        integer :: i3,j3
        integer :: i4,j4
        integer :: i5,j5
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
        Do i4=i3-1-Sh_m(4) ,Morb-4+Sh_m(5),-1
        nvecin(4)=i3-i4-1
        nvecin(5)=i4-1
        K=1+MCNK(i1-1,Morb-1) 
     $  +MCNK(i2-1,Morb-2)   
     $  +MCNK(i3-1,Morb-3)   
     $  +MCNK(i4-1,Morb-4)   
        j1=Npar+Morb-1 - (nvecin(1) - Sh_m(1) + Sh_p(1))
         j2= j1-1-(nvecin(2) - Sh_m(2) + Sh_p(2))
         j3= j2-1-(nvecin(3) - Sh_m(3) + Sh_p(3))
         j4= j3-1-(nvecin(4) - Sh_m(4) + Sh_p(4))
        J=1+MCNK(j1-1,Morb-1) 
     $  +MCNK(j2-1,Morb-2)   
     $  +MCNK(j3-1,Morb-3)   
     $  +MCNK(j4-1,Morb-4)   
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

        end subroutine GetCIJKL1body_5_OMP
