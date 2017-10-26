!C=============== Input: C^IC^JC_KC_L as a packed INTEGER number P=I+100*J+10000*K+1000000*L
!C===============        VECIN - incoming complex*16 vector
!C=============== Output: VECOUT - outcoming complex*16 vector: VECOUT=C^IC^JC_KC_L*VecIN
!C===============        < VECIN^| C^IC^JC_KC_L*VecIN>== rho_ijkl
!C=============== Lexicografical Ordering and numberging -> inverse combinadic (only zero's are adressed) 
       subroutine PrdCIJKL1body_M_OMP(MYID,VIN)
       USE SHARED_DIMS
       USE rR_hW
       USE CI_All
       USE Parallel_CI
       USE omp_lib
       USE CI_prod
       implicit NONE
!c=========================================================
       INTEGER ::  MYID,IPRC
!c=========================================================
       integer::  P,N,M 
       integer::  i,j,k,l,ii
       integer::  I_current_term,Nterms
       integer::  Sh_k,Sh_l

       integer::  cI,cJ,cK,cL
!c=========================================================
        COMPLEX*16 :: VIN(:)
        COMPLEX*16, ALLOCATABLE  :: VOUT(:,:)
        COMPLEX*16, DIMENSION(325) :: RIJ !till M<25 AIS 18JUL2014
        integer :: TID,this_thread
        COMPLEX*16, DIMENSION(100) :: sumvout1=Zoner

        COMPLEX*16 :: Escale,zdir,zinv
        INTEGER, DIMENSION(25) :: nvecin,nvecout
!c==============================================
       real*4::    start,finish,exec_time,cnk_time       
!c=====================================================================
       INTEGER ::  FromN,TillN
!c=====================================================================
!           call cpu_time(start)
               IPRC=MYID+1           
              FromN=CI_Proc_From(IPRC)
              TillN=CI_Proc_Till(IPRC)
!c=====================================================================
!c=====================================================================
         Nterms=MaxTrm1b
!c=====================================================================
!c======= MEMORY Allocation =====================================
!c       write(6,*)" MEM for  VIN 1B is ok!",SIZE(VIN)
 !             TID = OMP_GET_MAX_THREADS()
         TID=1
        ALLOCATE( VOUT(SIZE(VIN),TID) )
        IF(ALLOCATED(VOUT).eqv..FALSE.)& 
      write(6,*)" MEM for VOUT is NOT ok!",SIZE(VOUT)
!        VOUT=ZERO
        RIJ=(0d0,0d0) !AIS 18JUL2014 this command at line 63 created problem in OMP simulations!! Now it is at the correct place!!!

           vout=(0d0,0d0)
!c=====================================================================
!c=============== DO LOOP over the desired TRM_REQ
!c============== Nterms is the number of the required terms
           L=0
!!$OMP PARALLEL  DEFAULT(PRIVATE) &
!!$OMP shared(VIN,VOUT,RIJ,Nterms,TERM_INDEX_1B,&
!!$OMP TERM_REQ_1B,RESCALE_1B,HIJ,TID,FromN,TillN,&
!!$OMP Ind_CI_1b,Prefactors_1b,Nconf)
 !       this_thread= OMP_GET_THREAD_NUM()+1
                   this_thread=1
           !        VOUT(:,this_thread)=(0d0,0d0)
!!$OMP DO  REDUCTION(+:RIJ)
      DO I=FromN,TillN
!c================ Unpack cI cJ cK cL from P
        Do K=1,Nconf 
           RIJ(TERM_REQ_1B(I))=RIJ(TERM_REQ_1B(I))&
                     +Conjg(VIN(Ind_CI_1b(K,I-FromN+1)))*VIN(K)*&
                     DSQRT(Prefactors_1b(K,I-FromN+1))
!c=========== Conjugate C^_K C^_0 C_0 C_J   rho_0J_K0=rho^*_K0_0J
           VOUT(Ind_CI_1b(K,I-FromN+1),this_thread)=&
                         VOUT(Ind_CI_1b(K,I-FromN+1),this_thread)&
                         +VIN(K)* DSQRT(Prefactors_1b(K,I-FromN+1))&
                         *HIJ(TERM_REQ_1B(I))&
                         *DREAL(RESCALE_1B(TERM_REQ_1B(I)))
           VOUT(K,this_thread)=VOUT(K,this_thread)&
                              +VIN(Ind_CI_1b(K,I-FromN+1))*& 
                              DSQRT(Prefactors_1b(K,I-FromN+1))&
                              *DConjg(HIJ(TERM_REQ_1B(I)))&
                              *DIMAG(RESCALE_1B(TERM_REQ_1B(I)))

        EndDo
!c=========================================== 
 
      EndDo
!!$OMP END DO
!!$OMP END PARALLEL
      do i=1,Rdim
         ZRIJ(i)=RIJ(i)
      enddo
      CALL ZGEMV('N',Nconf,TID,dcmplx(1.d0,0d0),&
                 VOUT,Nconf,sumvout1,1,dcmplx(0d0,0d0),VIN,1)

 
      DEALLOCATE(VOUT)
!          call cpu_time(finish)
!          exec_time=finish-start
!          write(6,*)IPRC,TillN-FromN,"CPU PRD 1b",exec_time

      end subroutine PrdCIJKL1body_M_OMP
