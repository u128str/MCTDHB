!C=============== Input: C^IC^JC_KC_L as a packed INTEGER number P=I+100*J+10000*K+1000000*L
!C===============        VECIN - incoming complex*16 vector
!C=============== Output: VECOUT - outcoming complex*16 vector: VECOUT=C^IC^JC_KC_L*VecIN
!C===============        < VECIN^| C^IC^JC_KC_L*VecIN>== rho_ijkl
!C=============== Lexicografical Ordering and numberging -> inverse combinadic (only zero's are adressed) 
!cORG       subroutine GetCIJKL2body_I(VIN,WIJKL,Escale,ZRHO,P,N,M,Nconf)
!c       subroutine  .GetCIJKL2body_Par(N,M,Nconf,VIN,ZWIJKL,ZRIJKL,TRM_REQ,FromN,TillN)
       subroutine  PrdCIJKL2body_M_OMP(MYID,VIN)
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
!c=================== MPI
!c       INCLUDE 'mpif.h'
!c       INTEGER ::  ierr,MYID,numprocs
!c=========================================================
       integer :: P,N,M 
       integer :: i,j,jcnjg,k,l

       integer :: cI,cJ,cK,cL,Nadr
       INTEGER ::  FromN,TillN,ii

!c=========================================================
        COMPLEX*16 :: VIN(:)
        COMPLEX*16, ALLOCATABLE  :: VOUT(:,:)
!        COMPLEX*16  :: VOUT(SIZE(VIN))
        COMPLEX*16, DIMENSION(52975) :: RIJKL !Till morb<=25 AIS 18JUL2014
        COMPLEX*16, DIMENSION(100) :: sumvout1=Zoner
        integer :: TID,this_thread
        INTEGER   , DIMENSION(25) :: nvecin,nvecout
      integer :: I_current_term,Nterms
      real*8 ::  xbprefac        
      real*4 ::  start,finish, exec_time , cnk_time       
!c=========================================================
        COMPLEX*16 :: ZDOTC
         real*8 ::  DZNRM2
!      write(6,*)" NOiN", SUM(Conjg(VIN(:))*VIN(:))
!      write(6,*)" Ind_CI_2b", Ind_CI_2b(1,1)
!      write(6,*)" Prefactors_2b", Prefactors_2b(1,1)
             
!c==========================================================================================================================
!        call cpu_time(start)
!c==========================================================================================================================
              IPRC=MYID+1
!c==========================================================================================================================
              FromN=CI_Proc_From(IPRC)
              TillN=CI_Proc_Till(IPRC)
!c=====================================================================
              TID = 1 !OMP_GET_MAX_THREADS()
!      write(6,*)" Prd with OpenMP number of threads:", TID

        ALLOCATE( VOUT(SIZE(VIN),TID) )
        IF(ALLOCATED(VOUT).eqv..FALSE.)& 
      write(6,*)" MEM for VOUT is NOT ok!",SIZE(VOUT)
        RIJKL=(0d0,0d0) !AIS 18JUL2014 this command at line 63 created problem in OMP simulations!! Now it is at the correct place!!!
           L=0
          vout=(0d0,0d0)!AIS 15AUG2014  - for the same as above reason - needed only for GNU :-)
!!$OMP PARALLEL  DEFAULT(PRIVATE) &
!!$OMP shared(VIN,VOUT,RIJKL,Nterms,TERM_INDEX_2B,&
!!$OMP TERM_REQ_2B,RESCALE_2B,WIJKL,TID,FromN,TillN,&
!!$OMP Ind_CI_2b,Prefactors_2b,Nconf)
!c=====================================================================
        this_thread= 1 !OMP_GET_THREAD_NUM()+1
!        VOUT(:,this_thread)=(0d0,0d0)
!!$OMP DO  REDUCTION(+:RIJKL)
       Iloop: DO I=FromN,TillN
       Do K=1,Nconf  
            RIJKL(TERM_REQ_2B(I))=RIJKL(TERM_REQ_2B(I))&
                        +Conjg(VIN(Ind_CI_2b(K,I-FromN+1)))&
                        *VIN(K)* SQRT(Prefactors_2b(K,I-FromN+1))
            VOUT(Ind_CI_2b(K,I-FromN+1),this_thread)=&
             VOUT(Ind_CI_2b(K,I-FromN+1),this_thread)&
                   +VIN(K)* SQRT(Prefactors_2b(K,I-FromN+1))&
                   *WIJKL(TERM_REQ_2B(I))&
                   *DREAL(RESCALE_2B(TERM_REQ_2B(I)))
            VOUT(K,this_thread)=VOUT(K,this_thread)&
                 +VIN(Ind_CI_2b(K,I-FromN+1))&
                 *SQRT(Prefactors_2b(K,I-FromN+1))&
                 *Conjg(WIJKL(TERM_REQ_2B(I)))*&
                  DIMAG(RESCALE_2B(TERM_REQ_2B(I)))
      EndDo
      End Do Iloop
!!$OMP END DO
!!$OMP END PARALLEL
!c================================================================
!      write(6,*)Rdim, "inter",(i,SUM(Conjg(VOUT(:,i))*VOUT(:,i)),i=1,TID)
           do i=1,Rdim*(Rdim+1)/2
           ZRIJKL(i)=RIJKL(i)
           enddo
      CALL ZGEMV('N',Nconf,TID,dcmplx(1.d0,0d0),&
                 VOUT,Nconf,sumvout1,1,dcmplx(0d0,0d0),VIN,1)
!      write(6,*)" NORM", SUM(Conjg(VIN(:))*VIN(:))
!        stop "PRD"
!c================================================================
      DEALLOCATE(VOUT)
!c================================================================
!      call cpu_time(finish)
!      exec_time=finish-start
!      write(6,*)IPRC,TillN-FromN,"CPU PRD 2b",exec_time



      end subroutine PrdCIJKL2body_M_OMP
