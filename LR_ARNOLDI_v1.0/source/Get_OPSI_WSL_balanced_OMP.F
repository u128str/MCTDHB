c============== This subroutine computes 2body part of the right hand side part of the equations for DPsi/Dt -
c============== usig WSL- based algorithm - one needs to get an access
c============== to all elements of the 2-body density density matrix:
c============== RhoAll(I.J.K.L) I,J,K,L=1..Morb
Corg       SUBROUTINE  Get_OPSI_WSL_based (PSI,OPSI,ZRIJKL,
Corg     .             NDX,NDY,NDZ,Morb,
Corg     .             delta,xlambda0,
Corg     .             Fij,F,weight,ort,NtVec,
Corg     .             FromN,TillN)

        SUBROUTINE  Get_OPSI_WSL_balanced_OMP (PSI,OPSI,FromN,TillN,
     &                                                   WSLTime)
c        use PASS_ARG
        USE SHARED_DIMS
        USE W_INTERPARTICLE
        USE rR_hW
        USE CI_All
        USE omp_lib

        IMPLICIT NONE
c=================== MPI ================================
       INCLUDE 'mpif.h'
       INTEGER ::  ierr,MYID,numprocs
       real*8 :: tm,tmb
       real*8 :: WSLTime
c=================== END MPI ================================
       COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb)  :: PSI
       COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb)  :: OPSI
!       COMPLEX*16, DIMENSION(NDX*NDY*NDZ)  :: temp1,temp2,WSL
       COMPLEX*16,ALLOCATABLE :: temp1(:),temp2(:),WSL(:)
c============================== Timings  =============================
       real*4   start,finish, exec_time , cnk_time       
c=====================================================================
       INTEGER :: K,S,L,Q,I,n1,Nadr,ii,n11,SL,P,III
c========================= ============================================
c       INTEGER, INTENT(IN) ::  FromN,TillN
       INTEGER  ::  FromN,TillN
c============================ Data with 2-body density elements =============================
        COMPLEX*16 :: RhoAll,Z,ZDOTC,DZNRM2,w1,w2,tw1,tw2
        EXTERNAL RhoAll,Nadr,ZDOTC,DZNRM2
c============================ RhoAll(K,L,Q,S,ZRIJKL) complex function
c============================ providing required element of density
c============================ matrix, assuming that all unique elements of R_ijkl are packed as (i<=j)<=(k<=l)
c=====================================================================

       allocate(temp1(NDX*NDY*NDZ),stat=ierr)
       allocate(temp2(NDX*NDY*NDZ),stat=ierr)
       allocate(WSL(NDX*NDY*NDZ),stat=ierr)
       if(ierr /= 0)write(*,*)"alloc error Get_OPSI_WSL_balanced_OMP"
c=====================================================================
!          write(6,*)" I am in WSL"
!            stop
        tm=MPI_WTIME(ierr)
c=====================================================================
!$OMP PARALLEL DO PRIVATE(ii,k) SHARED (OPSI,ndx,ndy,ndz,Morb)
                Do ii=1, ndx*ndy*ndz
                Do k=1,Morb
        OPSI(ii,k)=(0d0,0d0)
                 Enddo
                 Enddo
!$END OMP PARALLEL DO
c            OPSI=Zero

!          write(6,*)" I ma in WSL norm:",FromN,TillN
c     . , ZDOTC(NDX*NDY*NDZ,PSI(:,1),1,PSI(:,1),1)
         IF(ENRG_EVAL.eqv..TRUE.) THEN
         WIJKL=Zero
         ENDIF


c=====================================================================
c        SLloop: DO SL=1,Morb*(Morb+1)/2
        SLloop: DO SL=FromN,TillN
c================ Unpack cL cS from P
        P=TERM_INDEX_1B(SL)
        L= INT(P/100)
        S= P-L*100
        CALL  Get_WSL_omp(WSL,psi(:,S),psi(:,L),WSLTime)
           Qloop: DO Q=1,Morb
!$OMP PARALLEL DO
!$OMP& PRIVATE(ii)  
!$OMP& SHARED (temp1,temp2,PSI,WSL,ndx,ndy,ndz,Q)
                Do ii=1, ndx*ndy*ndz
                    temp1(ii)=WSL(ii)*PSI(ii,Q)
                    temp2(ii)=Conjg(WSL(ii))*PSI(ii,Q)
                 Enddo
!$END OMP PARALLEL DO
           Kloop:  DO K=1,Morb
           w1=dcmplx(0.d0,0.d0)
           w2=dcmplx(0.d0,0.d0)
!$OMP PARALLEL 
!$OMP& PRIVATE(II) 
!$OMP& SHARED(w1,w2,temp1,temp2,PSI,WSL,ndx,ndy,ndz,K)
!$OMP DO REDUCTION(+:w1,w2)
                Do ii=1, ndx*ndy*ndz
                   tw1=DCONJG(PSI(ii,K))*temp1(ii)
                   w1=w1+tw1
                   tw2=DCONJG(PSI(ii,K))*temp2(ii)
                   w2=w2+tw2
                Enddo
!$OMP END DO
!$OMP END PARALLEL 

                   IF(ENRG_EVAL.eqv..TRUE.) THEN
                   n1= ABS(Nadr(K,S,Q,L,Morb))
                   n11= Nadr(K,S,Q,L,Morb)
          IF(n11.ge.0)  THEN 
          WIJKL(n1)=WIJKL(n1)+w1
          ENDIF
                    ENDIF
                    Z=RhoAll(K,S,Q,L)
!$OMP PARALLEL DO
!$OMP& PRIVATE(ii)  
!$OMP& SHARED (Z,temp1,OPSI,ndx,ndy,ndz,K)
                Do ii=1, ndx*ndy*ndz
                OPSI(ii,K)=OPSI(ii,K)+Z*temp1(ii) 
                Enddo
!$END OMP PARALLEL DO
c            OPSI(:,K)=OPSI(:,K)+Z*temp1 
                     IF(L.eq.S) cycle Kloop
c============================ CYCLE is to awoid double counting in the
c============================ case of W11
                    IF(ENRG_EVAL.eqv..TRUE.) THEN
                    n1= ABS(Nadr(K,L,Q,S,Morb))
                    n11=    Nadr(K,L,Q,S,Morb)
          IF(n11.ge.0) THEN
          WIJKL(n1)=WIJKL(n1)+w2
          ENDIF
                     ENDIF
                    Z=RhoAll(K,L,Q,S)
!$OMP PARALLEL DO
!$OMP& PRIVATE(ii)  
!$OMP& SHARED (Z,temp2,OPSI,ndx,ndy,ndz,K)
                Do ii=1, ndx*ndy*ndz
                OPSI(ii,K)=OPSI(ii,K)+Z*temp2(ii) 
                Enddo
!$END OMP PARALLEL DO
c            OMP_OPSI(:,K,this_thread)=OMP_OPSI(:,K,this_thread)+Z*temp2
             END DO Kloop
           END DO Qloop
c=============== Reindexing - I need it for parallelization =========
         END DO SLloop
        tm=MPI_WTIME(ierr)
        timeCPU(5)=timeCPU(5)+MPi_WTIME(ierr)-tm
         deallocate(temp1,temp2,WSL)
      end subroutine  Get_OPSI_WSL_balanced_OMP

