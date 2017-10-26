!>     Use natural geminals' analysis when ALL elements of the 
!>     reduced two-body denisty matrix \rho_ijkl 
!>     ZRIJKL(Morb*(Morb+1)/2) are computed!
        SUBROUTINE Diag_Geminal(time,NatGemOcc,jobtype)
        use PASS_ARG
        USE SHARED_DIMS
        USE rR_hW
        USE CI_ALL
        implicit NONE
       COMPLEX*16, ALLOCATABLE :: AllZRIJ1(:,:)
       COMPLEX*16, ALLOCATABLE :: WORK(:)

       REAL*8, ALLOCATABLE :: RWORK(:)
       REAL*8 NatGemOcc(Morb*(Morb+1)/2)

c=============================================================================
       INTEGER ::  ierr
c=============================================================================
          INTEGER :: I,J,K,icntr,Nadr
          INTEGER :: Iorb,Jorb,iflag,info
          INTEGER :: N,n1,n2,n3,LWORK,Mgem
          REAL*8 :: rfactor
          REAL*8 :: time
          INTEGER :: jobtype

          integer :: P,cK,cJ,cL,cI
c=============================================================================
          Mgem=Morb*(Morb+1)/2

        allocate(AllZRIJ1(Mgem,Mgem),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in Diag_Geminal"
        
        allocate(WORK(Mgem*Mgem),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in Diag_Geminal"

        allocate(RWORK(3*Mgem-2),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in Diag_Geminal"
        
        
c=============================================================================
c============= Here I reconstruct AllZRIJ(Mgem,Mgem)==========================

        AllZRIJ1=(0.d0,0d0)

        icntr=1
        DO Iorb=1,Mgem
        DO Jorb=Iorb,Mgem

        P=TERM_INDEX_2B(icntr)
c================ Unpack cI cJ cK cL from P
          cL= INT(P/1000000)
          cK= INT((P-cL*1000000)/10000)
          cJ= INT((P-cL*1000000-cK*10000)/100)
          cI= P-cL*1000000-cK*10000-cJ*100

!>      This subroutine returns the natural geminal occupation
!>      numbers for a given ZRIJKL reduced two-body denisty matrix
!>      strored as one-dimensional array. 

        rfactor=1.d0
        if (cI.ne.cJ) rfactor=rfactor*dsqrt(2.d0)
        if (cK.ne.cL) rfactor=rfactor*dsqrt(2.d0)

        AllZRIJ1(Iorb,Jorb)=ZRIJKL(icntr)*rfactor
c       write (6,*) P,cI,cJ,cK,cL,rfactor,Mgem
c       write (6,*) ZRIJKL(icntr)
c       write (6,*)

        icntr=icntr+1
        EndDO
        EndDO
c================== Natural Occupations and Natural Vectors ==================

        SELECT CASE (jobtype)
        CASE (1)
        CALL ZHEEV('N','U',Mgem,AllZRIJ1,Mgem,NatGemOcc,WORK,Mgem*Mgem,
     .                                  RWORK,INFO) ! Case 1: Evaluation
                                                    ! of eigenvalues
                                                    ! only.
        CASE (2)
        write (6,*) 'Diag_Geminal: Case (2) [eigenvectors+eigenvalues]' 
        write (6,*) '  still not supported.'
                    
        CASE DEFAULT
        write(6,*)" Something wrong in Diag_Geminal !!!!!!!!!"
                    END SELECT 




c      CALL ZHEEV('V','U',Mgem,AllZRIJ1,Mgem,NatGemOcc,WORK,Mgem*Mgem,
c    .                                    RWORK,INFO)
cccc       NatVec=AllZRIJ1
c       write (6,*) 'Sum of natural geminal occupations is:',
c    &  sum(NatGemOcc)
c    &  (NatGemOcc(i),i=1,Mgem)

         deallocate(AllZRIJ1,WORK,RWORK)
         END  SUBROUTINE Diag_Geminal
