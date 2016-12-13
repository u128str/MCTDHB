!====================================== MAIN get_LR subroutine =========================================
        SUBROUTINE Get_LinResp(time,PSI,VIN,pathLR)
          USE   LINEAR_RESPONSE
          USE CI_ALL
          USE Parallel_CI 
          USE W_INTERPARTICLE 
          USE DVR_ALL
          USE rR_hW
          USE SHARED_DIMS
        IMPLICIT NONE
        INCLUDE 'mpif.h'
!======LR variables=========================================
!         INTEGER :: dimL_orb,dimL,ND
        COMPLEX*16, DIMENSION(ND,Morb), INTENT(IN) :: PSI
        COMPLEX*16, DIMENSION(Nconf), INTENT(IN) :: VIN
        COMPLEX*16, DIMENSION(:), ALLOCATABLE :: VOUT
        COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: L_orb
        COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: L_full
        COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: h2
        COMPLEX*16, DIMENSION(:), ALLOCATABLE :: wo
        COMPLEX*16, DIMENSION(:), ALLOCATABLE :: norm_mat_orb
        COMPLEX*16, DIMENSION(:), ALLOCATABLE :: norm_mat_CI
        COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE :: u, v
        COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: Cu, Cv
        COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: COmat
        COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: OCmat
        REAL*8 :: timeh
        REAL*8   , DIMENSION(LR_maxsil-1)     :: EVALS
        REAL*8   , DIMENSION(LR_maxsil-1,LR_maxsil-1)  :: EV
        COMPLEX*16 :: CImat(Nconf,Nconf)
        REAL*8   , DIMENSION(Nconf)     :: EVALS_ex
        REAL*8   , DIMENSION(Nconf,Nconf)  :: EV_ex
        REAL*8   , DIMENSION(2*Nconf)     :: EVALS_d_ex
        REAL*8   , DIMENSION(2*Nconf,2*Nconf)  :: EV_d_ex
        INTEGER :: state_LR,k,l,ierr,jj
        character*10 pathLR
!================================================
        REAL*8 :: time
        COMPLEX*16 :: W,E
!=============================================================
!============= Dimensions ==============================================================================
           dimL_orb  =2*NDX*NDY*NDZ*Morb    
           dimL = dimL_orb + 2*Nconf    
           ND  =NDX*NDY*NDZ    
           LR_maxroot=MIN(LR_maxsil,LR_maxroot) !ADD BY STR to fix a problem in array ind
           numeig=MIN(LR_maxroot,dimL/2) !ADD BY STR to fix a problem in array ind
         if(numeig.ne.LR_maxroot) then
         write(*,*)"Number of required LR roots is changed from ",
     . LR_maxsil," to ",numeig
         endif
!============= ALLOCATION ==============================================================================
       ALLOCATE(VOUT(Nconf))
       ALLOCATE(L_orb(dimL_orb,dimL_orb))
       IF(ALLOCATED(L_orb).eqv..FALSE.)
     . write(6,*) " MEM for L_orb NOT is ok!", SIZE(L_orb)
       ALLOCATE(L_full(dimL,dimL))
       IF(ALLOCATED(L_full).eqv..FALSE.)
     . write(6,*) " MEM for L_full NOT is ok!", SIZE(L_full)
       ALLOCATE(h2(dimL_orb/2,dimL_orb/2))
       IF(ALLOCATED(h2).eqv..FALSE.)
     . write(6,*) " MEM for h2 NOT is ok!", SIZE(h2)
         allocate(wo(numeig),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in wo"
         allocate(norm_mat_orb(numeig),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in norm_mat_orb"
         allocate(norm_mat_CI(numeig),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in norm_mat_CI"
         allocate(u(ND,numeig,Morb),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in u"
         allocate(v(ND,numeig,Morb),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in v"
         allocate(Cu(Nconf,numeig),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in Cu"
         allocate(Cv(Nconf,numeig),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in Cv"
         allocate(COmat(2*Nconf,2*ND*Morb),stat=ierr)
         if(ierr /= 0)
     .   write(*,*)"allocation error in COmat",2*Nconf*2*ND*Morb,ierr
         allocate(OCmat(2*ND*Morb,2*Nconf),stat=ierr)
         if(ierr /= 0)
     .   write(*,*)"allocation error in OCmat",2*ND*Morb*2*Nconf,ierr
!c================0 STR added
        xlambda0=xlambda_0
!        CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,ierr)
!        CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,ierr)
!       CALL MNGR_CI_Part(numprocs) 
!        CALL MNGR_Orb_Part(NUMPROCS)
       CALL Get_InterPart(time)
       write(6,'(125a1)')("=",jj=1,118)
       write(6,'(125a1)')("=",jj=1,118)
       write(6,'(a30,F10.5,a26,i8,a28,i6)')
     . " LR analysis is at Time=",time,
     . " Size of LR-MCTDHB matrix =",dimL,
     . " number of roots to analyze =",numeig
       CALL Get_h_W(PSI,time)
       CALL Get_r_R(VIN,VOUT,time)
       CALL Get_h_W(PSI,time)
       CALL Get_r_R(VIN,VOUT,time)
       CALL Get_h_W(PSI,time)
       W=Zero
       E=SUM(AllHIJ*AllZRIJ)
       W=SUM(ZRIJKL*WIJKL*(DREAL(RESCALE_2B)))
     . +SUM(Conjg(ZRIJKL)*Conjg(WIJKL)*(DIMAG(RESCALE_2B)))
        Energy=Real(E+W)
       WRITE(6,'(a50,F26.16)') 
     . 'Energy in LR recomputed from W,h,rij,rijkl=', Real(E+W)
       WRITE(6,'(a20,e26.16,a20,e26.16)')
     . '      Im[E]=',Dimag(E),'      Im[W]=', Dimag(W)
       WRITE(6,*) 'LR xlambda0=', xlambda0
       write(6,'(125a1)')("=",jj=1,118)
       write(6,'(125a1)')("=",jj=1,118)
c=====================================================
c         CALL HPSI_LR(VIN,VOUT,Nconf)
         call Get_Full_Hij(H_IJ)
         call Get_Full_Wijkl(W_IJKL)
        CALL construct_matrix(PSI,VIN,L_full,Nconf,pathLR)
        write(6,*)"LR-MCTDHB matrix is sucessfully constructed DIM",dimL
        CALL diag_full_matrix(PSI,L_full,wo,u,v,Cu,Cv,Nconf,
     &                                 norm_mat_orb,norm_mat_CI,pathLR)
        write(6,*)"LR-MCTDHB matrix is sucessfully diagonalized"
        CALL evaluate(PSI,VIN,wo,u,v,Cu,Cv,Nconf,
     &                                 norm_mat_orb,norm_mat_CI,pathLR)
!         OPEN (32,FILE='CI.out')
!         DO i=1,Nconf
!               WRITE(32,'(I8,1X,2E16.8)') i,
!     &                          VIN(i)
!         END DO
!         CLOSE(32);
!
!==============DEALLOCATION=============================================================================
         DEALLOCATE(L_orb)
         DEALLOCATE(L_full)
         DEALLOCATE(h2)
         DEALLOCATE(wo)
         DEALLOCATE(u)
         DEALLOCATE(v)
         DEALLOCATE(Cu)
         DEALLOCATE(Cv)
         DEALLOCATE(COmat)
         DEALLOCATE(OCmat)
         write(6,*)"Okay LR is DONE for T",time
        END SUBROUTINE Get_LinResp 
