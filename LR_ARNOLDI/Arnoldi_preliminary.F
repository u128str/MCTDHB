C=============== Many-Body-WF-Analysis 
       subroutine Arnoldi_preliminary
       USE    PASS_ARG 
       USE   SHARED_DIMS
       USE   rR_hW
       USE   CI_All
!       USE   PROP_MB
       USE   W_INTERPARTICLE
       USE   DVR_ALL
       USE   Parallel_Orb
!       USE   CORRELATIONFUNCTIONS
!       USE ANALYZER
       USE  CI_prod
!       USE   LINEAR_RESPONSE
       USE   LR_ARNOLDI_MOD
       USE   LR_Rapha
       implicit NONE
       INCLUDE 'mpif.h'
c==========================================================
c                    CI Variables
c=============================================================
       INTEGER ::  n,MAXIT,MacroIter_Max,MacroIter,ii,jj
       INTEGER ::  ici=1,ici_prt
       INTEGER ::  cI,cJ,cK,CL,P,j, ierr
c=================== F90 F95
        SAVE
        COMPLEX*16, ALLOCATABLE :: VIN(:),VOUT(:) 
c        COMPLEX*16, ALLOCATABLE :: Rho_ij(:,:),Rho_ijkl(:,:,:,:)
        real*8, ALLOCATABLE :: NatGemOcc(:)
c====================== For SIL
        COMPLEX*16 :: zrho=ZERO,Z,Z1,zerr,zerrci,zerr_new
        COMPLEX*16 :: Escale,W,E
c==============================================
      real*8  xnrm,Error_SIL,xerr,eerr,Tmax,xerr_td,xi,xf
      LOGICAL  CNV,SIL,CIC,Fexist,Get_INTRPL
      real*8   FKTR,CNK,DZASUM,DZNRM2,DNRM2,DDOT
c===================== t E S T ========================================
      integer :: M,dilation,icntr,Iorb,Jorb,icnt,myid,numprocs
       COMPLEX*16, DIMENSION(100,100)  :: dOerr
       COMPLEX*16, DIMENSION(100)  :: dOerr_diag
c=============================================================
        COMPLEX*16, ALLOCATABLE ::  OPSI(:,:)
        COMPLEX*16, ALLOCATABLE ::  Rho_PR(:),FTPSI(:,:)
        COMPLEX*16, ALLOCATABLE ::  PSI(:,:)
        COMPLEX*16, ALLOCATABLE ::  qq(:), fftqq(:,:)
        REAL*8, ALLOCATABLE ::  momX(:)
        COMPLEX*16, ALLOCATABLE ::  PSI_INTERPOLATED(:,:)
       Real*8 :: time 
       REAL*8 :: tau,tau_new,ABMError,TolError_Total,t0
       REAL*8 :: time_next,print_next,tau_in,print_step
       REAL*8 :: time_cic, time_psi,scl,scl_cic
       COMPLEX*16 :: ZDOTC
       INTEGER ::  DM,time_points_org,iii
       character*20 jtp
       character*10 path,pathCI,pathORB,pathORBK,pathg1RR,pathg1KK
       character*10 pathLR
       character*100 pathLRtime
c===========================================================
       INTEGER :: I,Nproc,istep
c==================== Interpolation ============================
        COMPLEX*16, DIMENSION(:), ALLOCATABLE :: invec,outvec
        REAL*8 :: newgrid(NDX)
c=============================================================
      external FKTR,CNK
      external Expect_fiXfi
      external DZASUM,DZNRM2,DNRM2,ZDOTC,DDOT,ZSCAL,Integrator_ABM
c=============================================================
c===========================================================
         write(6,*)" Master task for ",Npar," particles and ", Morb,
     $" orbitals is  C(N+M-1,N) long:", Nconf
       ALLOCATE(VOUT(Nconf))
       IF(ALLOCATED(VOUT).eqv..FALSE.) 
     . write(6,*) " MEM for VOUT  NOT is ok!", SIZE(VOUT)
       ALLOCATE(VIN(Nconf))
       IF(ALLOCATED(VIN).eqv..FALSE.) 
     . write(6,*) " MEM for VIN  NOT is ok!", SIZE(VIN)
c=========================================================================
         allocate(momX(NDX*NDY*NDZ),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in momX"
         allocate(qq(NDX*NDY*NDZ),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in qq"
         allocate(fftqq(NDX*NDY*NDZ,Morb),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in fftqq"
         allocate(Rho_PR(NDX*NDY*NDZ),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in Pho_PR"
         allocate(OPSI(NDX*NDY*NDZ,Morb),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in OPSI"
         allocate(PSI(NDX*NDY*NDZ,Morb),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in PSI"
         allocate(PSI_INTERPOLATED(NDX*NDY*NDZ,Morb),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in PSI_INTERPOLATED"
         allocate(FTPSI(NDX*NDY*NDZ,Morb),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in FTPSI"
         allocate(NatGemOcc(Morb*(Morb+1)/2),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in Rho_ij"
!         allocate(Rho_ij(Morb,Morb),stat=ierr)
!         if(ierr /= 0)write(*,*)"allocation error in Rho_ij"
!         allocate(Rho_ijkl(Morb,Morb,Morb,Morb),stat=ierr)
!         if(ierr /= 0)write(*,*)"allocation error in Rho_ijkl"
!         allocate(H_ij(Morb,Morb),stat=ierr)
!         if(ierr /= 0)write(*,*)"allocation error in H_ij"

!         allocate(W_ijkl(Morb,Morb,Morb,Morb),stat=ierr)
!         if(ierr /= 0)write(*,*)"allocation error in W_ijkl"
         allocate(invec(NDX*NDY*NDZ),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in invec"
         allocate(outvec(NDX*NDY*NDZ),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in outvec"
c======================== G U E S S ===================================
       IF(ORB_Binr) then
       CALL Open_binary !Needed for restarts
       CALL Read_Header(778,'scan')
       CALL Read_Header(777,'scan')
       ENDIF
c===================== At this point all the binary files are open and====
c===================== ready to read the PSI and CI and ==================
 
             path='DATA'
         call mk_dir(path)
              pathORB='DATA'//path_sep//'orb_R'
         call mk_dir(pathORB)
              pathCI='DATA'//path_sep//'CIcnf'
         call mk_dir(pathCI)
!        IF(PRINT_DATA.eqv..TRUE.) call print_rho(time,PSI,pathORB,Rho_PR)
c=========================== Cleaning some old files
!      open(1234, file=trim(pathoRB)//'/'//trim('DNS_Minkovskii.dat'),
!     . status='unknown', iostat=ierr)
!      if ( ierr .eq. 0) close(1234, status='delete')
c============================Initialization of the MPI
       CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,ierr)
       CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,ierr)
!       CALL MNGR_CI_Part_Julian(NUMPROCS)
       CALL MNGR_CI_Part(NUMPROCS)
       CALL MNGR_Orb_Part(NUMPROCS)

c==========================================================================
      time=0d0
      time=Time_Bgn
      write(6,*)"          PSI Data exist from time:",Time_Bgn
      write(6,*)"          User wants LR analysis at:",T_From
      time=MAX(T_From,Time_Bgn)
      time_LR=time

!       CI_SCF=0 ! CI_SCF=0 -> CI PART
!       call MPI_BCAST(CI_SCF,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!       call Share_H_W

      CALL READ_BINARY(778, time, 1000,VIN,PSI)
      CALL READ_BINARY(777, time, 1000,VIN,PSI)
      write(*,*) myid,(REAL(VIN(I)),I=1,Nconf)

       write(*,*) myid,"ZRIJ",(ZRIJ(I),I=1,Rdim) 

!AIS 2017
        call read_parse_V_W_Psi !parsing V_W_Psi_string.in file
         write(6,*) "File with user potential was parsed" 


       xlambda0=xlambda_0
       CALL Get_InterPart(time)

       CALL Get_h_W(PSI,time)
       CI_SCF=0 ! CI_SCF=0 -> CI PART
       call MPI_BCAST(CI_SCF,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       call Share_H_W


       SIL=.FALSE.
       call MPI_BCAST(SIL,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
       write(*,*) myid,'AFTER SIL BCAST'
       write(*,*) myid,"Inv",InvZRIJ 
       write(*,*) "LINE 155",CI_SCF 
       CALL Get_r_R(VIN,VOUT,time)
       write(*,*) "HERE LINE 168"
       CI_SCF=1 ! CI_SCF=1 -> SCF PART
       call MPI_BCAST(CI_SCF,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       call Share_r_R
       ABMError=TolError_Total/1.0d0
       ENRG_EVAL=.FALSE.
c       ENRG_EVAL=.TRUE.
       CALL FUNCSTR(time,PSI,OPSI) 
       ENRG_EVAL=.FALSE.
       JOB_TYPE=999
       call MPI_BCAST(JOB_TYPE,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!       SIL=.TrUE.
!       call MPI_BCAST(SIL,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)


       CI_SCF=2 ! CI_SCF=1 -> SCF PART
       call MPI_BCAST(CI_SCF,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

       write(*,*) myid,"Inv",InvZRIJ 
       write(*,*) myid,"ZRIJ",AllZRIJ 
       write(*,*) myid,"OP_x",Op_x(1,1) 

       write(*,*) "HERE LINE 168"
      CALL Close_binary 



!======= BROADCASTING PSI, VIN and time, call picked up in "Prop_Main.F"
       PSI_LR=PSI
       VIN_LR=VIN 
       time_LR=time

         CALL  MPI_BCAST(VIN_LR,Nconf,MPI_DOUBLE_COMPLEX,0
     .            ,MPI_COMM_WORLD,ierr)

         CALL  MPI_BCAST(PSI_LR,NDX*NDY*NDZ*Morb,MPI_DOUBLE_COMPLEX,0
     .            ,MPI_COMM_WORLD,ierr)

         CALL  MPI_BCAST(time_LR,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)

         CALL  MPI_BCAST(xlambda_0,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)

       write(*,*) myid,Zmu 
c=========================================================================
!      DEALLOCATE(VIN)
      DEALLOCATE(qq,fftqq)
      DEALLOCATE(Rho_PR)
      DEALLOCATE(OPSI)
      DEALLOCATE(PSI)
      deallocate(PSI_INTERPOLATED)
!      DEALLOCATE(Rho_ij)
!      DEALLOCATE(Rho_ijkl)
      DEALLOCATE(NatGemOcc)
      return
      end subroutine Arnoldi_preliminary
c=========================================================================
c=======================================================================
       subroutine Get_Full_Rij(Rho_JK)
       USE   PASS_ARG
       USE   SHARED_DIMS
       USE   rR_hW
       USE   CI_All
!       USE   PROP_MB
       USE   W_INTERPARTICLE
       USE   DVR_ALL
       USE   Parallel_Orb
       integer :: i,P,cK,cJ,icntr,Iorb,Jorb
       complex*16 :: rho_jk(Morb,Morb)
c=====================================================================
c      write(6,*) "IN=",MaxTrm1b
!         ZRIJ=ZRIG+(3.0,5.0)
        Do I=1,MaxTrm1b
        P=TERM_INDEX_1B(I)
c================ Unpack cI cJ cK cL from P
        cK= INT(P/100)
        cJ= P-cK*100
        Rho_JK(cJ,cK)=ZRIJ(I) 
        IF(cK.ne.cJ) Rho_JK(cK,cJ)=Conjg(ZRIJ(I))
!        write(6,*) cJ,cK,"NEW Rho",Rho_JK(cJ,cK)
        EnDdo
!        write(6,*) "NEW", Rho_JK
        
C========== BOTH ARE CORRECT !!!!!!!!!!!!!!
        icntr=1
        DO Iorb=1,Morb
        DO Jorb=Iorb,Morb

        Rho_JK(Iorb,Jorb)=ZRIJ(icntr)
        IF(Jorb.ne.Iorb) Rho_JK(Jorb,Iorb)=Conjg(ZRIJ(icntr)) ! Correct
        icntr=icntr+1
        EndDO
        EndDO
!        write(6,*) "OLD", Rho_JK
c        write(6,*) "OLD Rho",Rho_JK
!          pause
      end subroutine Get_Full_Rij
c=======================================================================
       subroutine Get_Full_Rijkl(Rho_ijkl)
       USE   PASS_ARG
       USE   SHARED_DIMS
       USE   rR_hW
       USE   CI_All
!       USE   PROP_MB
       USE   W_INTERPARTICLE
       USE   DVR_ALL
       USE   Parallel_Orb
       integer :: i,P,cK,cJ,cL,cI
       complex*16 :: rho_jk(Morb,Morb),check
       complex*16 :: rho_ijkl(Morb,Morb,Morb,Morb)
c=====================================================================
c      write(6,*) "IN=",MaxTrm2b
         Do I=1,MaxTrm2b
          P=TERM_INDEX_2B(I)
c================ Unpack cI cJ cK cL from P
          cL= INT(P/1000000)
          cK= INT((P-cL*1000000)/10000)
          cJ= INT((P-cL*1000000-cK*10000)/100)
          cI= P-cL*1000000-cK*10000-cJ*100
        Rho_ijkl(cI,cJ,cK,cL)=ZRIJKL(I) 
        Rho_ijkl(cI,cJ,cL,cK)=ZRIJKL(I) 
        Rho_ijkl(cJ,cI,cK,cL)=ZRIJKL(I) 
        Rho_ijkl(cJ,cI,cL,cK)=ZRIJKL(I) 
      IF(cI.ne.cK) then
       Rho_ijkl(cK,cL,cI,cJ)=Conjg(ZRIJKL(I))
       Rho_ijkl(cK,cL,cJ,cI)=Conjg(ZRIJKL(I))
       Rho_ijkl(cL,cK,cI,cJ)=Conjg(ZRIJKL(I))
       Rho_ijkl(cL,cK,cJ,cI)=Conjg(ZRIJKL(I))
      else
      IF(cJ.ne.cL) then
        Rho_ijkl(cK,cL,cI,cJ)=Conjg(ZRIJKL(I))
        Rho_ijkl(cK,cL,cJ,cI)=Conjg(ZRIJKL(I))
        Rho_ijkl(cL,cK,cI,cJ)=Conjg(ZRIJKL(I))
        Rho_ijkl(cL,cK,cJ,cI)=Conjg(ZRIJKL(I))
      endif
      endif
      EndDo

      end subroutine Get_Full_Rijkl

      subroutine Get_Full_Wijkl(W_ijkl)
       USE   PASS_ARG
       USE   SHARED_DIMS
       USE   rR_hW
       USE   CI_All
!       USE   PROP_MB
       USE   W_INTERPARTICLE
       USE   DVR_ALL
       USE   Parallel_Orb
       integer :: i,q,s,r,P,cK,cJ,cL,cI
       complex*16 :: W_ijkl(Morb,Morb,Morb,Morb)
c=====================================================================
         Do I=1,MaxTrm2b
          P=TERM_INDEX_2B(I)
c================ Unpack cI cJ cK cL from P
          cL= INT(P/1000000)
          cK= INT((P-cL*1000000)/10000)
          cJ= INT((P-cL*1000000-cK*10000)/100)
          cI= P-cL*1000000-cK*10000-cJ*100
        IF(mod(cL+cK+cJ+cI,2)==0) THEN
          W_ijkl(cI,cJ,cK,cL)=WIJKL(I) 
          W_ijkl(cI,cJ,cL,cK)=WIJKL(I) 
          W_ijkl(cJ,cI,cK,cL)=WIJKL(I) 
          W_ijkl(cJ,cI,cL,cK)=WIJKL(I) 
          IF(cI.ne.cK) then
             W_ijkl(cK,cL,cI,cJ)=Conjg(WIJKL(I))
             W_ijkl(cK,cL,cJ,cI)=Conjg(WIJKL(I))
             W_ijkl(cL,cK,cI,cJ)=Conjg(WIJKL(I))
             W_ijkl(cL,cK,cJ,cI)=Conjg(WIJKL(I))
          else
          IF(cJ.ne.cL) then
             W_ijkl(cK,cL,cI,cJ)=Conjg(WIJKL(I))
             W_ijkl(cK,cL,cJ,cI)=Conjg(WIJKL(I))
             W_ijkl(cL,cK,cI,cJ)=Conjg(WIJKL(I))
             W_ijkl(cL,cK,cJ,cI)=Conjg(WIJKL(I))
          endif
          endif
        ELSE
          W_ijkl(cI,cJ,cK,cL)=WIJKL(I)!*0.5d0 
          W_ijkl(cI,cJ,cL,cK)=WIJKL(I)!*0.5d0 
          W_ijkl(cJ,cI,cK,cL)=WIJKL(I)!*0.5d0  
          W_ijkl(cJ,cI,cL,cK)=WIJKL(I)!*0.5d0  
          IF(cI.ne.cK) then
             W_ijkl(cK,cL,cI,cJ)=Conjg(WIJKL(I))!*0.5d0  
             W_ijkl(cK,cL,cJ,cI)=Conjg(WIJKL(I))!*0.5d0  
             W_ijkl(cL,cK,cI,cJ)=Conjg(WIJKL(I))!*0.5d0  
             W_ijkl(cL,cK,cJ,cI)=Conjg(WIJKL(I))!*0.5d0  
!       WRITE(6,*) 'h11check', cI,cJ,cK,cL
          else
          IF(cJ.ne.cL) then
             W_ijkl(cK,cL,cI,cJ)=Conjg(WIJKL(I))!*0.5d0  
             W_ijkl(cK,cL,cJ,cI)=Conjg(WIJKL(I))!*0.5d0  
             W_ijkl(cL,cK,cI,cJ)=Conjg(WIJKL(I))!*0.5d0  
             W_ijkl(cL,cK,cJ,cI)=Conjg(WIJKL(I))!*0.5d0  
          endif
          endif
        END IF
      EndDo

!      DO i=1,Morb
!      DO q=1,Morb
!      DO r=1,Morb
!      DO s=1,Morb
!       WRITE(6,*) 'h11check', i,q,r,s, W_ijkl(i,q,r,s)
!      END DO
!      END DO
!      END DO
!      END DO

      end subroutine Get_Full_Wijkl

      subroutine Get_Full_Hij(H_ij)
       USE   PASS_ARG
       USE   SHARED_DIMS
       USE   rR_hW
       USE   CI_All
!       USE   PROP_MB
       USE   W_INTERPARTICLE
       USE   DVR_ALL
       USE   Parallel_Orb
       integer :: icntr, Iorb, Jorb
       complex*16 :: H_ij(Morb,Morb)
           icntr=1  
        DO Iorb=1,Morb
           DO Jorb=Iorb,Morb
              H_ij(Iorb,Jorb)=HIJ(icntr)       
              H_ij(Jorb,Iorb)=Conjg(H_ij(Iorb,Jorb))        
              icntr=icntr+1                                     
           EndDO   
        EndDO    
      end subroutine Get_Full_Hij

! To get input vector invec on a new grid actual goal is to use it to translate wavefunction from DVR to FFT grid  the PSI 
      subroutine interpolation1D(invec,outvec,
     &                        oldgrid,newgrid,NDold,NDnew)
       USE   PASS_ARG
       USE   SHARED_DIMS
       USE   rR_hW
       USE   CI_All
!       USE   PROP_MB
       USE   W_INTERPARTICLE
       USE   DVR_ALL
       USE   Parallel_Orb
       implicit NONE
       integer :: NDold,NDnew, i
       complex*16 :: invec(NDold),outvec(NDnew)
       real*8 :: Rinvec(NDold),Routvec(NDnew)
       real*8 :: oldgrid(NDold),newgrid(NDnew)
       real*8 :: diff(NDold,NDold),dfac,xi,xf,newweight(NDnew)
       integer enomax,enomin
       
      enomax=5
      enomin=1
      dfac=50.0d0 
      xi=Time_xint
      xf=Time_xfnl
      call getpos(xi,xf,NDnew,newgrid) !NEW FFT grid is constructed
      newweight = sqrt(newgrid(2)-newgrid(1))
      write(6,*) "INTERPOL DIM old  new",NDold,NDnew
!============== Real part interpolation
!Preparing Psi_i for interpolation - dvr weight must be taken into account
      do i=1,NDold
      Rinvec(i)=Real(invec(i))/weight_X(i)
      enddo
      call fdiff(oldgrid,Rinvec,NDold,diff)   !Diffeences between F(x_i)
      do i=1,NDnew
      call enoi(oldgrid,NDold,diff,newgrid(i),
     &enomax,enomin,dfac,Routvec(i))
      outvec(i)=Routvec(i)
      enddo
!============== Imaginary part interpolation
      do i=1,NDold
      Rinvec(i)=Imag(invec(i))/weight_X(i)
      enddo
      call fdiff(oldgrid,Rinvec,NDold,diff)   !Diffeences between F(x_i)
      do i=1,NDnew
      call enoi(oldgrid,NDold,diff,newgrid(i),
     &enomax,enomin,dfac,Routvec(i))
      outvec(i)=outvec(i) +(0d0,1.0)*Routvec(i)
      enddo
      outvec=outvec*newweight ! multiplication of the PSi on sqrt(weight)

!       OPEN(unit=117,File='VEC.out')
!          DO i=1,NDnew
!             WRITE(117,'(F16.10,2F16.8, F16.10,2F16.8)') 
!     & newgrid(i), outvec(i),oldgrid(i),invec(i)
!          END DO
!       CLOSE(117)
      write(6,*) "NORM  IN:",ABS(sum(conjg(invec(:))*invec(:)))
     & ,"INTERPOL NORM OUT:",ABS(sum(conjg(outvec(:))*outvec(:)))
!         stop
      end subroutine interpolation1D
   

