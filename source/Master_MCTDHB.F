C=============== Assistent MPI subroutive -assists in H|psi> multiplication
       subroutine Master_MCTDHB
       USE    PASS_ARG
       USE   SHARED_DIMS
       USE   rR_hW
       USE   CI_All
       Use OMP_lib
          USE   W_INTERPARTICLE
          USE   DVR_ALL
          USE   Parallel_Orb
       implicit NONE
c=================== MPI 
       INCLUDE 'mpif.h'
       INTEGER ::  ierr,MYID,numprocs
c==========================================================
c                    CI Variables
c=============================================================
       INTEGER ::  n,MAXIT,MacroIter_Max,MacroIter,ii,jj
       INTEGER ::  ici=1,ici_prt
       INTEGER ::  cI,cJ,cK,CL,P,j
c=================== F90 F95
!        SAVE
        COMPLEX*16, DIMENSION(:), ALLOCATABLE :: VIN 
        COMPLEX*16, DIMENSION(:), ALLOCATABLE :: VIN0
        COMPLEX*16, DIMENSION(:), ALLOCATABLE :: VOUT
c====================== For SIL
        COMPLEX*16 :: zrho=ZERO,Z,Z1,zerr,zerrci,zerr_new
        COMPLEX*16 :: Escale
c==============================================
      real*8, dimension(30) :: exp_val_U_op=0.d0
      real*8, dimension(12) :: time_rec(12)=0.0d0
      real*8 :: time_ci1=0.0d0,time_ci2=0d0,time_ci3=0d0
      real*8 :: time_psi1=0.0d0,time_psi2=0d0,time_psi3=0d0
c==============================================
      real*8  xnrm,Error_SIL,xerr,eerr,Tmax,xerr_td,ox,oy,oz
      LOGICAL  CNV,SIL,Fexist,Fexist1
      real*8   FKTR,CNK,DZASUM,DZNRM2,DNRM2,DDOT
c===================== t E S T ========================================
       integer :: M,K,ind
       COMPLEX*16, DIMENSION(100,100)  :: dOerr
       COMPLEX*16, DIMENSION(100)  :: dOerr_diag
c=============================================================
        COMPLEX*16, ALLOCATABLE ::  Rho_PR(:),Rho_PR1(:),Rho_PR2(:)
        COMPLEX*16, ALLOCATABLE ::  PSI(:,:)
        COMPLEX*16, ALLOCATABLE ::  PSI1(:,:)
        COMPLEX*16, ALLOCATABLE ::  PSI2(:,:)
        COMPLEX*16, ALLOCATABLE ::  OPSI(:,:)
        COMPLEX*16, ALLOCATABLE ::  AuxPsi(:)
        COMPLEX*16, ALLOCATABLE ::  PSI_IN(:,:)
       LOGICAL :: NewStart
       REAL*8 :: time, tau,tau_new,ABMError,TolError_Total,t0,tauprn
       REAL*8 :: time_next,print_next,tau_in,print_step
       REAL*8 x0in, gr
       REAL*8 :: Energy_old=0.0d0,exc_tm=0d0,exc_tm_new=0d0
       COMPLEX*16 :: ZDOTC
       INTEGER ::  DM
       character*20 jtp
       character*10 path,pathCI,pathORB
       character*255 cwd
c===========================================================
       INTEGER :: I,Nproc,istep
c================================================
       integer :: GetInd
!now in /TM/ global       integer :: print_level
c=============================================================
      external FKTR,CNK
      external Expect_fiXfi
      external DZASUM,DZNRM2,DNRM2,ZDOTC,DDOT,ZSCAL,Integrator_ABM
c=============================================================
!        print_level=2
#ifdef MKLFFT
      write(6,*) '################################################'
      write(6,*) '###### I USE THE FFT MKL LIBRARY ###############'
      write(6,*) '################################################'
#endif
#ifdef FFTW
      write(6,*) '################################################'
      write(6,*) '###### I USE THE FFT FFTW LIBRARY ##############'
      write(6,*) '################################################'
#endif
#ifdef CUDACPP
      write(6,*) '################################################'
      write(6,*) '###### I USE THE FFT CUDA  CPP LIBRARY #########'
      write(6,*) '################################################'
#endif
#ifdef CUDAPGI
      write(6,*) '################################################'
      write(6,*) '###### I USE THE FFT CUDA  PGI LIBRARY #########'
      write(6,*) '################################################'
#endif



          exc_tm=MPI_WTIME(ierr)
          exc_tm_new=MPI_WTIME(ierr)
c=============================================================
         call MPI_COMM_RANK(MPI_COMM_WORLD,MYID,ierr)
         call MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,ierr)
         write(6,*)" Process ",MYID," of ", numprocs
         CALL MNGR_CI_Part(NUMPROCS)
         CALL MNGR_Orb_Part(NUMPROCS)  
          timeCPU=0.0d0
c=============================================================
c       Nconf=NINT(CNK(Npar+Morb-1,Npar))
c===========================================================
         write(6,*)" Master task for ",Npar," particles and ", Morb,
     $" orbitals is  C(N+M-1,N) long:", Nconf
       ALLOCATE(VIN(Nconf))
       ALLOCATE(VIN0(Nconf))
       ALLOCATE(VOUT(Nconf))
       IF(ALLOCATED(VIN).eqv..FALSE.) 
     . write(6,*) " MEM for VIN  NOT is ok!", SIZE(VIN)
       IF(ALLOCATED(VIN0).eqv..FALSE.) 
     . write(6,*) " MEM for VIN0  NOT is ok!", SIZE(VIN0)
       IF(ALLOCATED(VOUT).eqv..FALSE.) 
     . write(6,*) " MEM for VOUT  is NOT ok!", SIZE(VOUT)
c=========================================================================
         allocate(Rho_PR(NDX*NDY*NDZ),stat=ierr)
         allocate(Rho_PR1(NDX*NDY*NDZ),stat=ierr)
         allocate(Rho_PR2(NDX*NDY*NDZ),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in Pho_PR"
         allocate(PSI(NDX*NDY*NDZ,Morb),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in PSI"
         allocate(PSI1(NDX*NDY*NDZ,Morb),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in PSI1"
         allocate(PSI2(NDX*NDY*NDZ,Morb),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in PSI2"
         allocate(OPSI(NDX*NDY*NDZ,Morb),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in OPSI"
         allocate(AuxPsi(NDX*NDY*NDZ*Morb*18),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in AuxPsi"
         allocate(PSI_IN(NDX*NDY*NDZ,Morb),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in PSI_IN"
c=========================================================================
            write(6,*)"Memory allocation is ok"
c=========================================================================
      xlambda0=xlambda_0
      write(6,*) " LAMBDA_0=",xlambda0
c======================== G U E S S ===================================
       IF(ORB_Binr) then
       CALL Open_binary !Needed for restarts
       CALL Read_Header(778,'scan')
       CALL Read_Header(777,'scan')
       CALL READ_BINARY(778, Binary_Start_Point_t, 1000,VIN,PSI)
       CALL READ_BINARY(777, Binary_Start_Point_t, 1000,VIN,PSI)
       CALL Close_binary !Needed for restarts
       ELSE
       call read_parse_V_W_Psi !parsing parse_V_W__PSI.in file
       CALL Guess_PSI(PSI)
       CALL Guess_CI(VIN)
       ENDIF
c================================================================
        CALL Open_binary ! Binary files are rewritten !!!!!!!!!!!!!!!!
        CALL Write_Header(777) !NAMELIST IS OVERWRITTEN to CIc binary file
        CALL Write_Header(778) !NAMELIST IS OVERWRITTEN to PSI binary file
c================================================================
c         Initialization of the arrays
c=========Initial TIME can be overwritten here ===========================
        time=Time_Bgn
c=========================================================================
       CALL Get_h_W(PSI,time)
        write(6,*) " Get_h_W 1 is ok",SUM(PSI(:,1)*PSI(:,1))
       CALL Get_r_R(VIN,VOUT,time)
       CALL Get_h_W(PSI,time)
c=================== Second order integration scheme ====================
c==================
c==================
c==================  TolError_Total - global error of the WAVE-FUNCTION
c================== means approxmatelly two order of magnitude better for total energy 
c================== TolError_Total = 10-6 for Wave-Function 
c==================                = 10-8 for Total Energy
       TolError_Total=Time_TolError_Total
c=================== Second order integration scheme ====================
       PreFac=JOB_PreFac
c========================================================================
      IF(PreFac.eq.(-1.0d0,0d0))write(jtp,*)"Job->Relax. Forward"
      IF(PreFac.eq.(+1.0d0,0d0))write(jtp,*)"Job->Relax.Backward"
      IF(PreFac.eq.(0d0,-1.0d0))write(jtp,*)"Job->Prop.  Forward"
      IF(PreFac.eq.(0d0,+1.0d0))write(jtp,*)"Job->Prop. Backward"
c      IF(PreFac.eq.(+1.0d0,0d0)) write(jtp,*)"JOB-> Relaxation Backward"
c       IF(PreFac.eq.(-1.0d0,0d0).or.PreFac.eq.(1.0d0,0d0)) then
c        write(jtp,*)"JOB->   Relaxation"
c       ELSE
c        write(jtp,*)"JOB->  Integration" 
c       ENDIF
c=================== Second order integration scheme ====================
        tau=Time_tau
c========================================================================
         Tmax=600.0d0
         Tmax=Time_Max
         Tmax=Time_Fnl
c========================================================================
       tau_in=tau
       print_step=0.1d0
       print_step=Time_print_step
       ici_prt=1 ! prints 10 times less often then psi's
       ici_prt=Time_ici_prt ! prints 10 times less often then psi's
       print_next=print_step+time
c========================================================================
c========================================================================
      call Write_PSI_string(time,-1,Energy,PSI)
      call Write_CIc_string(time,-1,Energy,VIN)
      Call flush(778)
      Call flush(777)
      write(6,*)" Initial stuff written to bin is",time,-1,Energy
             path='media'
        call mk_dir(path)
             path='DATA'
        call mk_dir(path)
             pathORB='DATA'//path_sep//'orb_R'          
        call mk_dir(pathORB)
             pathCI='DATA'//path_sep//'CIcnf'             
        call mk_dir(pathCI)
       IF(PRINT_DATA.eqv..TRUE.) call print_rho(time,PSI,pathORB,Rho_PR)
!       IF(PRINT_DATA.eqv..TRUE.) call print_rho3d(time,PSI)
       IF(PRINT_DATA.eqv..TRUE.) call print_cic(time,VIN,pathCI)
c      call Open_binary()
c========================================================================
          write(6,'(a66)')(banner(i)//CHALL,i=1,44)
          write(6,'(a66)')(banner0(i)//CHALL,i=1,7)
       write(6,'(100a1)')("=",jj=1,100)
       write(6,'(a40,a20,a40)')
     . "========================================",jtp,
     . " started ==============================" 
       write(6,'(100a1)')("=",jj=1,100)
c========================================================================
c SK + STR 08/01/15
        open(unit=9,file='OP_PR.out',form='formatted')
        call wr_info(9)
       write(9,'(23(10X,a9,10X))')"time  ",
     .  "  <E>/N  ", "  <T>/N  ",
     .  "  <V>/N  ", "  <W>/N  ",
     .  "  <x>/N  ", "  <y>/N  ", "  <z>/N  ",
     .  "<x*x>/N  ", "<y*y>/N  ", " <z*z>/N ",
     .  " Var(x)/N", "Var(y)/N ", "Var(z)/N ",
     .  " <px>/N  ", "  <py>/N ", " <pz>/N  ",
     .  "<px*px>/N", "<py*py>/N", "<pz*pz>/N",
     .  "Var(px)/N", "Var(py)/N", "Var(pz)/N"


        open(unit=10,file='NO_PR.out',form='formatted')
       call wr_info(10)
          timeCPU=0.0d0
c========================================================================
                    MacroIter_Max=400000
c                    MacroIter_Max=2000
          exc_tm=MPI_WTIME(ierr)-exc_tm
         write(6,*)" Initialization time took",exc_tm,"[sec]"
         Do istep=1,MacroIter_Max
c========================================================================
         time_rec(1)=MPI_WTIME(ierr) 
c========================================================================

        tauprn=tau !STR March 2015 to fix printout
        IF((Abs(1.0d0-Abs(Dreal(PreFac)))).le.1.0d-8) tauprn=tau/2
        IF(print_level.ge.2) write(6,'(100a1)')("=",jj=1,100)
        IF(print_level.ge.1) then
       write(6,'(a20,a15,I8,a12,F14.6,a3,F14.6,a3,F14.6,a3,F14.6,a3)
     .')jtp,
     ."Iteration:",istep,"Time: [",Time_Bgn,"->",time," + ",tauprn
     .,"->",Time_Fnl," ]"
         ENDIF !PRINT_LEVEL
       CALL Get_h_W(PSI,time)
         IF(print_level.ge.2) then
       write(6,'(a35,F26.16,a9,i10,a10,i3,a1,i10,a1,i10,a1)')
     ."Input orbital energy E(t+0):",Energy,"( CI Dim:",Nconf,
     .")(ORB Dim:",Morb,"*",NDX*NDY*NDZ,"=",Morb*NDX*NDY*NDZ,")"
          ENDIF !PRINT_LEVEL
c========================================================================
         time_rec(2)=MPI_WTIME(ierr)
c========================================================================
        VIN0=VIN
        PSI_IN=PSI !needed for step rejection
c========================================================================
c       write(6,*)"====================================================="
c       write(6,'(a10,10F15.10)')"IN CI:",(Abs(VIN(i)),i=1,9)
c       write(6,*)"====================================================="
c========================================================================
c========================================================================
c            I    CI P_a_r_t |                  0->tau/2
c========================================================================
        t0=tau/2
c=========================================================================
       CI_SCF=0 ! CI_SCF=0 -> CI PART
       call MPI_BCAST(CI_SCF,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       call Share_H_W
c        CALL MTRX_Diag
c=========================================================================
       Error_SIL=TolError_Total/1.0d0
       MAXIT= 20
       CALL SIL_PROPG(t0,VIN,VOUT,Error_SIL,MAXIT)
       If (abs(t0-tau/2).gt.1.d-8.and.ABS(DIMAG(PreFac)).ge.0.99d0)then
       write (6,*) "Rejection at 1ST SIL", t0,tau/2 
       zerr=1.d0
       time=time+tau/2
       tau_new=0.8d0*tau
       goto 8765
       End If
       VIN=VOUT
c========================================================================
         time_rec(3)=MPI_WTIME(ierr)
c========================================================================
c=========================================================================
       IF(PreFac.eq.(-1.0d0,0d0).or.PreFac.eq.(1.0d0,0d0)) VIN0=VOUT
c========================================================================
       CALL Get_r_R(VIN0,VOUT,time)
        PSI1=PSI
c=========================================================================
         time_rec(4)=MPI_WTIME(ierr)
c========================================================================
        t0=tau/2
        IF(MORB.eq.1) t0=tau
c=========================================================================
            OPSI =ZERO
            AuxPsi=ZERO
            NewStart=.True.
c========================================================================
       CI_SCF=1 ! CI_SCF=1 -> SCF PART
       call MPI_BCAST(CI_SCF,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       call Share_r_R
       ABMError=TolError_Total/1.0d0
       ENRG_EVAL=.FALSE.
c       ENRG_EVAL=.TRUE.
       CALL FUNCSTR(time,PSI1,OPSI) 
c          print *, "Before ABM"
       ENRG_EVAL=.FALSE.
       call Integrator_ABM (PSI1,OPSI,AuxPsi,time,t0,NewStart,ABMError)!  PSI1 orbitals at tau/2 with rij(0) Rijkl(0)
       JOB_TYPE=999
       call MPI_BCAST(JOB_TYPE,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c       write(6,*)"======= I ORB is done "
c           pause
c========================================================================
         time_rec(5)=MPI_WTIME(ierr)
c=========================================================================
              IF(MORB.eq.1) THEN
            CALL Get_h_W(PSI,time+t0)   !compute energy at the end of timne step 25AUG 2014 
              PSI=PSI1
              VOUT=VIN0
              tau_new=tau
              zerr=TolError_Total
              zerrci=Zero
              time=time+tau/2 
              goto 1234
              ENDIF
c=========================================================================
c=========================================================================
              IF(PreFac.eq.(-1.0d0,0d0).or.PreFac.eq.(1.0d0,0d0)) THEN
              PSI=PSI1
              VOUT=VIN0
              tau_new=tau
         IF(print_level.ge.3) then
         write(6,'(3(a10,F20.16 ))')
     ."Time:",time," New tau:",tau," Old tau:",tau_new
         write(6,*)" CI norm:",DZNRM2(Nconf,VIN0,1)
          ENDIF !PRINT_LEVEL
c               pause
              goto 1234
              ENDIF
c=========================================================================
         time_rec(6)=MPI_WTIME(ierr)
c========================================================================
       CALL Get_r_R(VIN,VOUT,time+tau/2)                              ! rij(tau/2) and Rijkl(tau/2)
            AllZRIJ0=AllZRIJ
c=========================================================================
         time_rec(7)=MPI_WTIME(ierr)
c=========================================================================
c       write(6,'(a25,F26.16)')"In CI  energy E(t+tau/2):",Energy
c       write(6,*)"================================Next======= After SIL"
c     .,time+tau/2
c       write(6,'(a10,10F15.10)')"New NO's:",(Nocc(i),i=1,Morb)
c       write(6,*)"====================================================="
c=========================================================================
c=========================================================================
c========================================================================
c             II and III  Orbital  P_a_r_t             0->tau/2 Needed for Error evaluation  and for propagation !!!!!!
c========================================================================
        t0=tau/2
c=========================================================================
            OPSI =ZERO
            AuxPsi=ZERO
            NewStart=.True.
c========================================================================
c                         II  Orbital
c       write(6,'(a10,10F15.10)')"II PSI uses H's:",(hij(i),i=1,Morb+1)
c========================================================================
       CI_SCF=1 ! CI_SCF=1 -> SCF PART
       call MPI_BCAST(CI_SCF,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       call Share_r_R
       ABMError=TolError_Total/1.0d0
       ENRG_EVAL=.FALSE.
c       ENRG_EVAL=.TRUE.
       CALL FUNCSTR(time,Psi,OPSI) 
       ENRG_EVAL=.FALSE.
c        pause
       call Integrator_ABM (PSI,OPSI,AuxPsi,time,t0,NewStart,ABMError) !  PSI2 orbitals at  tau/2   with rij(tau/2) Rijkl(tau/2)
c       write(6,*)"======= II ORB is done "
       PSI2=PSI
c=========================================================================
         time_rec(8)=MPI_WTIME(ierr)
c========================================================================
c                         III  Orbital
c       write(6,'(a10,10F15.10)')"III PSI uses H's:",(hij(i),i=1,Morb+1)
       time=time+tau/2
c========================================================================
       NewStart=.False.
       ENRG_EVAL=.FALSE.
       CALL FUNCSTR(time,Psi,OPsi) 
       ENRG_EVAL=.FALSE.
       call Integrator_ABM (PSI,OPSI,AuxPsi,time,t0,NewStart,ABMError) !  PSI orbitals at   tau   with rij(tau/2) Rijkl(tau/2)
       JOB_TYPE=999
       call MPI_BCAST(JOB_TYPE,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c       write(6,*)"======= III ORB is done "
c       pause
c========================================================================
       CALL Get_h_W(PSI,time+tau/2)                                    !  h(tau) Wijkl(tau)
c       write(6,'(a25,F26.16)')"Orbital energy E(t+tau/2):",Energy
c========================================================================
c            II and III    CI P_a_r_t |                  tau/2->tau and tau/2->0
c========================================================================
        t0=tau/2
c=========================================================================
       CI_SCF=0 ! CI_SCF=0 -> CI PART
       call MPI_BCAST(CI_SCF,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       call Share_H_W
c========================================================================
       time_rec(9)=MPI_WTIME(ierr)
c=========================================================================
       Error_SIL=TolError_Total/1.0d0
       MAXIT= 20
       CALL SIL_PROPG(t0,VIN,VOUT,Error_SIL,MAXIT)
       If (abs(t0-tau/2).gt.1.d-8) then
       write (6,*) "Rejection at 2nd SIL", t0,tau/2
       zerr=1.d0
       tau_new=0.8d0*tau
       goto 8765
       End If

c       write(6,*)"======= II and III CI are done "
c================= VIN  has CI Vout(time-tau/2) -> V(0) needed for error estimations 
c================= VOUT has CI Vout(time+tau) the results of propagation 
c=========================================================================
c=========================================================================
c=========================================================================
c================= Error evaluations!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c=========================================================================
c          PSI1 PSI2 VIN0, VIN
c=========================================================================
        VIN=VIN0-VIN
        zerr=ZDOTC(Nconf,VIN,1,VIN,1)
        zerrci=zerr
        xerr=ABS(zerr)
        xerr=0d0
        PSI2=PSI2-PSI1
        DM=NDX*NDY*NDZ
         dOerr=Zero
          zerr_new=Zero
        Do ii=1,Morb
           Do jj=1,Morb

              zerr_new=zerr_new 
     .                 +ZDOTC(DM,PSI2(:,II),1,PSI2(:,JJ),1)*
     .                 AllZRIJ(jj,ii)
           EndDo
        EndDo

c=============================================================================
          xerr=ABS(zerr_new)
         zerr=(ABS(zerr_new)+ABS(zerrci))*Zoner
c!================ New Error coming from TD potential NOT EXACT!!!!!  1-Exp(dE*dT)
           xerr_td=1.0d0-Exp((ABS(Energy-Energy_old))*ABS(tau))
         tau_new=TD_Err_scale*tau*(TolError_Total/REAL(ABS(zerr)))**0.25d0
8765      IF(ABS(zerr).gt.2*TolError_Total) THEN
          write(6,*) "!!!!! STEp IS REJECTED !!!!!!!",
     . ABS(zerr),">",TolError_Total
          VIN=VIN0
          PSI=PSI_IN
          time=time-tau/2
          CALL Get_h_W(PSI,time)
          CALL Get_r_R(VIN,VOUT,time)
          tau=0.8d0*tau_new
          cycle
          ENDIF
c=========================================================================
       time_rec(10)=MPI_WTIME(ierr)
c=========================================================================
1234    continue
c=========================================================================
       time_rec(11)=MPI_WTIME(ierr)
c=========================================================================
       VIN=VOUT
       time=time+tau/2
       CALL Get_r_R(VIN,VOUT,time)                                     ! rij(tau) and Rijkl(tau)
c=========================================================================
       time_rec(12)=MPI_WTIME(ierr)
c=========================================================================
         IF(print_level.ge.1) then
       write(6,'(a35,F26.16,a9,i10,a12,G17.10,a17,I1)')
     ."OUT  CI  energy E(t+tau):",Energy,
     ."   N = ",Npar," l0*(N-1)=",xlambda0*(Npar-1)
     .," kind of W(x-x'):",Wxx_TYPE
         ENDIF !PRINT_LEVEL

       IF(PreFac.eq.(-1.0d0,0d0).or.PreFac.eq.(1.0d0,0d0))
     .zerr=Energy-Energy_old
       eerr=Energy-Energy_old
         IF(print_level.ge.2) then
       write(6,'(2(a35,F26.16))')"Delta E     : +/-",Energy-Energy_old
     .                          ,"Error due to  dE:",xerr_td
         ENDIF !PRINT_LEVEL
         IF(print_level.ge.2) then
       write(6,'(a35,F26.16)')"Error E     : +/-",TolError_Total
         ENDIF !PRINT_LEVEL
       Energy_old=Energy

         IF(print_level.ge.3) then
       write(6,'(a35,2F26.16,a5)')"CI  Error   : +/-",ABS(zerrci)
     .,100*ABS(zerrci)/TolError_Total,"%"
       write(6,'(a35,2F26.16,a5)')"ORB Error   : +/-",  xerr
     .,100*xerr/TolError_Total,"%"
         ENDIF !PRINT_LEVEL

         IF(print_level.ge.3) then
       write(6,'(a10,a70)')" ",
     ."========= Natural Occupation numbers (only first 8)==========="
         ENDIF !PRINT_LEVEL
         IF(print_level.ge.2) then
      write(6,'(a16,8(G14.7,a3))')
     . "New NO's: ",(Nocc(Morb-i+1)," | ",i=1,Min(Morb,8))
         ENDIF !PRINT_LEVEL
         IF(print_level.ge.3) then
       write(6,'(a10,a70)')" ",
     ."=====================                            ============="
         ENDIF !PRINT_LEVEL
c      IF(ABS(tau_new).ge.0.05d0) tau_new=0.05d0*tau_new/ABS(tau_new)
      IF(ABS(tau_new).ge.Time_intgr_max_step) 
     .  tau_new=Time_intgr_max_step*tau_new/ABS(tau_new)

         IF(print_level.ge.3) then
       write(6,'(F12.6,a30,F12.6,a5,F12.6)')
     . time,"====time print next=====",print_next," tau=",tau
         ENDIF !PRINT_LEVEL

      if(ABS(time-print_next).le.1.0d-8) then
      print_next=print_next+print_step
      IF(PRINT_DATA.eqv..TRUE.) call print_rho(time,PSI,pathORB,Rho_PR)
!       IF(PRINT_DATA.eqv..TRUE.) call print_rho3d(time,PSI)
      call Write_PSI_string(time,istep,Energy,PSI)
      Call flush(778)
      if(ici.eq.ici_prt) then
      IF(PRINT_DATA.eqv..TRUE.) call print_cic(time,VIN,pathCI)
      call Write_CIc_string(time,istep,Energy,VIN)
      Call flush(777)
       ici=0
      endif
       ici=ici+1
      Endif

      tau=tau_new
      time_next=time+tau
      if(time_next-print_next.gt.1.0d-8) then
      tau=print_next-time
      endif        
c============== Case integration till desired time Tmax
      if(time_next-Tmax.gt.1.0d-8) then
      tau=Tmax-time
      write(6,*)"!!!====== Last Point arroaching Tmax=",Tmax
      write(6,*)"!!!====== Time-Step is reduced to   =",tau
      endif        



       IF(PreFac.eq.(-1.0d0,0d0).or.PreFac.eq.(1.0d0,0d0)) tau=tau_new

       write(10,2222)time,"  ",(Nocc(i)/Npar,"  ",i=1,Morb),Energy
     .,"  ",TolError_Total,"  ",eerr/TolError_Total
     .,"  ",ABS(zerr),"  ",ABS(zerrci)
     .,"  ",ABS(zerr)-ABS(zerrci),"  "
c=======================================
c SK+STR 08/01/15
        call User_OP(time,PSI,exp_val_U_op)
       write(9,2226)time,"  ",(exp_val_U_op(i),"  ",i=1,22)


2222  format((78(G26.16,a3)))
2226  format((78(E26.16,a3)))
        call flush(10) 
         IF(print_level.ge.3) then
       write(6,'(a75,F16.10,a20,F16.10)')"Time:",time
     .,"            New tau=",tau
         IF(ABS(zerr).ge.epsilon(0.d0))THEN
       write(6,'(a39,E26.16,2(a12,E12.5))')
     ."Overall ||PSI|| error  of this step is=",REAL(ABS(zerr))
     .," Err  CI[%]",100*ABS(zerrci)/TolError_Total
     .," Err ORB[%]",100*(ABS(zerr)-ABS(zerrci))/TolError_Total
          ENDIF
          ENDIF !PRINT_LEVEL 
c=========================================================================
       time_ci1=time_rec(3)-time_rec(2)
       time_ci2=time_rec(10)-time_rec(9)
       time_ci3=time_rec(4)-time_rec(3)
     .         +time_rec(7)-time_rec(6)
     .         +time_rec(12)-time_rec(11)
         IF(print_level.ge.3) then
       write(6,'(a15,F16.10,3(a10,F16.10))')
     ."CI  Time:",time_ci1+time_ci2+time_ci3,
     ."  CI   I:",time_ci1,
     ."  CI  II:",time_ci2,
     ." oth. CI:",time_ci3
          ENDIF !PRINT_LEVEL 
c=========================================================================
       time_psi1=time_rec(5)-time_rec(4)
       time_psi2=time_rec(8)-time_rec(7)
       time_psi3=time_rec(9)-time_rec(8)

         IF(print_level.ge.3) then
       write(6,'(a15,F16.10,4(a10,F16.10))')
     ."ORB Time:",time_psi1+time_psi2+time_psi3,
     ."  ORB   I:",time_psi1,
     ."  ORB  II:",time_psi2,
     ."  ORB III:",time_psi3,
     ." oth. ORB:",time_rec(2)-time_rec(1)
          ENDIF !PRINT_LEVEL 
c=========================================================================
       exc_tm=exc_tm+time_rec(12)-time_rec(1)
       write(6,'(a15,F16.10,a20,F16.10)')
     ."Itr Time:",time_rec(12)-time_rec(1)," execution time:",exc_tm
c       write(6,'(a15,F16.10,a20,F16.10)')
c     ."CPU Kin Time:",timeCPU(1)," CPU Non-lin Time:",timeCPU(2)
c       write(6,'(a15,F16.10,a20,F16.10)')
c     ."CPUKin FFT Time:",timeCPU(5)," CPU ZCPTime:",timeCPU(6)
c       write(6,'(a15,F16.10,a20,F16.10)')
c     ."CPU FFT Time:",timeCPU(3)," CPU FFT-DO Time:",timeCPU(4)

         IF(print_level.ge.3) then
       write(6,'(a17,F16.10)')"CPU NON-LIN Time:",timeCPU(7)
c       write(6,'(a17,F16.10)')"CPU FFT     Time:",timeCPU(5)
c       write(6,'(a15,F16.10)')"CPU V       Time:",timeCPU(6)
       write(6,'(a17,F16.10)')"CPU+KIN+V   Time:",timeCPU(1)
       write(6,'(a17,F16.10)')"CPU Proj    Time:",timeCPU(2)
       write(6,'(a17,F16.10)')"CPU Func    Time:",timeCPU(3)
          ENDIF !PRINT_LEVEL 
c=========================== STR Feb 2015 stop-file
        inquire(file='stop_MCTDHB',exist=Fexist)
        if (Fexist) then
          write(6,*)"STOP was required see ./stop_MCTDHB.done"
        inquire(file='stop_MCTDHB.done',exist=Fexist1)
             if (Fexist1) then
        open(121, file='stop_MCTDHB.done', 
     .   status="old", position="append", action="write")
        write(121,*)"Stop at Iteration",istep
        call wr_info(121) 
        close(121)
             else
         open(121, file='stop_MCTDHB.done',status="new", action="write")
        write(121,*)"Stop at Iteration",istep
        call wr_info(121) 
        close(121)
             endif
      open(1234,file='stop_MCTDHB',status='unknown', iostat=ierr)
        close(1234, status='delete')
          goto 1235
         end if
c=========================================================================
      if(abs(time-Tmax).le.0.0000000001d0) goto 1235
      if(time.ge.Tmax) goto 1235
       ENDDO
c=================== END Second order integration scheme ====================
1235    continue
c             Do i=1,Morb
c                  Rho_PR= psi(:,i)
c                  call Get_d_PSI(Rho_PR,1)
c                  call Get_d_PSI(Rho_PR,1)
c                  Z1=SUM(conjg(PSI(:,i))*(-Rho_PR/2))
c          write(6,*)"Txx (i=",i,")= ",Z1 
c                  Rho_PR= psi(:,i)
c             IF(Time_DVRMETHODX.ne.4) THEN
c                  CALL Get_Op_PSI(Rho_PR,Op_X,NDX,1)
c             ELSE
c#ifdef SFX1D              
c             CALL Get_Op_PSI_ SFX1D (Rho_PR,1)  !FFT 1D in X direction
c#endif
c             ENDIF
c        write(6,*)"RTxx(i=",i,")= ",SUM(conjg(PSI(:,i))*(Rho_PR))
c        write(6,*)"h_ii (i=",i,")= ",Z1
c     . +SUM(conjg(PSI(:,i))*(VTRAP_EXT*PSI(:,i)))
c
c                  Rho_PR= psi(:,i)
c                  call Get_d_PSI(Rho_PR,2)
c                  call Get_d_PSI(Rho_PR,2)
c                  Z1=Z1+SUM(conjg(PSI(:,i))*(-Rho_PR/2))
c          write(6,*)"Tyy (i=",i,")= ",SUM(conjg(PSI(:,i))*(-Rho_PR/2))
c              Rho_PR= psi(:,i)
c             IF(Time_DVRMETHODY.ne.4) THEN
c                  CALL Get_Op_PSI(Rho_PR,Op_Y,NDY,2)
c             ELSE
c#ifdef SFX1D              
c             CALL Get_Op_PSI_ SFX1D (Rho_PR,2)  !FFT 1D in Y direction
c#endif
c             ENDIF
c        write(6,*)"RTyy(i=",i,")= ",SUM(conjg(PSI(:,i))*(Rho_PR))
c        write(6,*)"h_ii (i=",i,")= ",Z1
c     . +SUM(conjg(PSI(:,i))*(VTRAP_EXT*PSI(:,i)))
c
c                  Rho_PR= psi(:,i)
c                  call Get_d_PSI(Rho_PR,3)
c                  call Get_d_PSI(Rho_PR,3)
c                  Z1=Z1+SUM(conjg(PSI(:,i))*(-Rho_PR/2))
c          write(6,*)"Tzz (i=",i,")= ",SUM(conjg(PSI(:,i))*(-Rho_PR/2))
c                  Rho_PR= psi(:,i)
c             IF(Time_DVRMETHODZ.ne.4) THEN
c                  CALL Get_Op_PSI(Rho_PR,Op_Z,NDZ,3)
c             ELSE
c#ifdef SFX1D              
c             CALL Get_Op_PSI_ SFX1D (Rho_PR,3)  !FFT 1D in X direction
c#endif
c             ENDIF
c        write(6,*)"RTzz(i=",i,")= ",SUM(conjg(PSI(:,i))*(Rho_PR))
c
c                  CALL VTRAP_EXT_TD(time)
c                  Z1=Z1+SUM(conjg(PSI(:,i))*(VTRAP_EXT*PSI(:,i))) 
c
c          write(6,*)"h_ii (i=",i,")= ",Z1 
!          write(6,*)"Kinetic,V orb (",i,")= ", Z1 
!     . ABS(SUM(conjg(PSI(:,i))*(-Rho_PR/2+VTRAP_EXT*PSI(:,i)))) 
!     . ,ABS(SUM(conjg(PSI(:,i))*(-Rho_PR/2))), 
!     . ABS(SUM(conjg(PSI(:,i))*(VTRAP_EXT*PSI(:,i))))
!          write(6,*)"<x/dx> (",i,")= ", SUM(conjg(PSI(:,i))*(Rho_PR))
c                enddo
c                  Rho_PR1= psi(:,1)
c                  Rho_PR2= psi(:,1)
c                  call Get_d_PSI(Rho_PR1,1)
c                  call Get_d_PSI(Rho_PR2,2)
c        write(6,*)"Px(1)= ",ABS(SUM(conjg(PSI(:,1))*(Rho_PR1)))
c        write(6,*)"Py(1)= ",ABS(SUM(conjg(PSI(:,1))*(Rho_PR2)))
c                   
c         ind=1
c         Do K=1,NDZ
c         Do J=1,NDY
c         Do I=1,NDX
c         ox=ort_X(i)
c         oy=ort_y(j)
c         oz=ort_z(k)
c         Rho_PR(ind)=(0.0d0,1.0d0)*(Rho_PR2(ind)*ox-Rho_PR1(ind)*oy)
c         ind=ind+1
c         Enddo
c         Enddo
c         Enddo
c        write(6,*)"Lz(1)= ",SUM(conjg(PSI(:,1))*(Rho_PR))
c          goto 1236
         open(unit=101,file='Rij.out',form='formatted')
         write(101,*)"============== h_ij r_ij ========================"
        Do i=1,MaxTrm1b
        P=TERM_INDEX_1B(I)
c================ Unpack cI cJ cK cL from P
        cK= INT(P/100)
        cJ= P-cK*100
       write(101,'(2i3,a10,F26.16,a10,F16.8,a10,F26.16,a10,F16.10)')
     . cK,cJ," Re[hij]=",Real(HIJ(I))," Im[hij]=",Dimag(HIJ(I)),
     ." Re[Rij]=",DReal(ZRIJ(I)), "Im[Rij]=",DImag(ZRIJKL(I))
        EnDdo
        
        write(101,*)"============== W_ijkl R_ijkl ====================="
         Do i=1,MaxTrm2b
          P=TERM_INDEX_2B(I)
c================ Unpack cI cJ cK cL from P
          cL= INT(P/1000000)
          cK= INT((P-cL*1000000)/10000)
          cJ= INT((P-cL*1000000-cK*10000)/100)
          cI= P-cL*1000000-cK*10000-cJ*100
       write(101,'(4i3,a8,F26.16,a8,F16.10,a12,F26.16,a12,F16.10)')
     . cI,cJ,cK,CL," Re[W]= ",Real(WIJKL(I)/xlambda0),
     ." Im[W]= ",Dimag(WIJKL(I)/xlambda0),
     ." Re[Rijkl]= ",DReal(ZRIJKL(I))," Im[Rijkl]= ",DImag(ZRIJKL(I))
          EndDo
         close(101) 
1236     continue
c         CALL MTRX_Diag
c=========================================================================
       CI_SCF=9999 ! CI_SCF=> DONE!!!
       call MPI_BCAST(CI_SCF,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c=========================================================================
c=========================================================================
c=========================================================================
333      continue
c=========================================================================
c     Lagrange multipliers !AIS 17JUL2014
       OPEN(unit=117,File='muij.out')
!       OPEN(unit=118,File='murij.out')
!Is used in LR         ZmuR=MAtMul( MAtmul(Zmu,Transpose(-Rho_ij)) !Rho_ij is not yet defined
       DO i=1,Morb
          DO j=1,Morb
!             WRITE(117,'(I4,I4,2E16.8)') i, j, ZmuR(i,j)
             WRITE(117,'(I4,I4,2E16.8)') i, j, Zmu(i,j)
          END DO
       END DO
       CLOSE(117)
!       CLOSE(118)
      call wr_info(6)
      call wr_info(9)
      call wr_info(10)
!      WRITE(9,*) "# This computation has been done in: ",TRIM(cwd)
!      WRITE(9,*)  banner(9)
!      write(9,'(a7,i3,a7,i8,a7,a20)')
!     . "# Morb=",Morb," Npar=",Npar," Job=",jtp
!      write(9,'(a13,a150)')"# V(x_y_z&t)= ",trim(stringV)
!      write(9,'(a17,a150)')"# W(R=|r1-r2|&t)= ",trim(stringW)
!      WRITE(10,*) "# This computation has been done in: ",TRIM(cwd)
!      write(10,'(a7,i3,a7,i8,a7,a20)')
!     . "# Morb=",Morb," Npar=",Npar," Job=",jtp
!      WRITE(10,*)  banner(9)
!      write(10,'(a13,a150)')"# V(x_y_z&t)= ",trim(stringV)
!      write(10,'(a17,a150)')"# W(R=|r1-r2|&t)= ",trim(stringW)
      close(9) 
      close(10) 

CCCCCCCCCCCCCCCC write output for testing purposes
      I = OMP_GET_NUM_THREADS()
      call basic_info(time,exc_tm,numprocs,I)

      DEALLOCATE(VIN)
      DEALLOCATE(VIN0)
      DEALLOCATE(VOUT)
      DEALLOCATE(Rho_PR)
      DEALLOCATE(PSI)
      deallocate(PSI1)
      deallocate(PSI2)
      deallocate(PSI_IN)
      deallocate(OPSI)
      deallocate(AuxPsi)
!      deallocate(IData)
!      deallocate(CData)
!      deallocate(RData)
!      deallocate(LData)
      write(6,*)" Master    ",MYID," is DONE after"
     ., MPI_WTIME(ierr)-exc_tm_new
         CALL  Close_binary
      return
      end subroutine Master_MCTDHB
