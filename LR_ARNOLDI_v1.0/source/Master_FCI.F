C=============== Assistent MPI subroutive -assists in H|psi> multiplication
       subroutine Master_FCI
       USE    PASS_ARG
       USE   SHARED_DIMS
       USE   rR_hW
       USE   CI_All
       Use OMP_lib

          USE   W_INTERPARTICLE
          USE   DVR_ALL
c          USE   Integrator
!          USE   rR_hW
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
        SAVE
        COMPLEX*16, DIMENSION(:), ALLOCATABLE :: VIN 
        COMPLEX*16, DIMENSION(:), ALLOCATABLE :: VIN0
        COMPLEX*16, DIMENSION(:), ALLOCATABLE :: VOUT
c====================== For SIL
        COMPLEX*16 :: zrho=ZERO,Z,Z1,zerr,zerrci,zerr_new
        COMPLEX*16 :: Escale
c==============================================
      real*8, dimension(12) :: time_rec(12)=0.0d0
      real*8 :: time_ci1=0.0d0,time_ci2=0d0,time_ci3=0d0
      real*8 :: time_psi1=0.0d0,time_psi2=0d0,time_psi3=0d0
c==============================================
      real*8  xnrm,Error_SIL,xerr,eerr,Tmax,xerr_td
      LOGICAL  CNV,SIL,Fexist,Fexist1
      real*8   FKTR,CNK,DZASUM,DZNRM2,DNRM2,DDOT
      real*8   t_before_SIL,t_after_SIL
c===================== t E S T ========================================
c      integer, DIMENSION(Morb-1) :: Ivec
c      integer, DIMENSION(Morb)   :: Nvec
      integer :: M
c                    RHF Variables
       COMPLEX*16, DIMENSION(100,100)  :: dOerr
       COMPLEX*16, DIMENSION(100)  :: dOerr_diag
c=============================================================
c       COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb)  :: PSI,PSI1,PSI2
c       COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb)  :: OPSI
c       COMPLEX*16, DIMENSION(NDX*NDY*NDZ*Morb*8)  :: AuxPsi
       COMPLEX*16, ALLOCATABLE ::  Rho_PR(:)
        COMPLEX*16, ALLOCATABLE ::  PSI(:,:)
        COMPLEX*16, ALLOCATABLE ::  PSI1(:,:)
        COMPLEX*16, ALLOCATABLE ::  PSI2(:,:)
        COMPLEX*16, ALLOCATABLE ::  OPSI(:,:)
        COMPLEX*16, ALLOCATABLE ::  AuxPsi(:)
        COMPLEX*16, ALLOCATABLE ::  PSI_IN(:,:)
       LOGICAL :: NewStart
       REAL*8 :: time, tau,tau_new,ABMError,TolError_Total,t0
       REAL*8 :: time_next,print_next,tau_in,print_step
       REAL*8 x0in, gr
       REAL*8 :: Energy_old=0.0d0,exc_tm=0d0,exc_tm_new=0d0
       COMPLEX*16 :: ZDOTC
       INTEGER ::  DM
       character*20 jtp
       character*10 path,pathCI,pathORB
c===========================================================
       INTEGER :: I,Nproc,istep
c================================================
       integer :: GetInd
c=============================================================
      external FKTR,CNK
      external Expect_fiXfi
      external DZASUM,DZNRM2,DNRM2,ZDOTC,DDOT,ZSCAL,Integrator_ABM
c=============================================================
          exc_tm=MPI_WTIME(ierr)
          exc_tm_new=MPI_WTIME(ierr)
c=============================================================
c         call MPI_INIT(ierr)
         call MPI_COMM_RANK(MPI_COMM_WORLD,MYID,ierr)
         call MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,ierr)
         write(6,*)" Process ",MYID," of ", numprocs
         CALL MNGR_CI_Part(NUMPROCS)
         CALL MNGR_Orb_Part(NUMPROCS)  
c        CALL MNGR_Orb_Part(1)  
c=============================================================
c       Nconf=NINT(CNK(Npar+Morb-1,Npar))
c===========================================================
         write(6,*)" Master task for ",Npar," particles and ", Morb,
     $" orbitals is  C(N+M-1,N) long:", Nconf
c            pause
C             Nvec=0
C             Nvec(1)=100
C         write (6,*) GetInd(100,4,Nvec)
C             Nvec=0
C             Nvec(4)=100
C         write (6,*) GetInd(100,4,Nvec)
C             Nvec=0
c             Nvec(1)=8
c             Nvec(2)=8
c             Nvec(3)=4
c         write (6,*)Nvec," Conf N.", GetInd(20,3,Nvec)
c=============================================================
c       Nconf=NINT(CNK(Npar+Morb-1,Npar))
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
c       COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb)  :: PSI,PSI1,PSI2
c       COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb)  :: OPSI
c       COMPLEX*16, DIMENSION(NDX*NDY*NDZ*Morb*8)  :: AuxPsi
         allocate(Rho_PR(NDX*NDY*NDZ),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in Pho_PR"
         allocate(PSI(NDX*NDY*NDZ,Morb),stat=ierr)
c      write(6,*) " MEM for PSI  is ok!", SIZE(PSI),Morb*NDX
         if(ierr /= 0)write(*,*)"allocation error in PSI"
         allocate(PSI1(NDX*NDY*NDZ,Morb),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in PSI1"
         allocate(PSI2(NDX*NDY*NDZ,Morb),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in PSI2"
         allocate(OPSI(NDX*NDY*NDZ,Morb),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in OPSI"
         allocate(AuxPsi(NDX*NDY*NDZ*Morb*8),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in AuxPsi"
         allocate(PSI_IN(NDX*NDY*NDZ,Morb),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in PSI_IN"
c=========================================================================
            write(6,*)"Memory allocation is ok"
c=========================================================================
C         N=Npar
C         M=Morb
C         jj=140
Cc         pause
Cc        call GET_Ind_vs_ii(jj,N,M,Ivec,Nvec)
C        call GetCIJKL1body_5_OMP(0,VIN)
C        stop
c=========================================================================
c                  W O R K I N G    P A R T 
c=========================================================================
c=========================================================================
      xlambda0=xlambda_0
          write(6,*) " LAMBDA_0=",xlambda0
c             pause
c      xlambda0=0.0000000001d0
c======================== G U E S S ===================================
       IF(ORB_Binr) then
       CALL Open_binary !Needed for restarts
       CALL Read_Header(778,'scan')
       CALL Read_Header(777,'scan')
       CALL READ_BINARY(778, Binary_Start_Point_t, 1000,VIN,PSI)
       CALL READ_BINARY(777, Binary_Start_Point_t, 1000,VIN,PSI)
       CALL Close_binary !Needed for restarts
       ELSE
       CALL Guess_PSI(PSI)
       CALL Guess_CI(VIN)
       ENDIF
c       CALL Guess_PSI(PSI)
c       CALL Guess_CI(VIN)
!        write(6,*) " Guess CI is ok",VIN(1),VIN(Nconf)
!        write(6,*) " Geess CI is ok",SUM(CONJG(VIN)*VIN)
!        write(6,*) " Geess PSI is ok",SUM(CONJG(PSI(:,1))*PSI(:,1))
!        write(6,*) " Geess PSI is ok",SUM(CONJG(PSI(:,2))*PSI(:,2))
!        write(6,*) " Guess PSI is ok",SUM(CONJG(PSI(:,2))*PSI(:,1))
c================================================================
        CALL Open_binary ! Binary files are rewritten !!!!!!!!!!!!!!!!
        CALL Write_Header(777) !NAMELIST IS OVERWRITTEN to CIc binary file
        CALL Write_Header(778) !NAMELIST IS OVERWRITTEN to PSI binary file
c================================================================
c         Initialization of the arrays
c=========Initial TIME can be overwritten here ===========================
        time=Time_Bgn
!        print *, "Before Get_h_W 1"
       CALL Get_h_W(PSI,time)
       CALL Get_r_R(VIN,VOUT,time)
       CALL Get_h_W(PSI,time)

!       CAll Get_BH_integrals
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
c=================== Second order integration scheme ====================
        tau=Time_tau
c========================================================================
         Tmax=600.0d0
         Tmax=Time_Max
         Tmax=Time_Fnl
c========================================================================
       tau_in=tau
       print_step=Time_print_step
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
             pathORB='DATA/orb_R'
        call mk_dir(pathORB)
             pathCI='DATA/CIcnf'
        call mk_dir(pathCI)
       IF(PRINT_DATA.eqv..TRUE.) call print_rho(time,PSI,pathORB,Rho_PR)
       IF(PRINT_DATA.eqv..TRUE.) call print_cic(time,VIN,pathCI)
c========================================================================
       write(6,'(125a1)')("=",jj=1,125)
       write(6,*)"================================Propagation started "
       write(6,'(125a1)')("=",jj=1,125)
c========================================================================
        open(unit=10,file='NO_PR.out',form='formatted')
c========================================================================
                    MacroIter_Max=400000
          exc_tm=MPI_WTIME(ierr)-exc_tm
         write(6,*)" Initialization time took",exc_tm,"[sec]"
         Do istep=1,MacroIter_Max
c========================================================================
         time_rec(1)=MPI_WTIME(ierr) 
c========================================================================


       write(6,'(125a1)')("=",jj=1,125)
       write(6,'(a20,a15,I8,a12,F16.10,a3,F16.10,a3,F16.10,a3,F16.10,a3)
     .')jtp,
     ."Iteration:",istep,"Time: [",Time_Bgn,"->",time," + ",tau
     .,"->",Time_Fnl," ]"
       write(6,'(a35,F26.16,a9,i10,a10,i10,a1,i10,a1,i10,a1)')
     ."Input orbital energy E(t+0):",Energy,"( CI Dim:",Nconf,
     .")(ORB Dim:",Morb,"*",NDX*NDY*NDZ,"=",Morb*NDX*NDY*NDZ,")"
c========================================================================
         time_rec(2)=MPI_WTIME(ierr)
c========================================================================
        VIN0=VIN
c========================================================================
c            I    CI P_a_r_t |                  0->tau
c========================================================================
       t0=tau
c=========================================================================
       CI_SCF=0 ! CI_SCF=0 -> CI PART
       call MPI_BCAST(CI_SCF,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       call Share_H_W
c=========================================================================
       Error_SIL=TolError_Total/1.0d0
       t_before_SIL=t0
       MAXIT= 20
       CALL SIL_PROPG(t0,VIN,VOUT,Error_SIL,MAXIT)
       t_after_SIL=t0
       write (6,*) "t_before_SIL",t_before_SIL
       write (6,*) "t_after_SIL",t_after_SIL
       tau=t0
       VIN=VOUT
!       write (6,*) HIJ,"Time=",time,WIJKL
c========================================================================
       time_rec(10)=MPI_WTIME(ierr)
       time=time+tau
       CALL Get_r_R(VIN,VOUT,time)                                     ! rij(tau) and Rijkl(tau)
c=========================================================================
       time_rec(12)=MPI_WTIME(ierr)
c=========================================================================
       write(6,'(a35,F26.16,a9,i10,a12,G17.10,a17,I1)')
     ."OUT  CI  energy E(t+tau):",Energy,
     ."   N = ",Npar," l0*(N-1)=",xlambda0*(Npar-1)
     .," kind of W(x-x'):",Wxx_TYPE
       IF(PreFac.eq.(-1.0d0,0d0).or.PreFac.eq.(1.0d0,0d0))
     .zerr=Energy-Energy_old
       eerr=Energy-Energy_old
       write(6,'(2(a35,F26.16))')"Delta E     : +/-",Energy-Energy_old
     .                          ,"Error due to  dE:",xerr_td
       write(6,'(a35,F26.16)')"Error E     : +/-",TolError_Total
       Energy_old=Energy
       write(6,'(a35,2F26.16,a5)')"CI  Error   : +/-",ABS(zerrci)
     .,100*ABS(zerrci)/TolError_Total,"%"
       write(6,'(a35,2F26.16,a5)')"ORB Error   : +/-",  xerr
     .,100*xerr/TolError_Total,"%"
c         write(6,*)" CI norm: Vout",DZNRM2(Nconf,VOUT,1)
       write(6,'(a10,a70)')" ",
     ."===================== Natural Occupation numbers ============="
c       write(6,'(115a1)')("=",jj=1,115)
      write(6,'(a15,20(G19.12,a3))')"New NO's:",(Nocc(i)," | ",i=1,Morb)
      write(6,'(a15,20(G19.12,a3))')"<fi_WR|X|fi_WR>:",
     . (Dreal(SUM(PSI(:,i)*ort_X*dconjg(PSI(:,i))))," | ",i=1,Morb)
        call Expect_fiXfi(PSI)
      write(6,'(a15,20(G19.12,a3))')"<fi_NO|X|fi_NO>:",
     . (Dreal(NO_Expectation_x(i))," | ",i=1,Morb)
       write(6,'(a10,a70)')" ",
     ."=====================                            ============="

      tau_new=Time_tau

      IF(ABS(tau_new).ge.1.0d0) tau_new=1.0d0*tau_new/ABS(tau_new)

c       write(6,*)time,"====time print next=====",print_next,"tau",tau
c       pause 
      if(ABS(time-print_next).le.1.0d-8) then
      print_next=print_next+print_step
      IF(PRINT_DATA.eqv..TRUE.) call print_rho(time,PSI,pathORB,Rho_PR)
      call Write_PSI_string(time,istep,Energy,PSI)
      Call flush(778)
c      if(ici.eq.10) then
      if(ici.eq.ici_prt) then
      IF(PRINT_DATA.eqv..TRUE.) call print_cic(time,VIN,pathCI)
c       write(6,*)"Before wr",time,istep,Energy,VIN
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
cORG      tau=print_next-time
c      print_next=print_next+print_step
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
     .,(Dreal(NO_Expectation_x(i)),"  ",i=1,Morb)
     .,(Dreal(SUM(PSI(:,i)*ort_X*dconjg(PSI(:,i)))),"  ",i=1,Morb)
2222  format((78(G26.16,a3)))
        call flush(10) 
c          pause
c       tau=tau_in
c       tau=-tau
c       PreFac=-1*Prefac
c         write(6,*)" Time:",time," New tau",tau," Old tau",tau_new
c           pause
c       write(6,'(a75,F16.10,a20,F16.10)')"Time:",time
c     .,"            New tau=",tau
c         IF(ABS(zerr).ge.epsilon(0.d0))THEN
c       write(6,'(a39,E26.16,2(a12,E12.5))')
c     ."Overall ||PSI|| error  of this step is=",REAL(ABS(zerr))
c     .," Err  CI[%]",100*ABS(zerrci)/TolError_Total
c     .," Err ORB[%]",100*(ABS(zerr)-ABS(zerrci))/TolError_Total
c          ENDIF
c       write(6,*)

c=========================================================================
cc       time_ci1=time_rec(3)-time_rec(2)
c       time_ci2=time_rec(10)-time_rec(9)
c       time_ci3=time_rec(4)-time_rec(3)
c     .         +time_rec(7)-time_rec(6)
c     .         +time_rec(12)-time_rec(11)
c       write(6,'(a15,F16.10,3(a10,F16.10))')
c     ."CI  Time:",time_ci1+time_ci2+time_ci3,
c     ."  CI   I:",time_ci1,
c     ."  CI  II:",time_ci2,
c     ." oth. CI:",time_ci3
cc=========================================================================
c       time_psi1=time_rec(5)-time_rec(4)
c       time_psi2=time_rec(8)-time_rec(7)
c       time_psi3=time_rec(9)-time_rec(8)
cc     .          +time_rec(2)-time_rec(1)

c       write(6,'(a15,F16.10,4(a10,F16.10))')
c     ."ORB Time:",time_psi1+time_psi2+time_psi3,
c     ."  ORB   I:",time_psi1,
c     ."  ORB  II:",time_psi2,
c     ."  ORB III:",time_psi3,
c     ." oth. ORB:",time_rec(2)-time_rec(1)
cc=========================================================================
c       exc_tm=exc_tm+time_rec(12)-time_rec(1)
c       write(6,'(a15,F16.10,a20,F16.10)')
c     ."Itr Time:",time_rec(12)-time_rec(1)," execution time:",exc_tm
c=========================================================================
c      if(istep.ge.1) goto 1235
c      if(time.ge.165.36560d0) goto 1235
c      if(time.ge.132.29249d0) goto 1235
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
      if(time.ge.Tmax) goto 1235
c      if(Energy.le.10.84082095058792d0) goto 1235
       ENDDO
c=================== END Second order integration scheme ====================
1235    continue
c=========================================================================
       CI_SCF=9999 ! CI_SCF=> DONE!!!
       call MPI_BCAST(CI_SCF,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c=========================================================================
c=========================================================================
c=========================================================================
333      continue
c=========================================================================
      close(10) 
c      call close_binary()
c=========================================================================
c      IF(ALLOCATED(PSI)) print *,"PSI is allocced",SIZE(PSI)
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
c      write(6,*)" Assistent ",MYID," is DONE"
      write(6,*)" Master    ",MYID," is DONE after"
     ., MPI_WTIME(ierr)-exc_tm_new
c      deallocate(MAP)
c      DEALLOCATE(PSI)
         CALL  Close_binary
      return
      end subroutine Master_FCI
