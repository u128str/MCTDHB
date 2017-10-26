C=============== Many-Body-WF-Analysis 
       subroutine Properties
       USE    PASS_ARG 
       USE   SHARED_DIMS
       USE   rR_hW
       USE   CI_All
       USE   PROP_MB
       USE   W_INTERPARTICLE
       USE   DVR_ALL
       USE   Parallel_Orb
       USE   CORRELATIONFUNCTIONS
       USE ANALYZER
       USE  CI_prod
       USE   LINEAR_RESPONSE
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
        COMPLEX*16, DIMENSION(:), ALLOCATABLE :: VIN 
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
      integer :: M,dilation,icntr,Iorb,Jorb,icnt
       COMPLEX*16, DIMENSION(100,100)  :: dOerr
       COMPLEX*16, DIMENSION(100)  :: dOerr_diag
c=============================================================
        COMPLEX*16, ALLOCATABLE ::  OPSI(:,:)
        COMPLEX*16, ALLOCATABLE ::  Rho_PR(:),PSI(:,:),FTPSI(:,:)
        COMPLEX*16, ALLOCATABLE ::  qq(:), fftqq(:,:)
        REAL*8, ALLOCATABLE ::  momX(:)
        COMPLEX*16, ALLOCATABLE ::  PSI_INTERPOLATED(:,:)
       REAL*8 :: time, tau,tau_new,ABMError,TolError_Total,t0
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
         allocate(Rho_ij(Morb,Morb),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in Rho_ij"
         allocate(Rho_ijkl(Morb,Morb,Morb,Morb),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in Rho_ijkl"
         allocate(H_ij(Morb,Morb),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in H_ij"

         allocate(W_ijkl(Morb,Morb,Morb,Morb),stat=ierr)
         if(ierr /= 0)write(*,*)"allocation error in W_ijkl"
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
              pathORBK='DATA'//path_sep//'orb_K'
         call mk_dir(pathORBK)
              pathg1RR='DATA'//path_sep//'g1_RR'
         call mk_dir(pathg1RR)
              pathg1KK='DATA'//path_sep//'g1_KK'
         call mk_dir(pathg1KK)
              pathLR='DATA'//path_sep//'getLR'
         call mk_dir(pathLR)
!        IF(PRINT_DATA.eqv..TRUE.) call print_rho(time,PSI,pathORB,Rho_PR)
c=========================== Cleaning some old files
      open(1234, file=trim(pathoRB)//'/'//trim('DNS_Minkovskii.dat'),
     . status='unknown', iostat=ierr)
      if ( ierr .eq. 0) close(1234, status='delete')
c============================Initialization of the MPI
       CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,ierr)
       CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,ierr)
       CALL MNGR_CI_Part(NUMPROCS)
       CALL MNGR_Orb_Part(NUMPROCS)
c==========================================================================
      time=0d0
      time_points_org=T_POINTS
      time=Time_Bgn
      write(6,*)"          PSI Data exist from time:",Time_Bgn
      write(6,*)"          User wants analysis from:",T_From
      time=MAX(T_From,Time_Bgn)
      write(6,*)"          I'll proceed        from:",time
      write(6,*)"          PSI Data exist till time:",time_psi_MAX 
      write(6,*)"          User wants analysis till:",T_TILL
      write(6,*)"                       delta step :",Time_print_step
      time_psi_MAX=MIN(time_psi_MAX,MIN(T_TILL,Time_Fnl))
      write(6,*)"          I'll proceed        till:",time_psi_MAX
      time_cic=time
      write(6,*) "The CI Part analysis starts at:",time_cic
      scl=(time_psi_MAX-time)/Time_print_step+1
      scl_cic=(time_psi_MAX-time)/(Time_print_step*time_ici_prt)+1
      
      write(6,*)"Number of available orbital time-points:",INT(scl)
      write(6,*)"User wants total number  of time-points:",T_Points
      T_POINTS=MIN(INT(scl),T_POINTS)
      write(6,*)"           I'll proceed  with  T_Points:",T_Points
      
      write(6,*)"Number of available CI time-points:",INT(scl_cic)

CCCCCCCCCCC ASSIGN STEP FOR ORBITAL PART >> scl
      if(INT(scl).ge.T_points) then 
         if(T_points.gt.1) then
            scl=NINT((scl-1)/(T_points-1))
         else if (T_points.eq.1) then
            scl=0
         endif
         write(6,*) " I change Delta t to:",scl*Time_print_step
       else if (T_POINTS.gt.scl) then
         T_POINTS=scl+1
         scl=1.d0
         write(6,*) " More points wanted than there are - 
     .         I use the number of available points:",scl
       endif 
      write(6,*)" the orbitals have the time-step:"
     .,scl*Time_print_step

       scl_cic=scl
       xlambda0=xlambda_0
        Prefac=JOB_PreFac
!==============================================================================
       call read_parse_V_W_Psi !parsing V_W_Psi_string.in file
       write(6,*) "File with user potential was parsed"
      open(unit=10,file='GO_PR.out',form='formatted')
      open(unit=11,file='nat_occs.out',form='formatted')

      icnt=0
      DO iii=1,T_POINTS 
        if (time.gt.t_till+1.d-9) then 
           time=time_psi_MAX
        endif

        if ((time_psi_max-time_cic).lt.
     .        (time_ici_prt*Time_print_step)) then 
                 time_cic=time_psi_MAX
        endif
      write(6,*) "ORB-time:",time

!> here ZRIJ and ZRIJKL at time time_psi are now read in from the binary file PSI_bin 
      CALL READ_BINARY(778, time, 1000,VIN,PSI)
      call Get_Full_Rij(Rho_IJ)
      call Diag_DNS
      write(11,2222)time,"  ",energy,"  ",(Nocc(i),"  ",i=1,Morb)
      call Get_Full_Rijkl(Rho_IJKL)

      IF (PROJ_X.eqv..TRUE.) then
         CALL Get_VtrapProjection(PSI,VTRAP_EXT,Rho_IJ,time,DIR)
      endif
!      write(*,*) NatVec
      IF (L_Z.eqv..TRUE.) then 
         CALL Get_LZ(PSI,Nocc,NatVec,time,Rho_IJ) 
      ENDIF
      call Diag_Geminal(time,NatGemOcc,1)
       
      write(10,2222)time,"  ",(NatGemOcc(i)/(Npar*(Npar-1)),
     &           "  ",i=1,(Morb*(Morb+1)/2))
2222  format((99(F21.16,a3)))
c================== Here everything is ready for PRINT ====================
       IF(DATA_PSI.eqv..TRUE.)  call print_rho(time,PSI,pathORB,Rho_PR) ! PRINT DATA-file
c================== Here I  PRINT FILE with DATA to plot 1D MInkovskii-like plot T.vs.X ====================
       IF(DIM_MCTDHB.eq.1) then 
        if(Mnkvsk_cmd.eqv..TRUE.)  then 
        call print_rho(time,PSI,pathORB,Rho_PR) ! PRINT DATA-file

        inquire(file=trim(pathoRB)//'/'//trim('DNS_Minkovskii.dat'),
     .  exist=Fexist)
        if (Fexist) then
        open(121, file=trim(pathoRB)//'/'//trim('DNS_Minkovskii.dat'), 
     .   status="old", position="append", action="write")
         else
         open(121, file=trim(pathoRB)//'/'//trim('DNS_Minkovskii.dat'),
     .   status="new", action="write")
         end if
          do i=1,NDX
         write(121,'(4(F12.6,a3))')
     .   ort_x(i),"  ",REAL(VTRAP_EXT(i)),
     .   "  ",Real(Rho_PR(i)),"   ",time
          end do
         write(121,'(a30)')"            "
         close(121)
         endif
       ENDIF
!===================================================================================================================
!============== Here we can transfer orbitals psi_i from any DVR  to equidistant FFT grid  only 1D so far AIS 21JUL2014
          Get_INTRPL=.False.
       IF(DIM_MCTDHB.eq.1) then 
        IF(Get_INTRPL.eqv..True.) then
       write(6,*) "Doing interpolation in 1D"
       Do i=1,Morb
       invec=psi(:,i)
       call interpolation1D(invec,outvec,ort_X,newgrid,NDX,NDX)
       PSI_INTERPOLATED(:,i)=outvec
       enddo   
       write(6,*) "DONE interpolation in 1D"
        endif
       endif
!===============================================================================================
       IF(NO_x.eqv..TRUE.) then 
       write(6,*) "rho(R) at T=",time
          call get_1BodyDiag(time,PSI,Rho_ij,1,1,pathORB) !STR 2014 ok 1D,2D
       endif
       IF (lossops.eqv..TRUE.) then
          write(*,*) 'I now calculate the N=2 lossop expectation values'
          call lossops_two_bosons(time,PSI,rho_ijkl,border)
       endif 
!================= Prints CORRELATION g1(x,x) in Coordinate (X) Space 
       IF(DNS_XX.eqv..TRUE.) then
       IF(DIM_MCTDHB.ne.1) STOP "g1(R,R') in D>1 ? -- check by yourself"
       write(6,*) "g1(R,R') at T=",time
          call get_correlations(Time,PSI,1,1,rho_ij,rho_ijkl,pathg1RR)
       endif
!NO_K ================= To Print diagonal one-particle Density in Momentum (K) Space  DNS(K) via FFT
       IF(NO_k.eqv..TRUE.) then               
       write(6,*) "rho(K) at T=",time
!STR change to TiLTN        if (MKLFT.eqv..FALSE.) then
         if (DiLTN.eqv..TRUE.) then  !Dilation is active for perastrian FFT ONLY
      IF(DIM_MCTDHB.ne.1) STOP "Dilation D>1 is not yet implemented"
           call pedestrian_FT(Psi,FTPSI,1)
           if(icnt.eq.0) call get_dilation(dilation,FTPSI,rho_ij)
!                dilation=1 ! STR to pu off scaling ==dilation
           call pedestrian_FT(Psi,FTPsi,dilation) 
           call get_1BodyDiag(time,FTPSI,rho_ij,2,dilation,pathORBK) 
         else if ((DiLTN.eqv..FALSE.).and.
     .            (ZeroPadding2D.eqv..FALSE.)) then
           dilation=1
!           call MKL_FT(Psi,FTPsi)  !Axel's MKLFFT
           call  Get_FFTPSI(Psi,FTPsi,1) !STR FFT
           call  Get_Shift_Zero_K_FFT(FTPsi)
       
           call get_1BodyDiag(time,FTPSI,rho_ij,3,dilation,pathORBK) 
         else if ((DiLTN.eqv..FALSE.).and.
     .           (ZeroPadding2D.eqv..TRUE.)) then
              STOP "Dilation in 2D has toobe checked !!!!"
           dilation=1
           call get_2D_1BodyDiag_dilated(time,rho_ij,PSI,dilation2D) !not checked yet STR 2014
         endif
       endif
!================= Prints CORRELATION g1(k,k) in Momentum (K) Space 
           IF(DNS_KK.eqv..TRUE.) then
       IF(DIM_MCTDHB.ne.1) STOP "g1(K,K') in D>1 ? -- check by yourself"
              dilation=1
!              call pedestrian_FT(Psi,FTPSI,1) !ORG 
           call  Get_FFTPSI(Psi,FTPsi,1) !STR FFT
           call  Get_Shift_Zero_K_FFT(FTPsi)
              call get_correlations(Time,FTPsi,2,dilation,
     .                                  rho_ij,rho_ijkl,pathg1KK)
            endif
! NONESCAPE PROBABILITY
      IF(Pnot.eqv..TRUE.) then
!       subroutine density_nonescape(xi,xf,rho_jk,time,PSI)
        call density_nonescape(xstart,xend,Rho_ij,time,PSI)
      ENDIF       
! correlation functions in xini,xfin and kxini,kxfin with xpts/kpts points
           write(6,*)"before CORR2RESTR FFT "
      if (CORR2RESTR.eqv..TRUE.) then
!         call MKL_FT(Psi,FTPsi) 
           call  Get_FFTPSI(Psi,FTPsi,1) !STR FFT
           call  Get_Shift_Zero_K_FFT(FTPsi)
      call Get_Full_Rij(Rho_IJ)
      call Get_Full_Rijkl(Rho_IJKL)
      call get_corr_two_restricted(time,PSI,1,
     .                        1,rho_ij,rho_ijkl,
     .                        xini2,xfin2,xpts2,0.d0,0.d0,1) 
      endif
      if (CORR2RESTRMOM.eqv..TRUE.) then
!         call MKL_FT(Psi,FTPsi) 
           call  Get_FFTPSI(Psi,FTPsi,1) !STR FFT
           call  Get_Shift_Zero_K_FFT(FTPsi)
         call Get_Full_Rij(Rho_IJ)
         call Get_Full_Rijkl(Rho_IJKL)
         call get_corr_two_restricted(time,FTPSI,2,
     .                      1,rho_ij,rho_ijkl,
     .                      0.d0,0.d0,1,kxini2,kxfin2,kpts2) 
      endif
      if (CORR1RESTR.eqv..TRUE.) then    
!         call MKL_FT(Psi,FTPsi) 
           call  Get_FFTPSI(Psi,FTPsi,1) !STR FFT
           call  Get_Shift_Zero_K_FFT(FTPsi)
         call Get_Full_Rij(Rho_IJ)
         call Get_Full_Rijkl(Rho_IJKL)
         call get_corr_one_restricted(time,PSI,1,
     .                            1,rho_ij,rho_ijkl,xini1,xfin1,xpts1,
     .                            0.d0,0.d0,1) 
      endif
      if (CORR1RESTRMOM.eqv..TRUE.) then    
!         call MKL_FT(Psi,FTPsi) 
           call  Get_FFTPSI(Psi,FTPsi,1) !STR FFT
           call  Get_Shift_Zero_K_FFT(FTPsi)
         call Get_Full_Rij(Rho_IJ)
         call Get_Full_Rijkl(Rho_IJKL)
         call get_corr_one_restricted(time,FTPSI,2,
     .                            1,rho_ij,rho_ijkl,0.d0,0.d0,1,
     .                            kxini1,kxfin1,kpts1) 
      endif
      if (MOMSPACE2D.eqv..TRUE.) then
       if (ZeroPadding2D.eqv..TRUE.) then
!         call twoD_zero_padding(PSI,FTPSI,dilation2d) !not checked
       elseif (ZeroPadding2D.eqv..FALSE.) then
!         call MKL_FT(Psi,FTPsi)
           call  Get_FFTPSI(Psi,FTPsi,1) !STR FFT
           call  Get_Shift_Zero_K_FFT(FTPsi)
         dilation2d=1 
       endif
         call Get_Full_Rij(Rho_IJ)
         call Get_Full_Rijkl(Rho_IJKL)
         call get_corr_slice(time,FTPSI,2,rho_ij,rho_ijkl,
     . x1const,x1slice,y1const,y1slice,
     . x2const,x2slice,y2const,y2slice,dilation2d) 
      endif
      if (REALSPACE2D.eqv..TRUE.) then
         call Get_Full_Rij(Rho_IJ)
         call Get_Full_Rijkl(Rho_IJKL)
         call get_corr_slice(time,PSI,1,rho_ij,rho_ijkl,
     . x1const,x1slice,y1const,y1slice,
     . x2const,x2slice,y2const,y2slice,1) 
      endif
       
!       CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,ierr)
!       CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,ierr)
!       CALL MNGR_CI_Part(NUMPROCS)
!       CALL MNGR_Orb_Part(NUMPROCS)
c===================================Wsl  analysis=================
      write(6,*) "WSL-time:",time, get_WSL
      IF (get_WSL.eqv..TRUE.) THEN
         write(6,*)" WSL is active "
         write(6,*)"==================================================="
        xlambda0=xlambda_0
        Prefac=JOB_PreFac
!        CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYID,ierr)
!        CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,ierr)
!       CALL MNGR_CI_Part(NUMPROCS)
!       CALL MNGR_Orb_Part(NUMPROCS)
!       CALL Get_InterPart(time)
       zmu=0d0
       write(6,*)"Zmu  at Time",time,size(zmu)
       Prefac=JOB_PreFac
       ENRG_EVAL=.FALSE.
       write(6,*) "NRM", DREAL(ZDOTC(NDX,PSI(:,1),1,PSI(:,1),1))
       CALL FUNCSTR(time,PSI,OPSI) 
       ENRG_EVAL=.FALSE.
       CALL Get_h_W(PSI,time)
       write(6,*)"Wsl  at Time",time
       W=Zero
       E=SUM(AllHIJ*AllZRIJ)
       W=SUM(ZRIJKL*WIJKL*(DREAL(RESCALE_2B)))
     . +SUM(Conjg(ZRIJKL)*Conjg(WIJKL)*(DIMAG(RESCALE_2B)))
        Energy=Real(E+W)
        WRITE(6,*) 'CI Energy in Wsl', Real(E+W),E,W
!        WRITE(6,*) 'in Wsls Zrij', ALLZrij
c====================================================
!         Get_WSL(WSL,psiS,psiL) 
!         COMPLEX*16,DIMENSION(NDX*NDY*NDZ) :: psiS,psiL
!         COMPLEX*16,DIMENSION(NDX*NDY*NDZ) :: WSL
!========= HERE I PRINT WSL S L in the same format as in time.dat file but instead of rho(r,r) in 8 and 9 colomns  it now writes real(WSL) and Imag(WSL)
         call print_WSL(time,PSI,1,1)
         call print_WSL(time,PSI,2,2)
         write(6,*)" WSL is done"
         write(6,*)"==================================================="
c=====================================================
      ENDIF  ! WSL part

!====================================CI Part=================================
         CIC=.FALSE.
      write(6,*) "CI-time:",time_cic
    
      if ((DATA_CIc.eqv..TRUE.).and.
     .   (Time_cic.le.(Time_PSI_MAX+1.d-8))) then
C         write(6,*) "I'll read the CI-binary with t=:",time_cic
         CALL READ_BINARY(777, time_cic, 1000,VIN,PSI)
         call print_cic(time_cic,VIN,pathCI) ! PRINT DATA-fila
         CIC=.TRUE.
      endif
 

      time_cic=time_cic+(Time_print_step*time_ici_prt)

      if (time_cic.gt.t_till) then
         DATA_CIC=.FALSE.
c         write(6,*) 'LAST STEP FOR CI:',
c     .                time_cic-(Time_print_step*time_ici_prt)
      endif


!      time=time+Time_print_step*scl
!      icnt=icnt+1
c      write(*,*) "increased CI-time to:",time_cic,"and the 
c     .  orbital time to:",time
c===================================Linear Response analysis=================
c===================================Test: Diagonalization of CI-Matrix=======
      IF (get_LR.eqv..TRUE.) THEN
           
      IF (CI_PRD.eqv..FALSE.) THEN
      WRITE(6,*)"Recompile the MCTDHB package with CI_PRD=.TRUE."
      WRITE(6,*)"in source/MODULES_ALL_allocate.f90 and try again"
      stop 
      ENDIF
      IF (CIC.eqv..FALSE.) THEN
      WRITE(6,*)"DATA_CIc=",DATA_CIc
      WRITE(6,*)"CIC=",CIC
      WRITE(6,*)"put  DATA_CIc=.T. in properties.in  to fill CI: VIN"
      WRITE(6,*)" and execute  properties_MB.exe again"
      stop 
      ENDIF
!         write (pathLRtime, "(A1,F10.6)")pathLR,"T",time
!           pathLRtime='DATA/getLR/T'trim(pathLRtime)
!         call mk_dir(pathLRtime)
!          print *, pathLRtime
         CALL Get_LinResp(time,PSI,VIN,pathLR) !Call for LR
       write(6,*)"LR DONE"
      ENDIF
c=================================== End of Linear Response analysis=================

      time=time+Time_print_step*scl
      icnt=icnt+1
c      write(*,*) "increased CI-time to:",time_cic,"and the 
c     .  orbital time to:",time


      if (time.gt.t_till+1.d-8) then 
         goto 201
      endif
      end do
 201   continue
      close(10)
      close(11)

      CALL Close_binary 


c=========================================================================
      DEALLOCATE(VIN)
      DEALLOCATE(qq,fftqq)
      DEALLOCATE(Rho_PR)
      DEALLOCATE(OPSI,PSI)
      deallocate(PSI_INTERPOLATED)
      DEALLOCATE(Rho_ij)
      DEALLOCATE(Rho_ijkl)
      DEALLOCATE(NatGemOcc)
      return
      end subroutine Properties
c=========================================================================
c=======================================================================
       subroutine Get_Full_Rij(Rho_JK)
       USE   PASS_ARG
       USE   SHARED_DIMS
       USE   rR_hW
       USE   CI_All
       USE   PROP_MB
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
       USE   PROP_MB
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
       USE   PROP_MB
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
       USE   PROP_MB
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
       USE   PROP_MB
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
   

