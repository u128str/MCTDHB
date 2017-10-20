      MODULE SHARED_DIMS
         IMPLICIT NONE
         SAVE 
         INTEGER, PUBLIC  :: DIM_MCTDHB
         INTEGER, PUBLIC :: Morb 
         COMPLEX*16, PUBLIC  ::    JOB_PreFac
         character*4, public ::  GUESS
         character*3, public ::  MB_JOB_TYPE
         logical, public ::      ORB_Binr
         REAL*8  ::  Binary_Start_Point_t
         character*18, public ::  Time_Res_Orb_File_Name
         character*18, public ::  Time_Res_CIc_File_Name
         LOGICAL,      public ::  CIc_Rest,ORB_Rest
         LOGICAL,      public ::  ORB_DIAG
         LOGICAL,      public ::  FFT_2D ! if ==true tru 2D FFT is working
         INTEGER, PUBLIC :: NDX,NDY,NDZ
         INTEGER, PUBLIC  :: ABMiter

         INTEGER, PUBLIC :: CI_SCF ! ADDED in CI case

!c===================================================
         REAL*8, public  :: Time_Bgn,Time_Fnl,Time_MAX, &
         Time_Print_Step,Time_ici_prt,Time_tau,Time_TolError_Total
         REAL*8, public  :: Time_Min_Occ_Allowed
         INTEGER, public :: Time_minsil,Time_maxsil
         character*18, public :: Time_intgr    ! ABM,OMPABM,BS
         INTEGER, public :: Time_intgr_order   ! Integrator  order for ABM<=8, for BS <=16
         INTEGER :: Intgr_total_steps   ! Total number of steps needed to solve the problem in ORB part
         REAL*8  :: Time_Intgr_MAX_Step      ! Maximal step allowed 
         LOGICAL,      public ::  PRINT_DATA
         INTEGER, PUBLIC  :: print_level ! Level of printed output 0 min 3 max
         REAL*8  ::  TD_Err_scale
!c===================================================
         COMPLEX*16, PARAMETER :: ZERO=(0.0d0,0.0d0)
         COMPLEX*16, PARAMETER :: ZONER=(1.0d0,0.0d0)
         COMPLEX*16, PARAMETER :: ZONEI=(0.0d0,1.0d0)
         REAL*8    , PARAMETER :: PI=3.141592653589793238462643d0
         REAL*8  ::  Energy,timeCPU(10)
         COMPLEX*16 ::  Prefac
!==========================================
!==========================================
         LOGICAL, PUBLIC :: LZ
         REAL*8, PUBLIC  :: OMEGAZ
         INTEGER, public :: STATE   ! Integrator  order for ABM<=8, for BS <=16
!====================================================
!       CHARACTER*3 ::  CHALL  = CHAR(13)//CHAR(11)//CHAR(0) 
       CHARACTER*3 ::  CHALL  = ""
#if  WIN32
       character(len=1), parameter :: path_sep='\'
!      
!#elif _WIN32
!      character(len=1), parameter :: path_sep='/' character(len=1), parameter :: path_sep='\'
#else
       character(len=1), parameter :: path_sep='/'
!#error "path_sep not defined. Set this constant for your system."
#endif
       
         character*66, DIMENSION(7) :: banner0=(/&
&" ____    ____    ______  _________  ______    ____  ____ ______    ",&
&"|_   \  /   _|.'' ___  ||  _   _  ||_   _ `. |_   ||  _||_   _  \  ",&
&"  |   \/   | / .''   \_||_/ | | \_|  | | `. \ | |__| |    | |_)  | ",&
&"  | |\  /| | | |            | |      | |  | | |  __  |    |  __''. ",&
&" _| |_\/_| |_\ `.___.''\   _| |_    _| |_.''/_| |  | |_  _| |__) | ",&
&"|_____||_____|`.____ .''  |_____|  |______.''|____||____||_______/ ",&
&"        http://mctdhb.org     http://mctdhb-lab.com/               "/)



         character*66, DIMENSION(44) :: banner=(/&
&" #===============================================================# ",&
&" #               __  __  ___ _____ ___  _  _ ___                 # ",&
&" #  Scientific  |  \/  |/ __|_   _|   \| || | _ ) (2006-present) # ",&
&" #              | |\/| | (__  | | | |) | __ | _ \    Germany     # ",&
&" #    Package   |_|  |_|\___| |_| |___/|_||_|___/   Heidelberg   # ",&
&" #      http://mctdhb.org    http://mctdhb-lab.com/              # ",&
&" #===============================================================# ",&
&" #   The Multiconfigurational Time-Dependent Hartree For Bosons  # ",&
&" #                  Major Version 3.3  (2006-2017)               # ",&
&" #===============  BBB: Be superB with the mctdhbB ==============# ",&
&" #                        Founders:                              # ",&
&" #                                                               # ",&
&" #     Alexej I. Streltsov, Ofir E. Alon, Lorenz S. Cederbaum    # ",&
&" #===============================================================# ",&
&" #            Created, developed and designed by                 # ",&
&" #                                                               # ",&
&" #           Alexej I. Streltsov (Project Leader)                # ",&
&" #          Alexej.Streltsov@pci.uni-heidelberg.de               # ",&
&" #                                                               # ",&
&" #                     Contributors:                             # ",&
&" #            Lorenz S. Cederbaum, Ofir E. Alon                  # ",&
&" #      Kaspar Sakmann, Axel U. J. Lode, Julian Grond            # ",&
&" #     Oksana I. Streltsova, Shachar Klaiman, Raphael Beinke     # ",&
&" #===============================================================# ",&
&" #                       Citation:                               # ",&
&" #   When citing the MCTDHB Package in the literature,           # ",&
&" #   please cite at least one of the papers 1), 2), or 3)        # ",&
&" #   as well as the Package 4):                                  # ",&
&" #                                                               # ",&
&" #    1) A. I. Streltsov, O. E. Alon, and L. S. Cederbaum,       # ",&
&" #       Phys. Rev. Lett. 99, 030402 (2007).                     # ",&
&" #                                                               # ",&
&" #    2) O. E. Alon, A. I. Streltsov, and L. S. Cederbaum,       # ",&
&" #       Phys. Rev. A 77, 033613 (2008).                         # ",&
&" #                                                               # ",&
&" #    3) A. I. Streltsov, O. E. Alon, and L. S. Cederbaum,       # ",&
&" #       Phys. Rev. A 81, 022124 (2010).                         # ",&
&" #                                                               # ",&
&" #    4) The Multiconfigurational Time-Dependent Hartree         # ",&
&" #       for Bosons Package, http://mctdhb.org,                  # ",&
&" #       A. I. Streltsov,  et al                                 # ",&
&" #                                                               # ",&
&" #         Major Version 3.2, Heidelberg, (2006-2015)            # ",&
&" #===============================================================# "/)

         END MODULE SHARED_DIMS

         MODULE CI_ALL
         USE SHARED_DIMS
         IMPLICIT NONE
         SAVE
         INTEGER, PARAMETER :: MaxUserCnf = 100
        INTEGER, DIMENSION(MaxUserCnf) :: UserCnfN !Added to read from V_W_Psi_string.in file
        real*8, DIMENSION(MaxUserCnf) :: UserCnfW  !Added to read from V_W_Psi_string.in file
         INTEGER, PUBLIC  :: Npar,Nconf
         INTEGER,  ALLOCATABLE :: MCNK(:,:)
         integer, ALLOCATABLE  :: MAPI(:)
         INTEGER, ALLOCATABLE :: Nmax(:)
!c==========================================================
         contains
         logical function CI_ALL_init()
         integer :: ierr,I,J
         real*8  :: CNK
         external CNK
         UserCnfN=0
         UserCnfW=0d0
          I=Npar+Morb-1
          J=Npar
          J=Morb-1
          Nconf=NINT(CNK(Npar+Morb-1,Npar))
        allocate(MCNK(0:I,0:J),stat=ierr)
        if(ierr /= 0) then
        write(*,*)"allocation error for MCNK"
        CI_ALL_init=.false.
        else
!        if(ierr == 0)write(*,*)"allocation ok for MCNK"
        CI_ALL_init=.true.
        endif
        allocate(Nmax(Morb),stat=ierr)
        if(ierr /= 0) then
        write(*,*)"allocation error for Nmax"
        CI_ALL_init=.false.
        else
!       if(ierr == 0)write(*,*)"allocation ok for Nmax"
        CI_ALL_init=.true.
        endif
        end function
!c==========================================================
         END MODULE CI_ALL

         MODULE Parallel_Orb
         IMPLICIT NONE
         SAVE
         INTEGER, PARAMETER :: MaxProc = 500
         INTEGER, DIMENSION(MaxProc) :: Proc_Job
         INTEGER, DIMENSION(MaxProc) :: Proc_Iorb_Starts
         INTEGER, DIMENSION(MaxProc) :: Proc_Iorb_Finish
         INTEGER, DIMENSION(MaxProc) :: Proc_WSL_Starts
         INTEGER, DIMENSION(MaxProc) :: Proc_WSL_Finish
         INTEGER, DIMENSION(MaxProc) :: Proc_From
         INTEGER, DIMENSION(MaxProc) :: Proc_Till
         INTEGER, DIMENSION(MaxProc) :: Orb_Block,Orb_Displ
         INTEGER :: JOB_TYPE
         END MODULE Parallel_Orb

         MODULE Parallel_CI
         IMPLICIT NONE
         SAVE
         INTEGER, PARAMETER :: MaxProc = 500
         INTEGER, DIMENSION(MaxProc) :: CI_Proc_Job
         INTEGER, DIMENSION(MaxProc) :: CI_Proc_Iorb_Starts
         INTEGER, DIMENSION(MaxProc) :: CI_Proc_Iorb_Finish
         INTEGER, DIMENSION(MaxProc) :: CI_Proc_Jdim_Starts
         INTEGER, DIMENSION(MaxProc) :: CI_Proc_Jdim_Finish
         INTEGER, DIMENSION(MaxProc) :: CI_Proc_From
         INTEGER, DIMENSION(MaxProc) :: CI_Proc_Till
         INTEGER, DIMENSION(MaxProc) :: MYID_TRM
         INTEGER :: EXCEPTIONAL_ID ! for V2 paralell menager ID of the process which works on both 1b and 2b
         INTEGER :: CI_JOB_TYPE
         END MODULE Parallel_CI



         MODULE W_INTERPARTICLE
         USE SHARED_DIMS
         IMPLICIT NONE
         SAVE 
         REAL*8  :: xlambda0,xlambda_0
         INTEGER :: Wxx_TYPE
         INTEGER ::  Time_DVRMETHODX,Time_DVRMETHODY,Time_DVRMETHODZ
         REAL*8 ::  Time_xint,Time_xfnl,Time_mass
         REAL*8 ::  Time_yint,Time_yfnl
         REAL*8 ::  Time_zint,Time_zfnl
         LOGICAL,   public :: WTD
!c====================== Wxx_TYPE=0 W(r,r')=Delta(r-r')
!c====================== Wxx_TYPE=1 W(r,r')=W(r-r')=WOp_X*WOp_Y*WOp_Z
!c         REAL*8, DIMENSION(NDX,NDX) :: WOp_X
!c         REAL*8, DIMENSION(NDY,NDY) :: WOp_Y
!c         REAL*8, DIMENSION(NDZ,NDZ) :: WOp_Z
         REAL*8, ALLOCATABLE :: WOp_X(:,:)
         REAL*8, ALLOCATABLE :: WOp_Y(:,:)
         REAL*8, ALLOCATABLE :: WOp_Z(:,:)
!c====================== Wxx_TYPE=2 W(r,r')=W(r-r')
!c         Complex*16, DIMENSION(NDX*NDY*NDZ*2-1) :: W2xx
         Complex*16, ALLOCATABLE :: W2xx(:)
!c====================== Wxx_TYPE=3 W(r,r')=W(r,r') - general case one has to store ALL matrix of Wxx !!!!!!  
!c====================== DOES NOT WELL TESTED !!!!!!!!!!
         INTEGER, PUBLIC :: NtVec
!c         Complex*16, DIMENSION(NtVec,NtVec) :: Wxx
!c         REAL*8, DIMENSION(NtVec,NtVec) :: FIJ
!c         REAL*8, DIMENSION(NtVec) :: F
         Complex*16, ALLOCATABLE :: Wxx(:,:)
         REAL*8, ALLOCATABLE :: FIJ(:,:)
         REAL*8, ALLOCATABLE :: F(:)
!c====================== Wxx_TYPE=4 W(r,r')=W(r-r') - only for an equidistant grid !!!!!!  For WSl as well
!C WXX TYPE 4 is for the IMEST ALGORITHM
         Complex*16, ALLOCATABLE :: W3xx(:),W3xxFFT(:)
                 
!c==========================================================
         contains
         logical function W_INTERPARTICLE_init()
         integer :: ierr,NtVec
!c====================== Wxx_TYPE=1 W(r,r')=W(r-r')=WOp_X*WOp_Y*WOp_Z
        IF(Wxx_TYPE == 1) THEN
        allocate(WOp_X(NDX,NDX),stat=ierr)
        if(ierr /= 0) then
        write(*,*)"allocation error for WOp_X"
        W_INTERPARTICLE_init=.false.
        else
        if(ierr == 0)write(*,*)"allocation ok for WOp_X"
        W_INTERPARTICLE_init=.true.
        endif
        allocate(WOp_Y(NDY,NDY),stat=ierr)
        if(ierr /= 0) then
        write(*,*)"allocation error for WOp_Y"
        W_INTERPARTICLE_init=.false.
        else
        if(ierr == 0)write(*,*)"allocation ok for WOp_Y"
        W_INTERPARTICLE_init=.true.
        endif
        allocate(WOp_Z(NDZ,NDZ),stat=ierr)
        if(ierr /= 0) then
        write(*,*)"allocation error for WOp_Z"
        W_INTERPARTICLE_init=.false.
        else
        if(ierr == 0)write(*,*)"allocation ok for WOp_Z"
        W_INTERPARTICLE_init=.true.
        endif
        ENDIF
!c====================== Wxx_TYPE=2 W(r,r')=W(r-r')
        IF(Wxx_TYPE == 2) THEN 
        allocate(W2xx(NDX*NDY*NDZ*2-1),stat=ierr)
        if(ierr /= 0) then
        write(*,*)"allocation error for W2xx"
        W_INTERPARTICLE_init=.false.
        else
        if(ierr == 0)write(*,*)"allocation ok for W2xx"
        W_INTERPARTICLE_init=.true.
        endif
        ENDIF
!c====================== Wxx_TYPE=3 W(r,r')=W(r,r') - general case one has to store ALL matrix of Wxx !!!!!!  
        NtVec=NDX
        IF(Wxx_TYPE == 3) THEN 
        IF(DIM_MCTDHB == 1)  NtVec=NDX
        IF(DIM_MCTDHB == 2)  NtVec=NDX*NDY
        allocate(Wxx(NtVec,NtVec),stat=ierr)
        if(ierr /= 0) then
        write(*,*)"allocation error for Wxx"
        W_INTERPARTICLE_init=.false.
        else
        if(ierr == 0)write(*,*)"allocation ok for Wxx"
        W_INTERPARTICLE_init=.true.
        endif
        ENDIF
!cccccccc       ALLOCATION OF ARRAYS FOR IMEST
       IF(Wxx_TYPE == 4) THEN
        allocate(W3xx(NDX*NDY*NDZ),stat=ierr)
        allocate(W3xxFFT(NDX*NDY*NDZ),stat=ierr)
        if(ierr /= 0) then
            write(*,*)"allocation error for W3xx"
            W_INTERPARTICLE_init=.false.
        else
        if(ierr == 0) write(*,*)"allocation ok for W3xx"
            W_INTERPARTICLE_init=.true.
        endif
        ENDIF
        end function
!c==========================================================
         END MODULE W_INTERPARTICLE

         MODULE DVR_ALL
         USE SHARED_DIMS
         USE W_INTERPARTICLE
         IMPLICIT NONE
         SAVE 

        REAL*8, ALLOCATABLE :: weight(:),weight_X(:),ort_X(:)
        REAL*8, ALLOCATABLE :: weight_Y(:),ort_Y(:),weight_Z(:),ort_Z(:)
        REAL*8, ALLOCATABLE :: Op_X(:,:), Op_Y(:,:),Op_Z(:,:)
        REAL*8, ALLOCATABLE :: dif1matX(:,:), dif1matY(:,:),dif1matZ(:,:)
        COMPLEX*16, ALLOCATABLE :: Op_Xcplx(:,:),VTRAP_EXT(:)
!c==========================================================
        REAL*8, ALLOCATABLE :: ort_kx(:),ort_ky(:),ort_kz(:)
!c========== Parsing of the parse_V_W_PSI.in file ===========
       character(len = 10),  dimension(11) :: variables
       real*8,               dimension(11) :: variablesvalues
       character (len = 5)  :: statusflag
       character*1255 stringV,stringW,stringMOM
       character*1255, DIMENSION(100)  :: stringPSI
          contains
         logical function DVR_ALL_init()
         integer :: ierr

           variables(1) = 'x'
           variables(2) = 'y'
           variables(3) = 'z'
           variables(4) = 'x1'
           variables(5) = 'y1'
           variables(6) = 'z1'
           variables(7) = 'x2'
           variables(8) = 'y2'
           variables(9) = 'z2'
           variables(10) = 't'
           variables(11) = 'r'
          stringV='Using defaults from VTRAP_EXT_TD.F'
          stringW='Using Defaults from Get_InterParticle.F'
        stringPSI='Using Defaults from Guess_PSI.F'
        stringMOM='Imprint defaults K from Guess_PSI.F'
        allocate(weight(NDX*NDY*NDZ),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in DVR_ALL"
        allocate(weight_X(NDX),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in DVR_ALL"
        allocate(   ort_X(NDX),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in DVR_ALL"
        allocate(weight_Y(NDY),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in DVR_ALL"
        allocate(   ort_Y(NDY),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in DVR_ALL"
        allocate(weight_Z(NDZ),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in DVR_ALL"
        allocate(   ort_Z(NDZ),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in DVR_ALL"

        IF(Time_DVRMETHODX==4) allocate(Op_X(1,1),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in DVR_ALL"
        IF(Time_DVRMETHODX/=4) allocate(Op_X(NDX,NDX),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in DVR_ALL"

        IF(Time_DVRMETHODY==4) allocate(Op_Y(1,1),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in DVR_ALL"
        IF(Time_DVRMETHODY/=4) allocate(Op_Y(NDY,NDY),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in DVR_ALL"


        IF(Time_DVRMETHODZ==4) allocate(Op_Z(1,1),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in DVR_ALL"
        IF(Time_DVRMETHODZ/=4) allocate(Op_Z(NDZ,NDZ),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in DVR_ALL"

        IF(Time_DVRMETHODX==4) allocate(dif1matX(1,1),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in DVR_ALL"
        IF(Time_DVRMETHODX/=4) allocate(dif1matX(NDX,NDX),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in DVR_ALL"

        IF(Time_DVRMETHODY==4) allocate(dif1matY(1,1),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in DVR_ALL"
        IF(Time_DVRMETHODY/=4) allocate(dif1matY(NDY,NDY),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in DVR_ALL"

        IF(Time_DVRMETHODZ==4) allocate(dif1matZ(1,1),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in DVR_ALL"
        IF(Time_DVRMETHODZ/=4) allocate(dif1matZ(NDZ,NDZ),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in DVR_ALL"

!        allocate(Op_Xcplx(NDX,NDX),stat=ierr)
!        if(ierr /= 0)write(*,*)"allocation error in DVR_ALL"
        allocate(VTRAP_EXT(NDX*NDY*NDZ),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in DVR_ALL"
       
        allocate(   ort_kx(NDX),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in DVR_ALL"
        allocate(   ort_ky(NDY),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in DVR_ALL"
        allocate(   ort_kz(NDZ),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in DVR_ALL"

        if(ierr /= 0) then
        write(*,*)"allocation error for DVR_ALL"
        DVR_ALL_init=.false.
        else
        if(ierr == 0)write(*,*)"allocation ok for DVR_ALL"
        DVR_ALL_init=.true.
        endif
        end function
         END MODULE DVR_ALL


         MODULE rR_hW
         USE SHARED_DIMS
         IMPLICIT NONE
         SAVE 
!c         INTEGER, PARAMETER :: Rdim = Morb*(Morb+1)/2
!c         COMPLEX*16, DIMENSION(Rdim) :: ZRIJ,HIJ
!c         COMPLEX*16, DIMENSION(Rdim) :: ZRIJ1
!c         COMPLEX*16, DIMENSION(Morb,Morb) :: AllZRIJ,AllHIJ,InvZRIJ
!c         COMPLEX*16, DIMENSION(Rdim*(Rdim+1)/2) :: ZRIJKL,WIJKL
!c         COMPLEX*16, DIMENSION(Rdim*(Rdim+1)/2) :: ZRIJKL1

        COMPLEX*16, ALLOCATABLE :: ZRIJ(:),HIJ(:),ZRIJ1(:)
        COMPLEX*16, ALLOCATABLE :: AllZRIJ(:,:),AllHIJ(:,:),InvZRIJ(:,:)
        COMPLEX*16, ALLOCATABLE :: AllZRIJ0(:,:)
        COMPLEX*16, ALLOCATABLE :: ZRIJKL(:),WIJKL(:),ZRIJKL1(:),WIJKL1(:)

!c         REAL*8, DIMENSION(Morb) :: Nocc
!c         COMPLEX*16, DIMENSION(Morb,Morb) :: NatVec
!c         COMPLEX*16, DIMENSION(Morb)  :: NO_Expectation_x=Zero

        REAL*8, ALLOCATABLE :: Nocc(:)
        COMPLEX*16, ALLOCATABLE :: NatVec(:,:),NO_Expectation_x(:)

!c         INTEGER, DIMENSION(Rdim)             :: TERM_INDEX_1B
!c         COMPLEX*16,DIMENSION(Rdim)           :: RESCALE_1B
!c         INTEGER, DIMENSION(Rdim)             :: TERM_REQ_1B

        INTEGER, ALLOCATABLE :: TERM_INDEX_1B(:),TERM_REQ_1B(:)
        COMPLEX*16, ALLOCATABLE :: RESCALE_1B(:)

!c         INTEGER, DIMENSION(Rdim*(Rdim+1)/2)  :: TERM_INDEX_2B
!c         COMPLEX*16,DIMENSION(Rdim*(Rdim+1)/2):: RESCALE_2B
!c         INTEGER, DIMENSION(Rdim*(Rdim+1)/2)  :: TERM_REQ_2B

        INTEGER, ALLOCATABLE :: TERM_INDEX_2B(:),TERM_REQ_2B(:)
        COMPLEX*16, ALLOCATABLE :: RESCALE_2B(:)


         LOGICAL :: ENRG_EVAL
         INTEGER :: MaxTrm1b,MaxTrm2b  !!!!!!! ADDED later!!!!!!!!!
!c         COMPLEX*16, DIMENSION(Morb,Morb) :: ZMU
        COMPLEX*16, ALLOCATABLE :: ZMU(:,:)
        COMPLEX*16, ALLOCATABLE :: ZMUR(:,:) !AIS 17JUL2014


         INTEGER, public :: Rdim,Rdim1
!c============================================================
         contains
         logical function rR_hW_init()
         integer :: ierr
         Rdim = Morb*(Morb+1)/2
         Rdim1=Rdim*(Rdim+1)/2
        allocate(ZRIJ(Rdim),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in rR_hW"
        allocate(HIJ(Rdim),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in rR_hW"
        allocate(ZRIJ1(Rdim),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in rR_hW"
        allocate(AllZRIJ(Morb,Morb),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in rR_hW"
        allocate(AllHIJ(Morb,Morb),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in rR_hW"
        allocate(InvZRIJ(Morb,Morb),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in rR_hW"
        allocate(AllZRIJ0(Morb,Morb),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in rR_hW"

        allocate(ZRIJKL(Rdim1),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in rR_hW"
        allocate(WIJKL(Rdim1),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in rR_hW"
        allocate(WIJKL1(Rdim1),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in rR_hW"
        allocate(ZRIJKL1(Rdim1),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in rR_hW"

        allocate(Nocc(Morb),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in rR_hW"
        allocate(NatVec(Morb,Morb),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in rR_hW"
        allocate(NO_Expectation_x(Morb),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in rR_hW"
        NO_Expectation_x=Zero

        allocate(TERM_INDEX_1B(Rdim),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in rR_hW"
        allocate(RESCALE_1B(Rdim),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in rR_hW"
        allocate(TERM_REQ_1B(Rdim),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in rR_hW"

        allocate(TERM_INDEX_2B(Rdim1),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in rR_hW"
        allocate(RESCALE_2B(Rdim1),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in rR_hW"
        allocate(TERM_REQ_2B(Rdim1),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in rR_hW"

!c         COMPLEX*16, DIMENSION(Morb,Morb) :: ZMU
        allocate(ZMU(Morb,Morb),stat=ierr)
        if(ierr /= 0)write(*,*)"allocation error in rR_hW"
        allocate(ZMUR(Morb,Morb),stat=ierr)
        if(ierr /= 0)write(*,*)"ZmuR allocation error in rR_hW" !AIS 17JUL2014


        if(ierr /= 0) then
        write(*,*)"allocation error for rR_hW"
        rR_hW_init=.false.
        else
        if(ierr == 0)write(*,*)"allocation ok for rR_hW"
        rR_hW_init=.true.
        endif
        end function
        END MODULE rR_hW

!c==================== Module for CI_Production with allocatable arrays!!!!
         module CI_prod ! All arrays are allocated in Parallel_MNGR_CI_Part_v2.F
                        ! allarrays are deallocated in ParaMain.F at the
                        ! end. Logical CI_Production is used to
                        ! specifify filling of the arrays and using of
                        ! them
          LOGICAL ::  CI_Production,CI_Production_1b,CI_Production_2b
          integer nXdim,nYdim
          INTEGER, ALLOCATABLE  :: Ind_CI_1b(:,:)
          INTEGER, ALLOCATABLE  :: Ind_CI_2b(:,:)
          Real*8,  ALLOCATABLE  :: Prefactors_1b(:,:)
          Real*8,  ALLOCATABLE  :: Prefactors_2b(:,:)
!          LOGICAL :: CI_PRD=.TRUE. ! Production option, where the MAP of each CiCjCkCl is stored explesitelly PRO Prd
          LOGICAL :: CI_PRD=.FalsE. ! Ignores Production option, MAP is not constructed
          end module
!c==================== END of Module for CI_Prod
        MODULE USR_PAR
!*****************************
!    08/01/2015  SK & STR
! USER DEFINED PARAMETERS 
! 20 DOUBLE / 10 INTEGER / 10 LOGICAL / 10 STRING(32)
!*****************************
        DOUBLE PRECISION :: DPAR1=0.d0
        DOUBLE PRECISION :: DPAR2=0.d0
        DOUBLE PRECISION :: DPAR3=0.d0
        DOUBLE PRECISION :: DPAR4=0.d0
        DOUBLE PRECISION :: DPAR5=0.d0
        DOUBLE PRECISION :: DPAR6=0.d0
        DOUBLE PRECISION :: DPAR7=0.d0
        DOUBLE PRECISION :: DPAR8=0.d0
        DOUBLE PRECISION :: DPAR9=0.d0
        DOUBLE PRECISION :: DPAR10=0.d0
        DOUBLE PRECISION :: DPAR11=0.d0
        DOUBLE PRECISION :: DPAR12=0.d0
        DOUBLE PRECISION :: DPAR13=0.d0
        DOUBLE PRECISION :: DPAR14=0.d0
        DOUBLE PRECISION :: DPAR15=0.d0
        DOUBLE PRECISION :: DPAR16=0.d0
        DOUBLE PRECISION :: DPAR17=0.d0
        DOUBLE PRECISION :: DPAR18=0.d0
        DOUBLE PRECISION :: DPAR19=0.d0
        DOUBLE PRECISION :: DPAR20=0.d0
        INTEGER :: IPAR1=0
        INTEGER :: IPAR2=0
        INTEGER :: IPAR3=0
        INTEGER :: IPAR4=0
        INTEGER :: IPAR5=0
        INTEGER :: IPAR6=0
        INTEGER :: IPAR7=0
        INTEGER :: IPAR8=0
        INTEGER :: IPAR9=0
        INTEGER :: IPAR10=0
        CHARACTER(LEN=32) :: SPAR1=''
        CHARACTER(LEN=32) :: SPAR2=''
        CHARACTER(LEN=32) :: SPAR3=''
        CHARACTER(LEN=32) :: SPAR4=''
        CHARACTER(LEN=32) :: SPAR5=''
        CHARACTER(LEN=32) :: SPAR6=''
        CHARACTER(LEN=32) :: SPAR7=''
        CHARACTER(LEN=32) :: SPAR8=''
        CHARACTER(LEN=32) :: SPAR9=''
        CHARACTER(LEN=32) :: SPAR10=''
        LOGICAL :: LPAR1=.FALSE.
        LOGICAL :: LPAR2=.FALSE.
        LOGICAL :: LPAR3=.FALSE.
        LOGICAL :: LPAR4=.FALSE.
        LOGICAL :: LPAR5=.FALSE.
        LOGICAL :: LPAR6=.FALSE.
        LOGICAL :: LPAR7=.FALSE.
        LOGICAL :: LPAR8=.FALSE.
        LOGICAL :: LPAR9=.FALSE.
        LOGICAL :: LPAR10=.FALSE.
        END MODULE  
