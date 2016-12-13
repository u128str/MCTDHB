         MODULE PROP_MB
         IMPLICIT NONE
         SAVE 

!New NAME LIST to read from command line and then to fill following nemelists...
!This is active only when commad line argumens are provided
! prog   T_cmf=0d0, RhoX_cmd=.t. LR_cmd=.t.
        REAL*8  ::  T_cmd     = 0.0d0
        REAL*8  ::  T_from_cmd     = 0.0d0
        REAL*8  ::  T_till_cmd     = 0.0d0
        REAL*8  ::  T_step_cmd     = 0.0d0
        LOGICAL,      public ::  RhoX_cmd = .FALSE.
        LOGICAL,      public ::  RhoK_cmd = .FALSE.
        LOGICAL,      public ::  g1x_cmd = .FALSE.
        LOGICAL,      public ::  g1k_cmd = .FALSE.
        LOGICAL,      public ::  LR_cmd = .FALSE.
        LOGICAL,      public ::  Mnkvsk_cmd = .FALSE.
        LOGICAL,      public ::  AnlsALL_cmd = .FALSE.
        LOGICAL,      public ::  TimeALL_cmd = .FALSE. !STR 21. Jan 2015 to platt all-times data
        NAMELIST /CMDL/&
        T_cmd,&
        RhoX_cmd,&
        RhoK_cmd,&
        g1x_cmd,&
        g1k_cmd,&
        LR_cmd,&
        Mnkvsk_cmd,&
        AnlsALL_cmd,&
        TimeALL_cmd,&
        T_from_cmd,&
        T_till_cmd,&
        T_step_cmd
        
! OLD tails (shift to main program)
        REAL*8  ::  time_psi_MAX, time_cic_MAX 

        LOGICAL,      public ::  Total_energy = .FALSE.
        REAL*8 ,      public ::    Kdip       = 0.d0
        LOGICAL,      public ::  DATA_PSI = .FALSE.
        LOGICAL,      public ::  DATA_CIc = .FALSE.

        REAL*8  ::  T_From     = 0.d0
        REAL*8  ::  T_TILL     = 0.d0
        integer ::  T_points     = 1

       NAMELIST /ZERO_body/& 
       Total_energy,&
       Kdip,&
       DATA_PSI,&
       DATA_CIc,&
       T_From,&
       T_TILL,&
       T_points

        LOGICAL,      public ::  NO_x  = .FALSE.
        LOGICAL,      public ::  NO_k  = .FALSE.
        LOGICAL,      public ::  DiLTN = .FALSE.
        LOGICAL,      public :: Pnot   = .FALSE.
        REAL*8,       public :: xstart = 0.d0
        REAL*8,       public :: xend   = 0.d0

        NAMELIST /ONE_body/ NO_x, No_k,DiLTN,Pnot, xstart, xend

        LOGICAL,      public ::  DNS_XX= .FALSE.
        LOGICAL,      public ::  DNS_KK= .FALSE.
        LOGICAL,      public :: corr1restr= .FALSE.
        LOGICAL,      public :: corr2restr= .FALSE.
        LOGICAL,      public :: corr1restrmom= .FALSE.
        LOGICAL,      public :: corr2restrmom= .FALSE.
        REAL*8,       public :: xini1= 0.d0
        REAL*8,       public :: xfin1= 0.d0
        REAL*8,       public :: xini2= 0.d0
        REAL*8,       public :: xfin2= 0.d0
        REAL*8,       public :: kxini1= 0.d0
        REAL*8,       public :: kxfin1= 0.d0
        REAL*8,       public :: kxini2= 0.d0
        REAL*8,       public ::  kxfin2= 0.d0 
        integer,     public :: xpts2 =0
        integer,     public ::kpts2  =0
        integer,     public :: xpts1 =0
        integer,     public ::kpts1  =0
        LOGICAL,      public :: lossops= .FALSE.
        REAL*8 ,      public :: border =0.d0

       NAMELIST /TWO_body/ DNS_XX,DNS_KK,&
                 corr1restr,xini1,xfin1,xpts1,&
                 corr1restrmom,kxini1,kxfin1,kpts1,&
                 corr2restr,xini2,xfin2,xpts2,&
                 corr2restrmom,kxini2,kxfin2,kpts2,&
                 lossops,border

         REAL*8,       public :: x1slice= 0.d0
         REAL*8,       public :: y1slice= 0.d0
         REAL*8,       public :: x2slice= 0.d0
         REAL*8,       public :: y2slice= 0.d0        
         LOGICAL,      public :: x1const= .FALSE.
         LOGICAL,      public ::y1const = .FALSE.
         LOGICAL,      public ::x2const = .FALSE.
         LOGICAL,      public ::y2const = .FALSE.
         LOGICAL,      public :: REALSPACE2D= .FALSE.
         LOGICAL,      public ::MOMSPACE2D  = .FALSE.
         LOGICAL,      public :: ZeroPadding2D  = .FALSE.
         INTEGER,      public :: Dilation2D     = 1
         LOGICAL,      public :: PROJ_X = .FALSE.
         Character(len=1) , public :: DIR='X'
         LOGICAL,      public :: L_Z = .FALSE.

        NAMELIST /TWO_D/MOMSPACE2D,REALSPACE2D,&
                       x1const,x1slice,&
                       y1const,y1slice,&
                       x2const,x2slice,&
                       y2const,y2slice,&
                       ZeroPadding2D,Dilation2D,&
                       PROJ_X,DIR,L_Z

         LOGICAL,      public :: get_LR
         INTEGER,      public :: LR_maxsil=100
         INTEGER,      public :: LR_maxroot=100

        NAMELIST /LR/get_LR,LR_maxsil,LR_maxroot

        LOGICAL,      public :: get_WSL
        NAMELIST /WSL/get_WSL



        INTEGER, public ::  MYID,numprocs  

        END MODULE PROP_MB
