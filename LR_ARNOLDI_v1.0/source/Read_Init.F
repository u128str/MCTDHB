        subroutine READ_init
        use SHARED_DIMS
        use CI_ALL
        use Parallel_Orb
        use Parallel_CI
        use W_INTERPARTICLE
        use DVR_ALL
        use rR_hW
        use CI_prod 
        USE USR_PAR
        logical :: ttt
        
        NAMELIST /PRM/ MB_JOB_TYPE, Morb, Npar, JOB_PreFac,xlambda_0,
     . GUESS,Binary_Start_Point_t,
     . Time_Res_Orb_File_Name,Time_Res_CIc_File_Name,ORB_DIAG
        NAMELIST /ORB/ DIM_MCTDHB, NDX, NDY, NDZ, 
     . Wxx_TYPE,Time_DVRMETHODX,Time_DVRMETHODY,Time_DVRMETHODZ,
     . Time_xint,Time_xfnl,Time_yint,Time_yfnl,Time_zint,Time_zfnl,
     . Time_mass, WTD
        NAMELIST /TM/
     & Time_Bgn,Time_Fnl,Time_MAX,Time_Print_Step,Time_ici_prt,Time_tau
     & ,Time_TolError_Total,Time_minsil,Time_maxsil
     & ,Time_intgr,Time_intgr_order,Time_intgr_max_step
     & ,Time_Min_Occ_Allowed
     & ,PRINT_DATA,print_level,TD_Err_scale
     & ,LZ,OMEGAZ,STATE
!STR add print_level 10 March 2015
!*****************************
!    08/01/2015  SK & STR
! USER DEFINED PARAMETERS 
! 20 DOUBLE / 10 INTEGER / 10 LOGICAL / 10 STRING(32)
!*****************************
        NAMELIST /USRPAR/ DPAR1,DPAR2,DPAR3,DPAR4,DPAR5,DPAR6,DPAR7
     . ,DPAR8,DPAR9,DPAR10,DPAR11,DPAR12,DPAR13,DPAR14,DPAR15,DPAR16
     . ,DPAR17,DPAR18,DPAR19,DPAR20
     . ,IPAR1,IPAR2,IPAR3,IPAR4,IPAR5,IPAR6,IPAR7,IPAR8,IPAR9,IPAR10 
     . ,SPAR1,SPAR2,SPAR3,SPAR4,SPAR5,SPAR6,SPAR7,SPAR8,SPAR9,SPAR10 
     . ,LPAR1,LPAR2,LPAR3,LPAR4,LPAR5,LPAR6,LPAR7,LPAR8,LPAR9,LPAR10 
c===============SOME default setup
        WTD=.FALSE.  !cAIS 15JUL2014
c===============SOME default setup
        Time_intgr='ABM'
        Time_intgr_order=7
        Intgr_total_steps=0
        Time_intgr_max_step=0.001d0
        LZ=.FALSE.
        OMEGAZ=0.d0
        STATE=1
        print_level=2
c===============END of  default setup
       
        write(6,*) "Trying to Read input.in file"
!        open(8,file="input.in",status='OLD',recl=80, delim='APOSTROPHE')
        open(8,file="input.in",status='OLD', delim='APOSTROPHE')
        read(8,NML=PRM) 
        read(8,NML=ORB) 
        read(8,NML=TM) 
        read(8,NML=USRPAR) 

c======== SOME EXCEPTIONS        
c======== SOME EXCEPTIONS        
c======== SOME EXCEPTIONS        
        IF ((REAL(job_prefac).lt.1.d-10).and.
     &     (ABS((DIMAG(job_prefac))-1.d0).lt.1.d-10).and.
     &     (STATE.ne.1)) then
           write(*,*) "I will not do a propagation with the 2nd
     &  Krylov Vector, Change STATE or Job_prefac in input.in!!"
           stop
        endif
        IF (((Time_DVRMETHODX.ne.4).or.
     &      (Time_DVRMETHODY.ne.4)).and.(LZ.eqv..TRUE.)) then
           write(*,*) "ANGULAR MOMENTUM WORKS WITH FFTDVR ONLY!!!"
           stop
        endif

        close(8)
        end subroutine


