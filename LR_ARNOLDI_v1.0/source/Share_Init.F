        subroutine Share_init
        use SHARED_DIMS
        use CI_ALL
        use Parallel_Orb
        use Parallel_CI
        use W_INTERPARTICLE
        use DVR_ALL
        use rR_hW
        use CI_prod 
        USE USR_PAR
        implicit NONE
c================================== MPI ======================================
       INCLUDE 'mpif.h'
       INTEGER ::  ierr,MYID,numprocs
        logical :: ttt
c        logical :: CI_ALL_init
c       logical :: W_INTERPARTICLE_init,DVR_ALL_init,rR_hW_init
c       external CI_ALL_init,W_INTERPARTICLE_init,DVR_ALL_init,rR_hW_init
        
!         ALL these constants are defined in the
!         MODULES_ALL_allocate.f90
        NAMELIST /PRM/ Morb, Npar, JOB_PreFac,xlambda_0,
     .  GUESS,Binary_Start_Point_t,
     . Time_Res_Orb_File_Name,Time_Res_CIc_File_Name,ORB_DIAG
!         ALL these constants are defined in the
!         MODULES_ALL_allocate.f90
        NAMELIST /ORB/ DIM_MCTDHB, NDX, NDY, NDZ, 
     . Wxx_TYPE,Time_DVRMETHODX,Time_DVRMETHODY,Time_DVRMETHODZ,
     . Time_xint,Time_xfnl,Time_yint,Time_yfnl,Time_zint,Time_zfnl,
     . Time_mass, WTD  !AIS 15Jul2014
!         ALL these constants are defined in the
!         MODULES_ALL_allocate.f90
        NAMELIST /TM/
     & Time_Bgn,Time_Fnl,Time_MAX,Time_Print_Step,Time_ici_prt,Time_tau
     & ,Time_TolError_Total,Time_minsil,Time_maxsil
     & ,Time_intgr,Time_intgr_order,Time_intgr_max_step
     & ,Time_Min_Occ_Allowed
     & ,PRINT_DATA,print_level,TD_Err_scale
     & ,LZ,OMEGAZ
!STR add print_level 10 March 2015
c===============SOME default setup
c===============SOME default setup
c=============================================================================

       call MPI_COMM_RANK(MPI_COMM_WORLD,MYID,ierr)
       call MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,ierr)
c       write(6,*)"Get_r_R Process ",MYID," of ", numprocs
c=============================================================================
c========================= MPI PART ======================================
c========================== NAMELIST /PRM/ Morb, Npar, JOB_PreFac,xlambda_0,Time_Res_Orb_File_Name,Time_Res_CIc_File_Name,ORB_DIAG
       call MPI_BCAST(Morb,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(Npar,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(JOB_PreFac,1,
     &                        MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(xlambda_0,1,MPI_DOUBLE_PRECISION,
     &                                            0,MPI_COMM_WORLD,ierr)
          xlambda0=xlambda_0
       call MPI_BCAST(Guess,4,MPI_CHARACTER, 0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(Binary_Start_Point_t,1,MPI_DOUBLE_PRECISION,
     &                                            0,MPI_COMM_WORLD,ierr)
c       call MPI_BCAST(ORB_Binr,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(Time_Res_Orb_File_Name,18,MPI_CHARACTER,
     $  0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(Time_Res_CIc_File_Name,18,MPI_CHARACTER,
     $  0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(ORB_DIAG,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)

c===========NAMELIST /ORB/ DIM_MCTDHB, NDX, NDY, NDZ, NCdim, NRdim, NIdim,NLdim, Wxx_TYPE,Time_DVRMETHOD,Time_xint,Time_xfnl,Time_mass,WTD
      call MPI_BCAST(DIM_MCTDHB,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(NDX,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(NDY,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(NDZ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(Wxx_TYPE,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(Time_DVRMETHODX,1,
     . MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(Time_DVRMETHODY,1,
     . MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(Time_DVRMETHODZ,1,
     . MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(Time_xint,1,MPI_DOUBLE_PRECISION,
     &                                            0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(Time_xfnl,1,MPI_DOUBLE_PRECISION,
     &                                            0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(Time_yint,1,MPI_DOUBLE_PRECISION,
     &                                            0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(Time_yfnl,1,MPI_DOUBLE_PRECISION,
     &                                            0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(Time_zint,1,MPI_DOUBLE_PRECISION,
     &                                            0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(Time_zfnl,1,MPI_DOUBLE_PRECISION,
     &                                            0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(Time_mass,1,MPI_DOUBLE_PRECISION,
     &                                            0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(WTD,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr) !c AIS 15Jul2014 Time -dependent W(r,r',t)

c=========NAMELIST /TM/ & Time_Bgn,Time_Fnl,Time_MAX,Time_Print_Step,Time_ici_prt,Time_tau,Time_TolError_Total,Time_minsil,Time_maxsil,Time_Min_Occ_Allowed
       call MPI_BCAST(Time_Bgn,1,MPI_DOUBLE_PRECISION,
     &                                            0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(Time_Fnl,1,MPI_DOUBLE_PRECISION,
     &                                            0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(Time_MAX,1,MPI_DOUBLE_PRECISION,
     &                                            0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(Time_Print_Step,1,MPI_DOUBLE_PRECISION,
     &                                            0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(Time_ici_prt,1,MPI_DOUBLE_PRECISION,
     &                                            0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(Time_tau,1,MPI_DOUBLE_PRECISION,
     &                                            0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(Time_TolError_Total,1,MPI_DOUBLE_PRECISION,
     &                                            0,MPI_COMM_WORLD,ierr)

       call MPI_BCAST(Time_minsil,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(Time_maxsil,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

       call MPI_BCAST(Time_intgr,18,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(Time_intgr_order,1,MPI_INTEGER,
     .                                            0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(Intgr_total_steps,1,MPI_INTEGER,
     .                                            0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(Time_Intgr_MAX_Step,1,MPI_DOUBLE_PRECISION,
     &                                            0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(Time_Min_Occ_Allowed,1,MPI_DOUBLE_PRECISION,
     &                                           0,MPI_COMM_WORLD,ierr)

       call MPI_BCAST(PRINT_DATA,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(print_level,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(TD_Err_scale,1,MPI_DOUBLE_PRECISION,
     &                                           0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(LZ,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(OMEGAZ,1,MPI_DOUBLE_PRECISION,
     &                                           0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(STATE,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

c=========NAMELIST /USRPAR/
c!*****************************
c!    08/01/2015  SK & STR
c! USER DEFINED PARAMETERS 
c! 20 DOUBLE / 10 INTEGER / 10 LOGICAL / 10 STRING(32)
c!*****************************
c        NAMELIST /USRPAR/ DPAR1,DPAR2,DPAR3,DPAR4,DPAR5,DPAR6,DPAR7
c     . ,DPAR8,DPAR9,DPAR10,DPAR11,DPAR12,DPAR13,DPAR14,DPAR15,DPAR16
c     . ,DPAR17,DPAR18,DPAR19,DPAR20
c     . ,IPAR1,IPAR2,IPAR3,IPAR4,IPAR5,IPAR6,IPAR7,IPAR8,IPAR9,IPAR10 
c     . ,SPAR1,SPAR2,SPAR3,SPAR4,SPAR5,SPAR6,SPAR7,SPAR8,SPAR9,SPAR10 
c     . ,LPAR1,LPAR2,LPAR3,LPAR4,LPAR5,LPAR6,LPAR7,LPAR8,LPAR9,LPAR10 
c
       call MPI_BCAST(DPAR1,1,MPI_DOUBLE_PRECISION,
     &                                            0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(DPAR2,1,MPI_DOUBLE_PRECISION,
     &                                            0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(DPAR3,1,MPI_DOUBLE_PRECISION,
     &                                            0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(DPAR4,1,MPI_DOUBLE_PRECISION,
     &                                            0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(DPAR5,1,MPI_DOUBLE_PRECISION,
     &                                            0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(DPAR6,1,MPI_DOUBLE_PRECISION,
     &                                            0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(DPAR7,1,MPI_DOUBLE_PRECISION,
     &                                            0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(DPAR8,1,MPI_DOUBLE_PRECISION,
     &                                            0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(DPAR9,1,MPI_DOUBLE_PRECISION,
     &                                            0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(DPAR10,1,MPI_DOUBLE_PRECISION,
     &                                            0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(DPAR11,1,MPI_DOUBLE_PRECISION,
     &                                            0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(DPAR12,1,MPI_DOUBLE_PRECISION,
     &                                            0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(DPAR13,1,MPI_DOUBLE_PRECISION,
     &                                            0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(DPAR14,1,MPI_DOUBLE_PRECISION,
     &                                            0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(DPAR15,1,MPI_DOUBLE_PRECISION,
     &                                            0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(DPAR16,1,MPI_DOUBLE_PRECISION,
     &                                            0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(DPAR17,1,MPI_DOUBLE_PRECISION,
     &                                            0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(DPAR18,1,MPI_DOUBLE_PRECISION,
     &                                            0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(DPAR19,1,MPI_DOUBLE_PRECISION,
     &                                            0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(DPAR20,1,MPI_DOUBLE_PRECISION,
     &                                            0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(IPAR1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(IPAR2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(IPAR3,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(IPAR4,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(IPAR5,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(IPAR6,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(IPAR7,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(IPAR8,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(IPAR9,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(IPAR10,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(SPAR1,32,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr) 
       call MPI_BCAST(SPAR2,32,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr) 
       call MPI_BCAST(SPAR3,32,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr) 
       call MPI_BCAST(SPAR4,32,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr) 
       call MPI_BCAST(SPAR5,32,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr) 
       call MPI_BCAST(SPAR6,32,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr) 
       call MPI_BCAST(SPAR7,32,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr) 
       call MPI_BCAST(SPAR8,32,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr) 
       call MPI_BCAST(SPAR9,32,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr) 
       call MPI_BCAST(SPAR10,32,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr) 
       call MPI_BCAST(LPAR1,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(LPAR2,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(LPAR3,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(LPAR4,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(LPAR5,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(LPAR6,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(LPAR7,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(LPAR8,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(LPAR9,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(LPAR10,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)


        ttt=CI_ALL_init()
        ttt=W_INTERPARTICLE_init()
        ttt=DVR_ALL_init()
        ttt=rR_hW_init()
        write(6,*)"Namelist xlambda_0",xlambda_0
        write(6,*)"GUESS TYPE=",GUESS
c==========================================================
                 select case (Guess)
             CASE ('DATA', 'HAND','BINR')   
         write(6,*)" Guess mode is OK"
        SELECT CASE (GUESS)
             CASE ('DATA')   
         IF(Time_Res_Cic_File_Name/='') CIc_Rest=.True. 
         if(Time_Res_Cic_File_Name=='') then
         stop "CIc formated Data file has wrong name "
         ENDif
        write(6,*)"CI coeff Guess is from ",Time_Res_Cic_File_Name
        write(6,*)"Orbitals Guess is from ",Time_Res_Orb_File_Name
         IF(Time_Res_Orb_File_Name/='') ORB_Rest=.True. 
         IF(Time_Res_Orb_File_Name=='')   then
         stop "PSI formated Data file has wrong name "
         ENDif
        ORB_DIAG=.False.
        ORB_Binr=.False.
             CASE ('HAND')   
        write(6,*)"CI coeff and orbitals  are given by HANDs"
        write(6,*)"CI coeff Guess is from user_guesslib/Guess_CI.F"
        write(6,*)"Orbitals Guess is from user_guesslib/Guess_PSI.F"
        IF(ORB_DIAG) THEN
      write(6,*)"Orbitals are by diagonalizing h with trap potential!"
      write(6,*)"see user_guesslib/Guess_PSI.F:"
      write(6,*)" IF(ORB_DIAG) CALL Guess_Diag(PSI,Vexternal)"
        ENDIF
        ORB_Binr=.False.
        ORB_Rest=.False.
             CASE ('BINR')   
        write(6,*)"CI coeff Guess is from CIc_bin"
        write(6,*)"Orbitals Guess is from PSI_bin"
        ORB_Binr=.TRUE.
        write(6,*)"The time-point in binary files , where we take"
        write(6,*)"initial PSI and T_in=",Binary_Start_Point_t
        ORB_Rest=.False.
        ORB_DIAG=.False.
        END SELECT 

               case default     ! something else
       WRITE(6,*)'Something else found in char*4 Guess: "', Guess, '"'
                  stop "Problem with guess"
                 end  select
        end subroutine


