        subroutine Share_init_prop
        use SHARED_DIMS
        use CI_ALL
        use Parallel_Orb
        use Parallel_CI
        use W_INTERPARTICLE
        use DVR_ALL
        use rR_hW
        use CI_prod 
        implicit NONE
c================================== MPI ======================================
        logical :: ttt
        NAMELIST /PRM/ Morb, Npar, JOB_PreFac,xlambda_0,
     .  GUESS,Binary_Start_Point_t,
     . Time_Res_Orb_File_Name,Time_Res_CIc_File_Name,ORB_DIAG
        NAMELIST /ORB/ DIM_MCTDHB, NDX, NDY, NDZ,
     . Wxx_TYPE,Time_DVRMETHODX,Time_DVRMETHODY,Time_DVRMETHODZ,
     . Time_xint,Time_xfnl,Time_yint,Time_yfnl,Time_zint,Time_zfnl,
     . Time_mass,WTD
        NAMELIST /TM/
     & Time_Bgn,Time_Fnl,Time_MAX,Time_Print_Step,Time_ici_prt,Time_tau
     & ,Time_TolError_Total,Time_minsil,Time_maxsil,Time_Min_Occ_Allowed


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
          
        ORB_Binr=.TRUE.
        ORB_Rest=.False.
        ORB_DIAG=.False.
        end subroutine


