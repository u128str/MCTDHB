       subroutine Open_binary()
!  Opens the Binary files PSI_bin and CIc_bin
       integer :: unitnum1,unitnum2
       logical :: I_OPEN1,I_EXIST1,I_OPEN2,I_EXIST2
       character*7 :: filename1,filename2
!================== First stage check on existence ================
       filename1="CIc_bin"
       unitnum1=777
       filename2="PSI_bin"
       unitnum2=778

       INQUIRE (FILE=filename1, OPENED=I_OPEN1, EXIST=I_EXIST1,
     &  NUMBER=I_NUMBER) 
       INQUIRE (FILE=filename2, OPENED=I_OPEN2, EXIST=I_EXIST2, 
     &  NUMBER=I_NUMBER) 

        IF(I_EXIST1.eqv..TRUE.) THEN
        WRITE(6,*)" Binary CIc_bin already exist",I_OPEN1,I_EXIST1
        ENDIF
        IF(I_EXIST2.eqv..TRUE.) THEN
        WRITE(6,*)" Binary PSI_bin already exist",I_OPEN2,I_EXIST2
        ENDIF

      Open (unit=unitnum1, File=filename1, Status='unknown', 
     &  Form='unformatted', Access='sequential')
      write(6,*)" Binary file ",filename1," is properly opened"

      Open (unit=unitnum2, File=filename2, Status='unknown', 
     &  Form='unformatted', Access='sequential')
      write(6,*)" Binary file ",filename2," is properly opened"


      end SUBROUTINE Open_binary
!==================================================================
       subroutine Close_binary()
       close(777)
       close(778)
       write(6,*)" Binary files are properly closed"
       end SUBROUTINE Close_binary
!==================================================================
       
!==================================================================
       subroutine Write_Header(unitnum)
        use SHARED_DIMS
        use CI_ALL
        use Parallel_Orb
        use Parallel_CI
        use W_INTERPARTICLE
        use DVR_ALL
        use rR_hW
        use CI_prod 
!                    Writes Headers
       integer :: unitnum
       character*16 :: fname
       
!========== The first line of the header is namelist       
        NAMELIST /PRM/ Morb, Npar, JOB_PreFac,xlambda_0,GUESS,
     . Binary_Start_Point_t,
     . Time_Res_Orb_File_Name,Time_Res_CIc_File_Name,ORB_DIAG
        NAMELIST /ORB/ DIM_MCTDHB, NDX, NDY, NDZ, NCdim, NRdim, NIdim,
     . NLdim, Wxx_TYPE,Time_DVRMETHODX,Time_DVRMETHODY,Time_DVRMETHODZ,
     . Time_xint,Time_xfnl,Time_yint,Time_yfnl,Time_zint,Time_zfnl,
     . Time_mass, WTD
        NAMELIST /TM/
     & Time_Bgn,Time_Fnl,Time_MAX,Time_Print_Step,Time_ici_prt,Time_tau
     & ,Time_TolError_Total,Time_minsil,Time_maxsil,Time_Min_Occ_Allowed

        write(unitnum)  Morb, Npar, JOB_PreFac,xlambda_0,GUESS,
     . Binary_Start_Point_t,
     . Time_Res_Orb_File_Name,Time_Res_CIc_File_Name,ORB_DIAG
        write(unitnum) DIM_MCTDHB, NDX, NDY, NDZ, 
     . Wxx_TYPE,Time_DVRMETHODX,Time_DVRMETHODY,Time_DVRMETHODZ,
     . Time_xint,Time_xfnl,Time_yint,Time_yfnl,Time_zint,Time_zfnl,
     . Time_mass, WTD
        write(unitnum) 
     & Time_Bgn,Time_Fnl,Time_MAX,Time_Print_Step,Time_ici_prt,Time_tau
     & ,Time_TolError_Total,Time_minsil,Time_maxsil,Time_Min_Occ_Allowed
c===========================================================================
       if(unitnum.eq. 777) fname="OLD_Cic_input.in"
       if(unitnum.eq. 778) fname="OLD_PSI_input.in"
       write(6,*)"WRITES Headers are written to Binary file ",fname
c===========================================================================
c        stop "Header problem"
       end SUBROUTINE Write_Header
       
!==================================================================
       subroutine Read_Header(unitnum,rtype)
        use SHARED_DIMS
        use CI_ALL
        use Parallel_Orb
        use Parallel_CI
        use W_INTERPARTICLE
        use DVR_ALL
        use rR_hW
        use CI_prod 
!                    Reads  Headers
       integer :: unitnum
       character*16 :: fname
        character*4 :: rtype !'scan','full'
       
!========== The first line of the header is namelist       
        NAMELIST /PRM/ Morb, Npar, JOB_PreFac,xlambda_0,GUESS,
     . Binary_Start_Point_t,
     . Time_Res_Orb_File_Name,Time_Res_CIc_File_Name,ORB_DIAG
        NAMELIST /ORB/ DIM_MCTDHB, NDX, NDY, NDZ,
     . Wxx_TYPE,Time_DVRMETHODX,Time_DVRMETHODY,Time_DVRMETHODZ,
     . Time_xint,Time_xfnl,Time_yint,Time_yfnl,Time_zint,Time_zfnl,
     . Time_mass, WTD
        NAMELIST /TM/
     & Time_Bgn,Time_Fnl,Time_MAX,Time_Print_Step,Time_ici_prt,Time_tau
     & ,Time_TolError_Total,Time_minsil,Time_maxsil,Time_Min_Occ_Allowed
c===========================================================================
c        read(unitnum) PRM
c        read(unitnum) ORB 
c        read(unitnum) TM 
c       write(6,*)"READS Headers from Binary file 0"
          select case (rtype)
           case('full')
       read(unitnum,end=10,err=10)Morb,Npar,JOB_PreFac,xlambda_0,GUESS,
     . Binary_Start_Point_t,
     . Time_Res_Orb_File_Name,Time_Res_CIc_File_Name,ORB_DIAG
c       write(6,*)"READS Headers from Binary file 1"
        read(unitnum,end=10,err=10)DIM_MCTDHB,NDX,NDY,NDZ,
     . Wxx_TYPE,Time_DVRMETHODX,Time_DVRMETHODY,Time_DVRMETHODZ,
     . Time_xint,Time_xfnl,Time_yint,Time_yfnl,Time_zint,Time_zfnl,
     . Time_mass, WTD
c       write(6,*)"READS Headers from Binary file 2"
        read(unitnum,end=10,err=10) 
     & Time_Bgn,Time_Fnl,Time_MAX,Time_Print_Step,Time_ici_prt,Time_tau
     & ,Time_TolError_Total,Time_minsil,Time_maxsil,Time_Min_Occ_Allowed
c       write(6,*)"READS Headers from Binary file"
c===========================================================================
       if(unitnum.eq. 777) fname="OLD_Cic_input.in"
       if(unitnum.eq. 778) fname="OLD_PSI_input.in"
!       open(776,file=fname, status='unknown', recl=90,
!     &                   delim='APOSTROPHE')
       open(776,file=fname, status='unknown', 
     &                   delim='APOSTROPHE')
        write(776,NML=PRM)
        write(776,NML=ORB)
        write(776,NML=TM)
        close(776)
       write(6,*)"FULL Header from old binary file is in ",fname
       return
            case ('scan')
       read(unitnum,end=10,err=10) 
       read(unitnum,end=10,err=10)
       read(unitnum,end=10,err=10)
       write(6,*)"SCAn Headers are skiped ok",unitnum
            return
            end select
c===========================================================================
10      write(6,*)  "Binary file is corrupted or empty, change Guess"
        stop
       end SUBROUTINE Read_Header
!==================================================================
!==================================================================
       subroutine Write_PSI_string(time,iter,enrg,PSI)
        use SHARED_DIMS
        use CI_ALL
        use Parallel_Orb
        use Parallel_CI
        use W_INTERPARTICLE
        use DVR_ALL
        use rR_hW
        use CI_prod 
        real*8 ::  time,enrg
        integer :: iter
       COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb)  :: PSI
!                    Writes Headers
!========== Typical line of the   PSI file     
! time; iteration; energy;rhoij,rhoijkl
        write(778)  time, iter, energy
     &   ,ZRIJ,ZRIJKL,PSI,VTRAP_EXT
c        write(6,*) "TO binary ", time, iter, energy
c===========================================================================
       end SUBROUTINE Write_PSI_string
c===========================================================================
       subroutine Write_CIc_string(time,iter,enrg,VIN)
       USE   SHARED_DIMS
       USE   rR_hW
       USE   CI_All
       USE   W_INTERPARTICLE
       USE   DVR_ALL
       real*8 ::  time,enrg
       integer :: iter,i
       COMPLEX*16, DIMENSION(Nconf) :: VIN
!                    Writes Headers
!========== Typical line of the   PSI file     
! time; iteration; energy;rhoij,rhoijkl
!       write(6,*) "TObinary before ",VIN
c        write(777)  time, iter, energy,(VIN(i),i=1,Nconf)
        write(777)  time, iter, energy,VIN
!       write(6,*) Nconf,"TObinary", time, iter, energy,VIN(1),VIN(Nconf)
       end SUBROUTINE Write_CIc_string
!==================================================================
       subroutine Read_BIN_string(unitnum,time,iter,enrg,state)
        real*8 ::  time, enrg
        integer :: iter,state,unitnum
                             ! time; iteration; energy
        read(unitnum,end=10,err=20)  time, iter, enrg
         return
c===========================================================================
10       state=0
      write(6,*)"Last record From binary ",unitnum, time,iter,enrg,state
          return
20       state=-1
      write(6,*)"Wrong binary cannot read",
     . unitnum, time, iter, enrg, state
          return
       end SUBROUTINE Read_BIN_string
c===========================================================================
c===========================================================================
       subroutine Analyse_BINARY(unitnum,
     .                      required_time,required_iter) ! Call after the header is read 
        use SHARED_DIMS
        use CI_ALL
        use Parallel_Orb
        use Parallel_CI
        use W_INTERPARTICLE
        use DVR_ALL
        use rR_hW
        use CI_prod 
       integer :: unitnum,iter,required_iter
       character*16 :: fname
       real*8 :: tm,required_time,enrg
c       CALL Open_binary !Needed for restarts
c       CALL Read_Header(778)
        state=1
        DO 20 WHILE (state .GT.0)
!        write(6,*)" I am reading binary stuff before",state
        select case (unitnum)
        case (778)
        call Read_BIN_string(778,tm,iter,enrg,state)
        write(6,*)"From PSI_bin file",tm,iter,enrg
        case (777)
        call Read_BIN_string(777,tm,iter,enrg,state)
        write(6,*)"From CIc_bin file",tm,iter,enrg
!        call Read_Cic_string(tm,iter,state)
        end select 
        IF(ABS(tm-required_time).le.1.0d-8) state=-2
!        write(6,*)"From ",unitnum," binary file",tm,iter,enrg
 20     IF (state.le.0) write(6,*)
     &   "Binary reading was stopped on",tm,iter,enrg
         required_time=tm
!        write(6,*)"From PSI_bin file",tm,iter,enrg

!       CALL  Close_binary !Needed for restarts
       end SUBROUTINE Analyse_BINARY
c===========================================================================
       subroutine READ_BINARY(unitnum,
     .                      required_time,required_iter,VIN,PSI) ! Call after the header is read 
       use PASS_ARG
       USE SHARED_DIMS
       USE rR_hW
       USE CI_ALL
       USE Parallel_Orb
       USE DVR_ALL !needed by VTRAP_EXT
       COMPLEX*16, DIMENSION(Nconf) :: VIN
       COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb)  :: PSI
       integer :: unitnum,iter,required_iter,i,nstate
       character*16 :: fname
       real*8 :: tm,required_time,enrg
        nstate=1
        DO 20 WHILE (nstate .GT.0)
        select case (unitnum)
        case (778)
c        write(778)  time, iter, energy
c     &   ,ZRIJ,ZRIJKL,PSI,VTRAP_EXT
        read(778,end=11,err=11) tm,iter,enrg ,ZRIJ,ZRIJKL,PSI,VTRAP_EXT
!         write(6,*) "PSI reading at",tm, iter , enrg
        case (777)
!        read(777)  tm, iter, enrg ,(VIN(i),i=1,Nconf)
        read(777,end=12,err=12)  tm, iter, enrg ,VIN
c        write(6,*) "CI reading at",tm, iter , enrg, VIN(1),VIN(Nconf)
        end select 
!        IF(ABS(tm-required_time).le.10000*epsilon(0d0)) nstate=-2
        IF(ABS(tm-required_time).le.0.00001) nstate=-2
!         write(6,*) "PSI reading at",tm, required_time
 20     IF (nstate.le.0) write(6,'(a43,F10.5,a5,i4,a8,F20.16)')
     & "Guess from binary: reading was stopped at T=",tm,
     & "Iter=",iter,"Energy=",enrg
         required_time=tm
         energy=enrg
          return
11      write(6,*)"PSI_bin does not have records at time",
     . required_time
        stop 'Problem in Binary_IO'
12      write(6,*)"CIc_bin does not have records at time",
     . required_time
        stop 'Problem in Binary_IO'
       end SUBROUTINE READ_BINARY
c===========================================================================
       subroutine  Binary_guess
       USE   CI_All
       USE   CI_prod
       USE   W_INTERPARTICLE
       USE   DVR_ALL
       USE   rR_hW
       implicit NONE
       INTEGER ::  ierr,nstate,iter
       real*8   time_start,time_finish,tm1,tm2,enrg,Init_Start_Point_t
       logical :: Restart
c================ Read Restart files or their intialization.
       IF(GUESS=='BINR') THEN !CASE involving operations with binary files
       Init_Start_Point_t=Binary_Start_Point_t
       CALL Open_binary !Needed for restarts
       CALL Read_Header(778,'full')
       CALL Read_Header(777,'full')
       tm1=Init_Start_Point_t
       call Analyse_BINARY(778, tm1, 1000)
       tm2=Init_Start_Point_t
       call Analyse_BINARY(777, tm2, 1000)
       CALL  Close_binary 
       write(6,*)"init.in: Binary_Start_Point_t=",Init_Start_Point_t
       write(6,*)"in old   Binary_Start_Point_t=",Binary_Start_Point_t
       write(6,*)"PSI_bin: Binary_Start_Point_t=",tm1
       write(6,*)"CIc_bin: Binary_Start_Point_t=",tm2
       if(abs(tm1-tm2).le.1.0d-6) then 
       write(6,*)"I CAN take as a Start_Point_t=",tm2
!     & ,Binary_Start_Point_t,tm1
             if(abs(Init_Start_Point_t-tm1).le.1.0d-8) then
          call READ_init ! reread init.in to owerwrite the old namelists from binary file
         Binary_Start_Point_t=Init_Start_Point_t
      write(6,*)"I WORK with Binary_Start_Point_t=",Binary_Start_Point_t
            return
             endif
       endif

       if(abs(Init_Start_Point_t-tm2).ge.1.0d-8) then
       write(6,*)"Change the Binary_Start_Point_t=",tm2,"in init.in"
       stop
       endif
       if(abs(Init_Start_Point_t-tm1).ge.1.0d-8) then
       write(6,*)"Change the Binary_Start_Point_t=",tm1,"in init.in"
       stop
       endif

       ENDIF

       end SUBROUTINE Binary_guess
