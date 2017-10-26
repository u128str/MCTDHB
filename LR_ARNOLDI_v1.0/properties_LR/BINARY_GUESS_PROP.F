c===========================================================================
       subroutine  Binary_guess_prop
       USE   CI_All
       USE   CI_prod
       USE   W_INTERPARTICLE
       USE   DVR_ALL
       USE   rR_hW
       USE   PROP_MB
       implicit NONE
       INTEGER ::  ierr,nstate,iter
       real*8   time_start,time_finish,tm1,tm2,enrg,Init_Start_Point_t
       logical :: Restart
c================ Read Restart files or their intialization.
       Init_Start_Point_t=1000000d0
       CALL Open_binary !Needed for restarts
       CALL Read_Header(778,'full')
       CALL Read_Header(777,'full')
       tm1=Init_Start_Point_t
       call Analyse_BINARY(778, tm1, 1000)
       tm2=Init_Start_Point_t
       call Analyse_BINARY(777, tm2, 1000)
      time_psi_MAX=tm1
      time_cic_MAX=tm2
       CALL  Close_binary 
       write(6,*)"init.in: Binary_Start_Point_t=",Init_Start_Point_t
       write(6,*)"in old   Binary_Start_Point_t=",Binary_Start_Point_t
       write(6,*)"PSI_bin: Binary_Start_Point_t=",tm1
       write(6,*)"CIc_bin: Binary_Start_Point_t=",tm2
       if(abs(tm1-tm2).le.100*epsilon(0d0)) then 
       write(6,*)"I CAN take as a Start_Point_t=",tm2
!     & ,Binary_Start_Point_t,tm1
!             if(abs(Init_Start_Point_t-tm1).le.100*epsilon(0d0)) then
!          call READ_init ! reread init.in to owerwrite the old namelists from binary file
!         Binary_Start_Point_t=Init_Start_Point_t
!      write(6,*)"I WORK with Binary_Start_Point_t=",Binary_Start_Point_t
            return
!             endif
       endif

c       if(abs(Init_Start_Point_t-tm2).ge.100*epsilon(0d0)) then
c       write(6,*)"Change the Binary_Start_Point_t=",tm2,"in init.in"
c       stop
c       endif
c       if(abs(Init_Start_Point_t-tm1).ge.100*epsilon(0d0)) then
c       write(6,*)"Change the Binary_Start_Point_t=",tm1,"in init.in"
c       stop
c       endif


       end SUBROUTINE Binary_guess_prop
