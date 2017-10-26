        subroutine READ_init_prop
        use Prop_MB
        implicit none
        integer :: nargs,ist,j,k
        CHARACTER(LEN=255) :: string,string0,srting1  ! for command line arguments
        CHARACTER(LEN=255) :: message           ! for I/O error messages


!========================================================================        
         nargs=command_argument_count()
         write(6,*)" Number of command line argumenst", nargs

       if(nargs>0)then
!    read command line arguments to CMDL name list
         CALL GET_COMMAND(string)
         CALL GET_COMMAND_ARGUMENT(1,string0)

         j=index(string,"T_cmd",BACK = .TRUE.)
         k=index(string,"-h",BACK = .TRUE.)
         write(6,*)'T_cmd at',j
         if((j.eq.0).and.(k.eq.0)) then
         write(6,*)'ERROR in reading command line: ',string0
         write(6,*)' I will Try to read properties.in file'
         write(6,*)"prog  T_cmd=0.2d0 RhoX_cmd=.t. RhoK_cmd=.t. "
         write(6,*)"      g1x_cmd=.t. g1k_cmd=.t.  LR_cmd=.t. "
         write(6,*)"      Mnkvsk_cmd=.t. TimeALL_cmd=.t."
         write(6,*)"prog T_cmd=0.0, AnlsALL_cmd=.t." 
         write(6,*)" to produce ALL the data except LR at all the times"
         write(6,*)"!!! prog T_cmd=0.0 must be as the first argument"
         nargs=0 
         goto 11
         endif
         string0=string(j:len(string))
         string='&CMDL '//trim(string0)//' /' 
         write(6,*)'Full string is :',string
         write(6,*)'Parameters are:',string0
         read(string,NML=CMDL,iostat=ist)
         if((ist.ne.0).or.(k.ne.0)) then
         write(6,*)'ERROR in reading command line: ',ist
         write(6,*)'usage:' 
         write(6,*)"prog  T_cmd=0.2d0 RhoX_cmd=.t. RhoK_cmd=.t. "
         write(6,*)"      g1x_cmd=.t. g1k_cmd=.t.  LR_cmd=.t. "
         write(6,*)"      Mnkvsk_cmd=.t. TimeALL_cmd=.t."
         write(6,*)"prog T_cmd=0.0, AnlsALL_cmd=.t." 
         write(6,*)" to produce ALL the data except LR at all the times"
         write(6,*)"!!! prog T_cmd=0.0 must be as the first argument"
         stop 'PROPERTIES: problem with command line arguments'
         endif
! HERE is the   NAMELIST /CMDL/& T_cmd,& RhoX_cmd,& RhoK_cmd,& g1x_cmd,& g1k_cmd,& LR_cmd
        
         WRITE(*,NML=CMDL)                      
          
          T_From=T_cmd
          T_TILL=T_cmd
          T_points=1

          if(RhoX_cmd.eqv..TRUE.) then
          DATA_PSI=.True.  ! To create DATA/orb_R/time.dat
          DATA_CIc=.False.  !To create DATA/CIcnf/conf.dat
          NO_x=.True.  !To create DATA/orb_R/x-density.dat
          endif
          if(RhoK_cmd.eqv..TRUE.) then
          DATA_PSI=.True.  !To create DATA/orb_R/time.dat
          DATA_CIc=.False.  !To create DATA/CIcnf/conf.dat
!          DNS_XX=.True.  !To create DATA/g1_RR/x-correlations.dat
          NO_k=.True.    !To create DATA/orb_K/k-density.dat
          endif
          if(g1x_cmd.eqv..TRUE.) then
          DATA_PSI=.True.  !To create DATA/orb_R/time.dat
          DATA_CIc=.False.  !To create DATA/CIcnf/conf.dat
          DNS_XX=.True.  !To create DATA/g1_RR/x-correlations.dat
          endif
          if(g1k_cmd.eqv..TRUE.) then
          DATA_PSI=.True.  !To create DATA/orb_R/time.dat
          DATA_CIc=.False.  !To create DATA/CIcnf/conf.dat
          NO_k=.True.    !To create DATA/orb_K/k-density.dat
          DNS_KK=.True.  !To create DATA/g1_KK/k-correlations.dat
          endif
          if(LR_cmd.eqv..TRUE.) then
          DATA_PSI=.True.  !To create DATA/orb_R/time.dat
          DATA_CIc=.True.  !To create DATA/CIcnf/conf.dat
          get_LR=.True.  !To activate Linear Response computation at time=T_cmd 
          LR_maxsil=230  !Number of LR roots ....
          endif
          if(Mnkvsk_cmd.eqv..TRUE.) then  !Creates a single  DATA/orb_R/DNS_Minkovskii.dat to plot Space.vs.Time evolution of the density
          DATA_PSI=.False.  !To create DATA/orb_R/time.dat
          DATA_CIc=.False.  !To create DATA/CIcnf/conf.dat
          T_From=-10000.0d0
          T_TILL=+10000.0d0
          T_points=1000000
          endif
          if(AnlsALL_cmd.eqv..TRUE.) then  !Creates all DATA to plot
          Mnkvsk_cmd=.TRUE.
          DATA_PSI=.True.  !To create DATA/orb_R/time.dat
          DATA_CIc=.False.  !To create DATA/CIcnf/conf.dat
          NO_x=.True.      !To create DATA/orb_R/x-density.dat
          NO_k=.True.      !To create DATA/orb_K/k-density.dat
          DNS_XX=.True.    !To create DATA/g1_RR/x-correlations.dat
          DNS_KK=.True.    !To create DATA/g1_KK/k-correlations.dat
          T_From=-10000.0d0
          T_TILL=+10000.0d0
          T_points=1000000
          endif
          if(TimeALL_cmd.eqv..TRUE.) then  ! Selects all the time-points STR 21 Jan 2015
          T_From=-10000.0d0
          T_TILL=+10000.0d0
          T_points=1000000
          endif


         end if


!========================================================================        
11       continue
       if(nargs==0)then
        write(6,*) "Trying to Read properties.in file"
        open(8,file="properties.in",
     &         status='OLD',recl=80, delim='APOSTROPHE')
        read(8,NML=ZERO_body) 
        read(8,NML= ONE_body) 
        read(8,NML= TWO_body) 
!        write(*,*) PROJ_X
        read(8,NML= TWO_D)
!        write(*,*) PROJ_X
        read(8,NML= LR)
        read(8,NML= WSL)
        write(6,*) "Read of properties.in file is OK"
        close(8)
        endif

        end subroutine  READ_init_prop


