C=============== Main MPI File calls only two subprogramms
       program ParaMain
       USE   CI_All
       USE   CI_prod
       USE   W_INTERPARTICLE
       USE   DVR_ALL
       USE   rR_hW
       implicit NONE
c=================== MPI 
       INCLUDE 'mpif.h'
       INTEGER ::  MYID,numprocs,ierr,iter,i
       real*8   time_start,time_finish,tm1,tm2,enrg,Init_Start_Point_t
       real*4, dimension(2) :: TArray 
       real*4 tResult
        character(len=50) :: MCTDHB_V
        character*4 :: rtype ! can be 'scan' or 'full'
       logical :: Restart
c==========================================================
          MCTDHB_V="Current version 3.2.32 l Heidelberg Mar (2015)"
c================ MPI ===================================
         call MPI_INIT(ierr)
         call MPI_COMM_RANK(MPI_COMM_WORLD,MYID,ierr)
         call MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,ierr)
         write(6,*)" Process ",MYID," of ", numprocs
c==========================================================
      time_start=MPI_WTIME(ierr)
c==========================================================
c================Master node  Reads initial arrays namelists and so on...
       IF(MYID.eq.0) call READ_init
c================ Read Restart files or their intialization.
       IF(MYID.eq.0) call Binary_guess !CASE involving operations with binary files
c==========================================================
       call Share_init !Master node shares the namelist of the current job
       write(6,*) "Initial data are shared"
       call Init_CI
       write(6,*) "CI initialization is done"
       call Get_Operators
       write(6,*) "DVR initialization is done"
       call read_parse_V_W_Psi !parsing V_W_Psi_string.in file
       write(6,*) "File with user potential was parsed"
       call Get_InterPart(Time_Bgn)
       write(6,*)"Interparticle interaction was computed at T=",Time_Bgn
!          MCTDHB_V="Current version 3.2.14  w7h Heidelberg Aug (2014)"
          write(6,*)"MCTDHB_V=",MCTDHB_V
          write(6,*)banner(9)
          write(banner(9),'(a8,a50,a8)') banner(9),MCTDHB_V,"#"
          banner(43)=banner(9)
!          write(6,*)banner(9)
!          stop 'Version' 

       IF(MYID.eq.0) then 
        SELECT CASE (MB_JOB_TYPE)
        CASE ("ALL")
        CALL MASTER_MCTDHB
        CASE ("FCI")
        CALL MASTER_FCI
        CASE DEFAULT
        Write (6,*) "This option", MB_JOB_TYPE, " is not supported"
        stop
        END SELECT
       ENDIF
       IF(MYID.ne.0) CALL ASSIST_ALL

          time_finish=MPI_WTIME(ierr)
       write(6,*) 
     ."Proc",MYID,"MPI Execution time [sec]",time_finish-time_start
         call MPI_FINALIZE(ierr)
          IF(MYID.eq.0) then
          write(6,'(a66)')(banner0(i)//CHALL,i=1,7)
          write(6,'(a66)')(banner(i)//CHALL,i=7,10)
          endif
!          write(6,*)"OK1"
c================= CI Prod
           IF(ALLOCATED(Ind_CI_1b)) deallocate(Ind_CI_1b)
           IF(ALLOCATED(Prefactors_1b)) deallocate(Prefactors_1b)
           IF(ALLOCATED(Ind_CI_2b)) deallocate(Ind_CI_2b)
           IF(ALLOCATED(Prefactors_2b)) deallocate(Prefactors_2b)
c=================
           IF(ALLOCATED(MCNK)) deallocate(MCNK)
           IF(ALLOCATED(Nmax)) deallocate(Nmax)
           IF(ALLOCATED(WOP_X)) deallocate(WOP_X)
           IF(ALLOCATED(WOP_Y)) deallocate(WOP_Y)
           IF(ALLOCATED(WOP_Z)) deallocate(WOP_Z)
           IF(ALLOCATED(W2xx)) deallocate(W2xx)
           IF(ALLOCATED(Wxx)) deallocate(Wxx)
           IF(ALLOCATED(FIJ)) deallocate(FIJ)
           IF(ALLOCATED(F)) deallocate(F)
           IF(ALLOCATED(weight)) deallocate(weight)
           IF(ALLOCATED(weight_X)) deallocate(weight_X)
           IF(ALLOCATED(weight_Y)) deallocate(weight_Y)
           IF(ALLOCATED(weight_Z)) deallocate(weight_Z)
           IF(ALLOCATED(ORT_X)) deallocate(ORT_X)
           IF(ALLOCATED(ORT_Y)) deallocate(ORT_Y)
           IF(ALLOCATED(ORT_Z)) deallocate(ORT_Z)
           IF(ALLOCATED(Op_X)) deallocate(Op_X)
           IF(ALLOCATED(Op_Y)) deallocate(Op_Y)
           IF(ALLOCATED(Op_Z)) deallocate(Op_Z)
           IF(ALLOCATED(dif1matX)) deallocate(dif1matX)
           IF(ALLOCATED(dif1matY)) deallocate(dif1matY)
           IF(ALLOCATED(dif1matZ)) deallocate(dif1matZ)
           deallocate(ort_kx,ort_ky,ort_kz)
           deallocate(ZRIJ,HIJ,ZRIJ1,AllZRIJ,AllHIJ,InvZRIJ,AllZRIJ0)
           deallocate(ZRIJKL,WIJKL,ZRIJKL1)
           deallocate(Nocc,NatVec,NO_Expectation_x)
           deallocate(TERM_INDEX_1B,TERM_REQ_1B, RESCALE_1B)
           deallocate(TERM_INDEX_2B,TERM_REQ_2B, RESCALE_2B,ZMU)
       stop "COMPUTATION IS FINISHED ok! BBB: Be superB with the mctdhB"
      end program ParaMain
