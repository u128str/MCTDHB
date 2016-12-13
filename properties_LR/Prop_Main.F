C=============== Main MPI File calls only two subprogramms
       program Prop_Main
       USE   CI_All
       USE   CI_prod
       USE   W_INTERPARTICLE
       USE   DVR_ALL
       USE   rR_hW
       implicit NONE
c=================== MPI
       INCLUDE 'mpif.h' 
       INTEGER ::  iter
       INTEGER ::  MYID,numprocs,ierr,i

c       integer*4 :: ierr
c       real*4   start,finish, exec_time 
       real*8   time_start,time_finish,tm1,tm2,enrg,Init_Start_Point_t
       real*8 :: Time0
       real*4, dimension(2) :: TArray 
       real*4 tResult
       character*4 :: rtype ! can be 'scan' or 'full'
       logical :: Restart
c==========================================================
c================ MPI ===================================
         call MPI_INIT(ierr)
         call MPI_COMM_RANK(MPI_COMM_WORLD,MYID,ierr)
         call MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,ierr)
         write(6,*)" Process ",MYID," of ", numprocs
c=========================== Command line arguments ==============
c==========================================================
       call READ_init_prop
c======= Read Restart files or their intialization.
       call Binary_guess_prop !CASE involving operations with binary files
       call Share_init_prop !Master node shares the namelist of the current job
       call Init_CI ! Needed to activate Ter_ReQ for densitity matrix ....
       call Get_Operators

       Time0=time_finish
       call Get_InterPart(Time0)

       CALL Properties

       call MPI_FINALIZE(ierr) 
c==========================================================
         stop
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
           deallocate(ZRIJ,HIJ,ZRIJ1,AllZRIJ,AllHIJ,InvZRIJ,AllZRIJ0)
           deallocate(ZRIJKL,WIJKL,ZRIJKL1)
           deallocate(Nocc,NatVec,NO_Expectation_x)
           deallocate(TERM_INDEX_1B,TERM_REQ_1B, RESCALE_1B)
           deallocate(TERM_INDEX_2B,TERM_REQ_2B, RESCALE_2B,ZMU)
c           IF(ALLOCATED()) deallocate()
        stop 
      end program Prop_Main
