C=============== Main MPI File calls only two subprogramms
       program LR_Arnoldi_Main
       USE   CI_All
       USE   CI_prod
       USE   W_INTERPARTICLE
       USE   DVR_ALL
       USE   rR_hW
       USE SHARED_DIMS
       USE LR_RAPHA
       USE LR_ARNOLDI_MOD
       USE PARALLEL_CI
       USE ARNOLDI_MOD

c       USE   PROP_MB 
       implicit NONE
c=================== MPI
       INCLUDE 'mpif.h' 
       INTEGER ::  iter, status_open, st
       INTEGER ::  MYID,numprocs,ierr,i,j,counter
!       INTEGER ::  ierr,i,j,counter


c       integer*4 :: ierr
c       real*4   start,finish, exec_time 
       real*8   time_start,time_finish,tm1,tm2,enrg,Init_Start_Point_t
       real*8 :: Time0,a,b
       real*4, dimension(2) :: TArray 
       real*4 tResult
       character*4 :: rtype ! can be 'scan' or 'full'
       logical :: Restart
       Complex*16 :: comp
       COMPLEX*16, allocatable :: VOUT(:) 
       COMPLEX*16, allocatable :: VIN(:) 
       REAL*8 :: time
       COMPLEX*16, allocatable :: PSI(:,:)
       COMPLEX*16, allocatable :: something(:)
c==========================================================
c================ MPI ===================================
         call MPI_INIT(ierr)
         call MPI_COMM_RANK(MPI_COMM_WORLD,MYID,ierr)
         call MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,ierr)
         write(6,*)" Process ",MYID," of ", numprocs
c=========================== Command line arguments ==============
c==========================================================
       call READ_init_LR_Arnoldi
       call READ_init
c======= Read Restart files or their intialization.
       call Binary_guess_prop !CASE involving operations with binary files
       call Share_init_prop !Master node shares the namelist of the current job
       call Init_CI ! Needed to activate Ter_ReQ for densitity matrix ....
!AIS 2017
        call read_parse_V_W_Psi !parsing V_W_Psi_string.in file
         write(6,*) "File with user potential was parsed" 
 

!       if(construction.eqv..FALSE.) then
     
!         call Get_Operators
!         call DIAG_ARNOLDI_MPI()

!       else  


       call Get_Operators

       xlambda0=xlambda_0
       Time0=time_finish
       call Get_InterPart(Time0)



       allocate(VIN_LR(Nconf)) 
       allocate(VIN(Nconf)) 
       allocate(InvZRIJ_LR(Morb,Morb)) 
       allocate(PSI_LR(NDX*NDY*NDZ,Morb))
       allocate(PSI(NDX*NDY*NDZ,Morb))
       if(myid.eq.0) CALL Arnoldi_preliminary
       IF(MYID.ne.0) CALL ASSIST_ALL



      if(myid.ne.0)  then 
         CALL MPI_BCAST(VIN_LR,Nconf,MPI_DOUBLE_COMPLEX,0
     .            ,MPI_COMM_WORLD,ierr)

         CALL MPI_BCAST(PSI_LR,NDX*NDY*NDZ*Morb,MPI_DOUBLE_COMPLEX,0
     .            ,MPI_COMM_WORLD,ierr)

         CALL MPI_BCAST(time_LR,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)

         CALL  MPI_BCAST(xlambda_0,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
      end if

       VIN=VIN_LR
       PSI=PSI_LR
       time=time_LR  

      if(myid.eq.0) write(*,*) "PSI, VIN and time available on all PEs"

!       if(construction.eqv..FALSE.) then
     
!         call Get_Operators
!         call DIAG_ARNOLDI_MPI(PSI,VIN)
      
!       else

!===== Get Everything for C_IJ production on every PE
       allocate(Prefactors_1b_Arnoldi(Nconf,MaxTrm1b), stat=ierr)
          if(ierr /= 0) then
          write(*,*)"allocation error for Prefactors_1b_Arnoldi"
          else
        if(ierr == 0)write(*,*)"allocation ok for Prefactors_2b_Arnoldi"
          endif
       allocate(Prefactors_2b_Arnoldi(Nconf,MaxTrm2b), stat=ierr)
          if(ierr /= 0) then
          write(*,*)"allocation error for Prefactors_2b_Arnoldi"
          else
        if(ierr == 0)write(*,*)"allocation ok for Prefactors_2b_Arnoldi"
          endif
       allocate(Ind_CI_1b_Arnoldi(Nconf,MaxTrm1b), stat=ierr)
          if(ierr /= 0) then
          write(*,*)"allocation error for Ind_CI_1b_Arnoldi"
          else
        if(ierr == 0)write(*,*)"allocation ok for Ind_CI_1b_Arnoldi"
          endif
       allocate(Ind_CI_2b_Arnoldi(Nconf,MaxTrm2b), stat=ierr)
          if(ierr /= 0) then
          write(*,*)"allocation error for Ind_CI_2b_Arnoldi"
          else
        if(ierr == 0)write(*,*)"allocation ok for Ind_CI_2b_Arnoldi"
          endif
      
       Ind_CI_1b_Arnoldi=1
       Ind_CI_2b_Arnoldi=1
       Prefactors_1b_Arnoldi=0d0 
       Prefactors_2b_Arnoldi=0d0 

 
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       write(*,*) myid,"STR 137 ALLOCATED",allocated(Ind_CI_1b)
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       write(*,*) myid,CI_Proc_From(myid+1),CI_Proc_TILL(myid+1)

       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       do i=1, numprocs
         write(*,*) myid,"MYID_Trm",MyID_TRM(i)
         write(*,*) myid,CI_Proc_From(i),CI_Proc_TILL(i)
       end do
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)




!===== Set CI arrays according to availability on each PE
       do i=1, numprocs
          if(myid.eq.i-1) then
 
            if(MyID_TRM(myid+1).eq.1) then
              counter=1 
              do j=CI_Proc_FROM(myid+1), CI_PROC_TILL(myid+1)
                  Prefactors_1b_Arnoldi(:,j)=
     .                 Prefactors_1b(:,counter)
                  IND_CI_1b_Arnoldi(:,j)=
     .                 IND_CI_1b(:,counter)

                  counter=counter+1 
              end do 
            end if

            if(MyID_TRM(myid+1).eq.2) then
              counter=1 
              do j=CI_Proc_FROM(myid+1), CI_PROC_TILL(myid+1)
                  Prefactors_2b_Arnoldi(:,j)=
     .                 Prefactors_2b(:,counter)
                  IND_CI_2b_Arnoldi(:,j)=
     .                 IND_CI_2b(:,counter)
                  counter=counter+1 
              end do 

            end if

            if(MyID_TRM(myid+1).eq.3) then
              counter=1 
              do j=CI_Proc_FROM(myid+1), MaxTrm1b
                  Prefactors_1b_Arnoldi(:,j)=
     .                 Prefactors_1b(:,counter)
                  IND_CI_1b_Arnoldi(:,j)=
     .                 IND_CI_1b(:,counter)

                  counter=counter+1 
              end do 

              counter=1 
              do j=1, CI_PROC_TILL(myid+1)
                  Prefactors_2b_Arnoldi(:,j)=
     .                 Prefactors_2b(:,counter)
                  IND_CI_2b_Arnoldi(:,j)=
     .                 IND_CI_2b(:,counter)

                  counter=counter+1 
              end do 

            end if
             
          end if 
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)

       end do


       write(*,*) myid, size(Ind_CI_1b) 
       write(*,*) myid, size(Ind_CI_2b)



!====  BROADCAST CI arrays all to all
        do i=1, numprocs
 
           if(myid_Trm(i).eq.1) then
             do j=CI_PROC_FROM(i), CI_PROC_TILL(i)

            call MPI_BCAST(Prefactors_1b_Arnoldi(:,j),Nconf,
     .          MPI_DOUBLE,i-1,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(Ind_CI_1b_Arnoldi(:,j),Nconf,
     .          MPI_INTEGER,i-1,MPI_COMM_WORLD,ierr)
             
             end do 

           end if  

           if(myid_Trm(i).eq.2) then
             do j=CI_PROC_FROM(i), CI_PROC_TILL(i)

            call MPI_BCAST(Prefactors_2b_Arnoldi(:,j),Nconf,
     .          MPI_DOUBLE,i-1,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(Ind_CI_2b_Arnoldi(:,j),Nconf,
     .          MPI_INTEGER,i-1,MPI_COMM_WORLD,ierr)
             
             end do 

           end if  

           if(myid_Trm(i).eq.3) then
             do j=CI_PROC_FROM(i), MaxTrm1b

            call MPI_BCAST(Prefactors_1b_Arnoldi(:,j),Nconf,
     .          MPI_DOUBLE,i-1,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(Ind_CI_1b_Arnoldi(:,j),Nconf,
     .          MPI_INTEGER,i-1,MPI_COMM_WORLD,ierr)
             
             end do 
             do j=1, CI_Proc_TILL(i)

            call MPI_BCAST(Prefactors_2b_Arnoldi(:,j),Nconf,
     .          MPI_DOUBLE,i-1,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(Ind_CI_2b_Arnoldi(:,j),Nconf,
     .          MPI_INTEGER,i-1,MPI_COMM_WORLD,ierr)
             
             end do 

           end if  
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)


        end do





!      READ Zmu
       call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
       if(myid.eq.0) then  
          open(unit=218, file='muij.out', 
     .      status='old', action='READ', iostat=status_open)
            do 
              read(218,*,iostat=st) i,j,a,b
!STR              write(6,*) "STR ==== read mu", i,j,a,b
              if(st.ne.0) exit
              Zmu(i,j)=Cmplx(a,b,kind=8)
            end do 
          close(218)
       end if
       call MPI_BARRIER(MPI_COMM_WORLD,ierr) 

!      BROADCAST Zmu
       do i=1, Morb 
         do j=1, Morb 
           CALL MPI_BCAST(Zmu(i,j),1,MPI_DOUBLE_COMPLEX,0,
     .             MPI_COMM_WORLD,ierr)
         end do
       end do


       call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
       if(myid.eq.0) then  
          open(unit=218, file='Energy.out', 
     .      status='old', action='READ', iostat=status_open)
              read(218,*) Energy
          close(218)
       end if
!       do i=1,numprocs
         CALL MPI_BCAST(Energy,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
         if(myid.eq.0) write(*,*) Energy
!       end do 
       call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
 

!      READ InvZRIJ
       call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
       if(myid.eq.0) then  
          open(unit=218, file='InvZRIJ.out', 
     .      status='old', action='READ', iostat=status_open)
            do 
              read(218,*,iostat=st) i,j,a,b
              if(st.ne.0) exit
              InvZRIJ(i,j)=Cmplx(a,b,kind=8)
            end do 
          close(218)
       end if
       call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
!      BROADCAST InvZRIJ
       do i=1, Morb 
         do j=1, Morb 
           CALL MPI_BCAST(InvZRIJ(i,j),1,MPI_DOUBLE_COMPLEX,0,
     .                      MPI_COMM_WORLD,ierr)
         end do
       end do
       write(*,*) myid,"InvZRIJ",InvZRIJ
       call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
!       call get_mu_test(Zmu)
!       !Zmu=Matmul(Zmu,-InvZRij)



!      READ ZRIJ
       call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
       if(myid.eq.0) then  
          open(unit=218, file='ZRIJ.out', 
     .      status='old', action='READ', iostat=status_open)
            do  
              read(218,*,iostat=st) i,a,b
              if(st.ne.0) exit
              ZRIJ(i)=Cmplx(a,b,kind=8)
            end do
          close(218)
       end if
!      BROADCAST ZRIJ
       CALL MPI_BCAST(i,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       do j=1, i
         CALL MPI_BCAST(ZRIJ(j),1,MPI_DOUBLE_COMPLEX,0,
     .                      MPI_COMM_WORLD,ierr)
       end do 
       call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
       write(*,*) myid,"ZRIJ",ZRIJ
       call MPI_BARRIER(MPI_COMM_WORLD,ierr) 


!      READ ZRIJKL
       call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
       if(myid.eq.0) then  
          open(unit=218, file='ZRIJKL.out', 
     .      status='old', action='READ', iostat=status_open)
            do   
              read(218,*,iostat=st) i,a,b
              if(st.ne.0) exit 
              ZRIJKL(i)=Cmplx(a,b,kind=8)
            end do
          close(218)
       end if
!      BROADCAST ZRIJKL
       CALL MPI_BCAST(i,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       do j=1, i
         CALL MPI_BCAST(ZRIJKL(j),1,MPI_DOUBLE_COMPLEX,0,
     .                      MPI_COMM_WORLD,ierr)
       end do 
       call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
       write(*,*) myid,"ZRIJKL",ZRIJKL
       call MPI_BARRIER(MPI_COMM_WORLD,ierr) 

    
!      READ WIJKL
       call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
       if(myid.eq.0) then  
          open(unit=218, file='WIJKL.out', 
     .      status='old', action='READ', iostat=status_open)
            do 
              read(218,*,iostat=st) i,a,b
              if(st.ne.0) exit
              WIJKL(i)=Cmplx(a,b,kind=8)
            end do
          close(218)
       end if
!      BROADCAST WIJKL
       CALL MPI_BCAST(i,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       do j=1, i
         CALL MPI_BCAST(WIJKL(j),1,MPI_DOUBLE_COMPLEX,0,
     .                      MPI_COMM_WORLD,ierr)
       end do 
       call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
       if(myid.eq.6) write(*,*) i,WIJKL
       call MPI_BARRIER(MPI_COMM_WORLD,ierr) 


       write(*,*) "STR 2" 

!      stop

       if(task.eq.2) then
     
!         call Get_Operators
         call DIAG_ARNOLDI_MPI(PSI,VIN)
      
       else

  
       write(*,*) "STR before LR_ARNOLDI" 
         call LR_ARNOLDI(time,PSI,VIN)

       end if 

       call MPI_FINALIZE(ierr) 
c==========================================================
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
      end program LR_Arnoldi_Main
