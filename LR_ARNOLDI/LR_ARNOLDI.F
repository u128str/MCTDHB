       subroutine LR_ARNOLDI(time,PSI,VIN) 
       USE ARNOLDI_MOD
       USE LR_ARNOLDI_MOD
       USE W_INTERPARTICLE
       USE PASS_ARG
       USE SHARED_DIMS
       USE rR_hW
       USE CI_ALL
       USE DVR_ALL
       USE LR_RAPHA
       USE CI_Prod
      USE omp_lib
      IMPLICIT NONE
      include 'mpif.h'
c      include 'debug.h'
c      include 'stat.h'
c 
c     %---------------%
c     | MPI INTERFACE |
c     %---------------%
 
      integer           comm, nprocs, final, nloc, done
c     %-----------------------------%
c     | Define maximum dimensions   |
c     | for all arrays.             |
c     | MAXN:   Maximum dimension   |
c     |         of the A allowed.   |
c     | MAXNEV: Maximum NEV allowed |
c     | MAXNCV: Maximum NCV allowed |
c     %-----------------------------%
c
      integer*8           ldv, int1, int2, maxdim_prec
      parameter         ( 
     &                    maxdim_prec=500000)
c
c     %--------------%
c     | Local Arrays |
c     %--------------%
c
      integer           iparam(11), ipntr(14), mjob(6)
      logical           select(maxncv)
      logical,allocatable ::        sel(:)
      logical           wantt,wantz,tempp_log 
      Complex*16
     &                  ax(maxn), d(maxncv), 
     &                  workd(3*maxn), 
     &                  workev(3*maxncv), resid(maxn), 
     &                  workl(3*maxncv*maxncv+5*maxncv),
     &                  vec(maxn), alpha, beta 
      Complex*16, allocatable ::  v(:,:)
      Complex*16, allocatable ::  matrix(:,:)
      Complex*16, allocatable ::  matrix2(:,:)
      Complex*16, allocatable ::  v_final(:,:)
      Complex*16, allocatable ::  mat_val(:)
      Complex*16, allocatable ::  vec_new(:)
      Complex*16, allocatable ::  vec2(:)
      Complex*16, allocatable ::  vec2_new(:)
      Complex*16, allocatable ::  resf(:)
      Complex*16, allocatable ::  hess(:,:)
      Complex*16, allocatable ::  hess2(:,:)
      Complex*16, allocatable ::  z(:,:)
      Complex*16, allocatable ::  final_evctrs(:,:)
      Complex*16, allocatable ::  vl(:,:)
      Complex*16, allocatable ::  vr(:,:)
      Complex*16, allocatable ::  evals(:)
      Complex*16, allocatable ::  evals_sel(:)
      Complex*16, allocatable ::  eigenw_temp(:)
      Complex*16, allocatable ::  work_hess(:)
      Complex*16, allocatable ::  work_zhseqr(:)
      Complex*16, allocatable ::  Pb_vec(:)

      Complex*16, allocatable ::  KSL_hilf(:)
      Complex*16, allocatable ::  KSL_hilf_p0(:)
      Complex*16, allocatable ::  KSL_conjg_hilf(:)
      Complex*16, allocatable ::  KSL_conjg_hilf_p0(:)

      Complex*16, allocatable ::  KSL_mat_vals_coo(:)
      Complex*16, allocatable ::  vals_store(:)
      Complex*16, allocatable ::  KSL_conjgmat_vals_coo(:)
      Complex*16, allocatable ::  KSL_mat_vals_csr_inter(:)
      Complex*16, allocatable ::  KSL_conjgmat_vals_csr_inter(:)
      integer*4, allocatable    ::  KSL_mat_rows_coo(:)
      integer*4, allocatable    ::  KSL_conjgmat_rows_coo(:)
      integer*4, allocatable    ::  KSL_mat_cols_coo(:)
      integer*4, allocatable    ::  KSL_conjgmat_cols_coo(:)
      integer*4, allocatable    ::  KSL_mat_cols_csr_inter(:)
      integer*4, allocatable    ::  KSL_conjgmat_cols_csr_inter(:)
      integer*4, allocatable    ::  KSL_mat_ind_csr_inter(:)
      integer*4, allocatable    ::  KSL_conjgmat_ind_csr_inter(:)
      integer*4, allocatable    ::  rows_store(:)
      integer*4, allocatable    ::  cols_store(:)
      integer*4, allocatable    ::  send_counts(:)
      integer*4, allocatable    ::  displ(:)


      Complex*16, allocatable ::  left_Loc_row(:)
      Complex*16, allocatable ::  Loo_row(:)
!      Complex*16, allocatable ::  KSL_hilf2_p0(:)
!      Complex*16, allocatable ::  KSL_conjg_hilf2_p0(:)
      Double precision, allocatable ::  rwork_hess(:)
      integer*4, allocatable    ::  mat_csr1(:)
      integer*4, allocatable    ::  mat_csr2(:)
      integer*4, allocatable    ::  proj_index(:)
      integer, allocatable    ::  ifaill(:)
      integer, allocatable    ::  ifailr(:)
      integer, allocatable    ::  pind(:)
      logical,allocatable     ::  sel_temp(:)
      Double precision 
     &                  rwork(maxncv), rd(maxncv,3)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character         bmat*1, which*2
      character(50)     file1,file2
      parameter         (file1='LR_MATRIX_Original.dat')
      integer           myid, numprocs, ido, n, nx, lworkl, info, j,
     &                  ierr, nconv, maxitr, ishfts, mode, i, k,
     &                  status_open, rows, sum, iter, iter2,infoo,
     &                  conv_count,lwork, ll, s ,q, ndigit, logfil,
     &                  mcaupd, counter, ii,kk,jj,ss,lll,qq,
     &                  last_row, start_from_upper, start_from_lower,
     &                  start_from_i,start_from_j, orb_pos,
     &                  start_counter, jmax, m, qqmax, kkmax, ssmax,
     &                  lllmax, mm, pos, pos2, st, ind, kkmax_conjg, 
     &                  run, last_jj, last_ll, last_kk, last_ii,
     &                  till, from, p  
      Integer*8         element_count 
      Integer*4         ind1, ind2, ind1_conjg, ind2_conjg, pp
      Integer*4         KSL_maxnonzero, rows_proj
      Complex*16
     &                  sigma,tempp
      Double precision
     &                  tol, a, b, ticks,  start, finish,time, 
     &                  start_gen, finish_gen, start_tot, finish_tot 
      logical           rvec, restart, exists, skip_upper_part
c
c     %----------------------------------------------%
c     | Local Buffers needed for MPI communication |
c     %----------------------------------------------%
c
      Complex*16
     &                  mv_buf(maxn)
c
c     %-----------------------------%
c     | BLAS & LAPACK routines used |
c     %-----------------------------%
c
      Double precision
     &                  pdznorm2
      external          pdznorm2, zaxpy 


!      Real*8, INTENT(IN)                                 :: time_LR   
      COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb), INTENT(IN) :: PSI
!      COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb) :: h2_PSI,h2_PSI_conjg
      COMPLEX*16, DIMENSION(Nconf), INTENT(IN) :: VIN
      COMPLEX*16, DIMENSION(2*Morb*NDX*NDY*NDZ) :: L_orb_row
      COMPLEX*16, DIMENSION(Morb,2*Morb*NDX*NDY*NDZ) :: L_orb_rows
      COMPLEX*16, DIMENSION(Morb,2*Morb*NDX*NDY*NDZ) :: L_orb_rows_2
      COMPLEX*16, DIMENSION(2*Morb*NDX*NDY*NDZ) :: COmat_row,COmat_row_2
      COMPLEX*16, DIMENSION(2*Morb*NDX*NDY*NDZ) :: COmat_bare_row
      COMPLEX*16, DIMENSION(Nconf) :: VOUT,Vhelp,Vhelp2
!      COMPLEX*16, DIMENSION(2*Nconf,nprocs) :: CI_row 
      COMPLEX*16, ALLOCATABLE :: CI_row(:,:) 
      COMPLEX*16, DIMENSION(2*Nconf) :: OCmat_row 
      COMPLEX*16, DIMENSION(Morb,2*Nconf) :: OCmat_rows 
      COMPLEX*16, DIMENSION(Morb,2*Nconf) :: OCmat_rows_2 
      COMPLEX*16, DIMENSION(2*Nconf) :: COmat_bare_col
      COMPLEX*16, DIMENSION(2*(Morb*NDX*NDY*NDZ+Nconf)) :: L_row
      COMPLEX*16, DIMENSION(NDX*NDY*NDZ) :: h2_out
      COMPLEX*16, allocatable :: KSL_temp(:,:)
      COMPLEX*16, allocatable :: KSL_temp_conjg(:,:)
!      COMPLEX*16, allocatable :: Tkin_2D(:,:)
      COMPLEX*16, DIMENSION(2*Morb,2*(Morb*NDX*NDY*NDZ)) 
     &                                    :: L_row_top
      COMPLEX*16, DIMENSION(2,2*(Morb*NDX*NDY*NDZ+Nconf)) 
     &                                    :: L_row_bottom
      COMPLEX*16, DIMENSION(2,2*(Morb*NDX*NDY*NDZ+Nconf)) 
     &                                    :: L_row_bottom_2
      LOGICAL                        :: bare   
      INTEGER           ::  NC, check, status, s_status(MPI_STATUS_SIZE)
      COMPLEX*16, DIMENSION(NDX*NDY*NDZ) :: Tkin_col
!      DOUBLE COMPLEX, allocatable, DIMENSION(:,:) :: LR_mat
      DOUBLE COMPLEX, allocatable, DIMENSION(:) :: something
      REAL*8 :: tol1_KSL, tol1_proj, t1, t2 , szgb
      COMPLEX*16, DIMENSION(2*NDX*NDY*NDZ*Morb) :: L_orb_col

      dimL_orb  =2*NDX*NDY*NDZ*Morb
      dimL = dimL_orb + 2*Nconf
      ND  = NDX*NDY*NDZ
      NC  = Nconf

      ldv=maxn
      allocate(v(ldv,maxncv), stat=ierr)    
          if(ierr /= 0) then
        write(*,*)"allocation error for (v(ldv,maxncv)", ldv, maxncv 
            szgb= ldv*maxncv*16*2.0/1024/1024/1024
       write(*,*)"Size v(ldv,maxncv) Gigabytes:", szgb
          stop "Err alloc v"
          else
        if(ierr == 0)write(*,*)"allocation ok for (v(ldv,maxncv) "
          endif

c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c 

      call CPU_TIME(start_tot)

!      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )

 
c
      ndigit = -3
      logfil = 6
      mcaupd = 1
c
c     %--------------------------------------------------%
c     | The number NX is the number of interior points   |
c     | in the discretization of the 2-dimensional       |
c     | convection-diffusion operator on the unit        |
c     | square with zero Dirichlet boundary condition.   | 
c     | The number N(=NX*NX) is the dimension of the     |
c     | matrix.  A standard eigenvalue problem is        |
c     | solved (BMAT = 'I').  NEV is the number of       |
c     | eigenvalues to be approximated.  The user can    |
c     | modify NX, NEV, NCV, WHICH to solve problems of  |
c     | different sizes, and to get different parts of   |
c     | the spectrum.  However, The following            |
c     | conditions must be satisfied:                    |
c     |                   N <= MAXN                      |
c     |                 NEV <= MAXNEV                    |
c     |           NEV + 2 <= NCV <= MAXNCV               | 
c     %--------------------------------------------------% 
c
      nx    = dimL 
      n     = nx!*nx 
!      nev   = 40
!      ncv   = 80
!      min_conv = 25
c
c     %--------------------------------------%
c     | Set up distribution of data to nodes |
c     %--------------------------------------%
c
      numprocs=nprocs
      nloc = (nx / nprocs)

c      if ( mod(nx, nprocs) .gt. myid ) nloc = nloc + 1


c
      if ( nloc .gt. maxn ) then
         print *, ' ERROR with _NDRV1: NLOC is greater than MAXN '
         go to 9000
      else if ( nev .gt. maxnev ) then
         print *, ' ERROR with _NDRV1: NEV is greater than MAXNEV '
         go to 9000
      else if ( ncv .gt. maxncv ) then
         print *, ' ERROR with _NDRV1: NCV is greater than MAXNCV '
         go to 9000
      end if
      bmat  = 'I'
      which = 'SM'
c
c     %---------------------------------------------------%
c     | The work array WORKL is used in ZNAUPD as         | 
c     | workspace.  Its dimension LWORKL is set as        |
c     | illustrated below.  The parameter TOL determines  |
c     | the stopping criterion. If TOL<=0, machine        |
c     | precision is used.  The variable IDO is used for  |
c     | reverse communication, and is initially set to 0. |
c     | Setting INFO=0 indicates that a random vector is  |
c     | generated to start the ARNOLDI iteration.         | 
c     %---------------------------------------------------%
c
      lworkl  = 3*ncv**2+5*ncv 
      tol    = -1.0e-10 
!      tol1   = 1.0e-08 
!      tol2   = 1.0e-08 
      ido    = 0

      restart=.FALSE.
      if(restart.eqv..FALSE.) then
        info   = 0  !random initial resid vector
      else
        info   = 0 
      end if 

c
c     %---------------------------------------------------%
c     | This program uses exact shift with respect to     |
c     | the current Hessenberg matrix (IPARAM(1) = 1).    |
c     | IPARAM(3) specifies the maximum number of Arnoldi |
c     | iterations allowed.  Mode 1 of ZNAUPD is used     |
c     | (IPARAM(7) = 1). All these options can be changed |
c     | by the user. For details see the documentation in |
c     | ZNAUPD.                                           |
c     %---------------------------------------------------%
c
      ishfts = 1
      maxitr = 300000
      mode   = 1
c
      iparam(1) = ishfts
      iparam(3) = maxitr 
      iparam(7) = mode 
      iparam(5) = 0 

      alpha=cmplx(1.0,0.0,kind=8)
      beta=cmplx(0.0,0.0,kind=8)


!==========================================================================================

      


      call CPU_TIME(start)

!AIS 2017
        call read_parse_V_W_Psi !parsing V_W_Psi_string.in file
         write(6,*) "File with user potential was parsed" 


       xlambda0=xlambda_0
!       CALL Get_InterPart(time)
       write(*,*) "xlambda_0",xlambda_0
       write(*,*) "xlambda0",xlambda0
      

       CALL VTRAP_EXT_TD(time)
       do i=1, NDX
         write(*,*) i,VTRAP_EXT(i)
       end do
         
!       CALL Get_h_W(PSI,time)
       do i=1, nprocs 
        if(myid.eq.i-1)  write(*,*) "myid",myid,InvZRIJ
            call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
       
       end do 
            call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
!       stop 


!======  Check if orbitals are Real
          orb_real=.TRUE.
          check=0 
          do j=1, Morb
            do i=1, ND
              if(Dimag(PSI(i,j)).gt.1e-12) check=check+1
            end do 
          end do 
          if(check.gt.0d0) orb_real=.FALSE.
          write(*,*) "ORBITALS REAL?",orb_real      




!======   SPECIAL PRECONSTRUCTION CALLS that are time-consuming
!=============  SAVE TKIN and KSL/WSL IF POSSIBLE TO SAVE TIME
        if(restart.eqv..TRUE.) then
          if(dimL.gt.maxdim_prec) preconstr=.TRUE.
          if(dimL.le.maxdim_prec) preconstr=.FALSE.
        else
          preconstr=.TRUE.
        end if

        if(preconstr.eqv..FALSE.) goto 434


         Tkin_prec=.TRUE.
         h2_prec=.TRUE.
         WSL_prec=.TRUE.
         KSL_prec=.TRUE.
         V_CIJ_prec=.FALSE.
         left_Loo_prec=.TRUE.
         left_Loc_prec=.TRUE.
         proj_prec=.TRUE.
 
         if(xlambda_0.eq.0d0) then
            WSL_prec=.FALSE.
            KSL_prec=.FALSE.
            V_CIJ_prec=.FALSE.
         end if 
  

          write(*,*) "Starting with preconstruction part..."

          if (Tkin_prec) then


!             INQUIRE(FILE="Tkin_1D.dat", EXIST=exists) 


!         1D kinetic energy 
           if(dim_mctdhb.eq.1) then

!             if(exists.eqv..TRUE.) then
!               open(unit=256, file='Tkin_1D.dat', 
!     .         status='old', action='readwrite')
!             else
!               open(unit=256, file='Tkin_1D.dat', 
!     .         status='replace', action='readwrite')
!             end if

              allocate(Tkin_1D(NDX*NDY*NDZ,NDX*NDY*NDZ))  !on all PEs
              call CPU_TIME(start)
                 if(Time_DVRMETHODX.eq.4) then !FFT CASE
                   counter=1
                   do while (counter.le.NDX-nprocs+1)
                         Tkin_col=0d0
                         call T_FFT_GENERAL(Tkin_col,counter+myid)
                         do i=1, nprocs
                           if(myid.eq.i-1) 
     .                        Tkin_1D(:,counter+myid)=Tkin_col(:)
                           call MPI_BCAST(Tkin_1D(:,counter+i-1),
     .                             ND,MPI_DOUBLE_COMPLEX,
     .                             i-1,MPI_COMM_WORLD,ierr)
                         end do
                         write(*,*) "FFT col",counter+myid
!                            do K=1, NDX
!                             if(abs(Tkin_1D(k,counter))
!     .                                          .gt.tol1) then 
!                              write(256,'(2I10,10x,2F20.10)') k,
!     .                         counter,
!     .                         Real(Tkin_1D(k,counter))
!     .                        ,Dimag(Tkin_1D(k,counter))
!                             end if
!                            end do
                         counter=counter+nprocs
                   end do
                   do while (counter.le.NDX.and.counter.gt.NDX-nprocs+1)
                         Tkin_col=0d0
                         call T_FFT_GENERAL(Tkin_col,counter)
                         Tkin_1D(:,counter)=Tkin_col(:)
                         if(myid.eq.0) write(*,*) "FFT col",counter
!                            do K=1, NDX
!                             if(abs(Tkin_1D(k,counter))
!     .                                          .gt.tol1) then 
!                              write(256,'(2I10,10x,2F20.10)') k,
!     .                         counter,
!     .                         Real(Tkin_1D(k,counter))
!     .                        ,Dimag(Tkin_1D(k,counter))
!                             end if
!                            end do
                         counter=counter+1
                   end do
                 else
                   Tkin_1D=Op_x
                 end if
            call CPU_TIME(finish)
            write(*,*) "Tkin 1D done",finish-start
!            close(256)
           end if



!         2D kinetic energy
           if(dim_mctdhb.eq.2) then
                if(.not.allocated(Tkin_2D).and.myid.eq.0) then
                   allocate(Tkin_2D(NDX*NDY*NDZ,NDX*NDY*NDZ))
                   Tkin_2D=0d0
                end if     
        
            if(.not.allocated(Tkin_col_p0))
     .       allocate(Tkin_col_p0(ND)) 

               start_counter=1

!            if(myid.eq.0) then 

!             INQUIRE(FILE="Tkin_2D.dat", EXIST=exists) 

!             if(exists.eqv..TRUE.) then
!               open(unit=256, file='Tkin_2D.dat', 
!     .         status='old', action='readwrite')

!               jmax=1
!               do 
!                 read(256,*,iostat=st) i,j,a,b
!                 if(st.ne.0) exit
!                 Tkin_2D(i,j)=Cmplx(a,b,kind=8)
!                 if(j.gt.jmax) jmax=j 
!               end do  

!               start_counter=jmax

!             else
!               open(unit=256, file='Tkin_2D.dat', 
!     .         status='replace', action='readwrite')

!               start_counter=1

!             end if  
!            end if

            CALL MPI_BCAST(start_counter,1,MPI_INTEGER,0,
     .             MPI_COMM_WORLD,ierr)
            call MPI_BARRIER(MPI_COMM_WORLD,ierr) 

!            if(start_counter.eq.NDX*NDY) goto 445

             call CPU_TIME(start)
                 if(Time_DVRMETHODX.eq.4.and.Time_DVRMETHODY.eq.4) then !FFT CASE
                    counter=start_counter
                       do while (counter.le.NDX*NDY-nprocs+1)
                         Tkin_col=0d0
                         call T_FFT_GENERAL(Tkin_col,counter+myid)

!                         For p=0 
                          if(myid.eq.0) then
                           Tkin_2D(:,counter)=Tkin_col(:)
!                            do K=1, NDX*NDY
!                             if(abs(Tkin_2D(k,counter))
!     .                                          .gt.tol1) then 
!                              write(256,'(2I10,10x,2F20.10)') k,
!     .                         counter,
!     .                         Real(Tkin_2D(k,counter))
!     .                        ,Dimag(Tkin_2D(k,counter))
!                             end if
!                            end do
                          end if

!                        For all other PEs 
                         do i=2, nprocs
                          Tkin_col_p0=0d0
                          if(myid.eq.(i-1)) Tkin_col_p0=Tkin_col
                          if(myid.eq.0) CALL MPI_RECV(Tkin_col_p0(:),
     .                      ND,MPI_DOUBLE_COMPLEX,i-1,0,MPI_COMM_WORLD,
     .                      s_status,ierr)
                          if(myid.eq.i-1) CALL MPI_SEND(Tkin_col_p0(:),
     .                    ND,MPI_DOUBLE_COMPLEX,0,0,MPI_COMM_WORLD,ierr)

                          if(myid.eq.0) then
                           Tkin_2D(:,counter+i-1)=Tkin_col_p0(:)
!                            do K=1, NDX*NDY
!                             if(abs(Tkin_2D(k,counter+i-1))
!     .                                          .gt.tol1) then 
!                              write(256,'(2I10,10x,2F20.10)') k,
!     .                         counter+i-1,
!     .                         Real(Tkin_2D(k,counter+i-1))
!     .                        ,Dimag(Tkin_2D(k,counter+i-1))
!                             end if
!                            end do
                          end if
                          call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
                         end do 
                         write(*,*) "FFT col",counter+myid
                         counter=counter+nprocs
                       end do

                       if(counter.gt.NDX*NDY-nprocs+1.and.
     .                        counter.le.NDX*NDY.and.myid.eq.0) then
                           do while(counter.le.NDX*NDY)
                             Tkin_col=0d0
                             call T_FFT_GENERAL(Tkin_col,counter)
                             Tkin_2D(:,counter)=Tkin_col(:)
!                             if(myid.eq.0) then
!                               do K=1, NDX*NDY
!                                if(abs(Tkin_2D(k,counter)).gt.tol1) then
!                                 write(256,'(2I10,10x,2F20.10)') k,
!     .                             counter
!     .                            ,Real(Tkin_2D(k,counter))
!     .                            ,Dimag(Tkin_2D(k,counter))
!                                end if
!                               end do
!                             end if
                             write(*,*) "FFT col",counter
                             counter=counter+1 
                           end do 
                       end if 
                       call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
                 else  !CASE OF HO-DVR
                   call get_Tkin_2D(Tkin_2D)
                 end if
445         continue
            close(256)
            call CPU_TIME(finish)
            write(*,*) "Tkin 2D done",finish-start
           end if
           call MPI_BARRIER(MPI_COMM_WORLD,ierr) 

          end if


!      Projector prec
          IF(proj_prec.eqv..TRUE.) then
       call CPU_TIME(t1)
           
            tol1_proj=tol1 

              if(.not.allocated(proj_index)) 
     .             allocate(proj_index(dimL_orb/2+1))

            IF(myid.eq.0) then

              if(.not.allocated(proj_vals_coo)) 
     .             allocate(proj_vals_coo(Morb*ND*ND))
              if(.not.allocated(proj_cols_coo)) 
     .             allocate(proj_cols_coo(Morb*ND*ND))
              if(.not.allocated(proj_rows_coo)) 
     .             allocate(proj_rows_coo(Morb*ND*ND))

!             INQUIRE(FILE="proj.dat", EXIST=exists) 

!             if(exists.eqv..TRUE.) then
!               open(unit=256, file='proj.dat', 
!     .         status='old', action='readwrite')

!               k=1
!               last_ii=1  
!               proj_index=1
!               do 
!                 read(256,*,iostat=st) i,j,a,b
!                 if(st.ne.0) exit
!                 proj_vals_coo(k)=Cmplx(a,b,kind=8)
!                 proj_cols_coo(k)=j
!                 proj_rows_coo(k)=i
!                 if(i.ne.last_ii) then
!                   proj_index(i)=k  
!                   last_ii=i
!                 end if 
!                 k=k+1
!               end do  
!               proj_nonzero=k-1
!               proj_index(dimL_orb/2+1)=proj_nonzero+1

!             else
!               open(unit=256, file='proj.dat', 
!     .         status='replace', action='readwrite')
!             end if 

            END IF

            call MPI_BCAST(exists,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
 

            if(.not.allocated(proj_vals_inter)) 
     .             allocate(proj_vals_inter(ND))
            if(.not.allocated(proj_cols_inter)) 
     .             allocate(proj_cols_inter(ND))
            if(.not.allocated(proj_rows_inter)) 
     .             allocate(proj_rows_inter(ND))





            if(.not.allocated(pind)) allocate(pind(nprocs)) 

!            if(exists.eqv..TRUE.) goto 4000 

            if(.not.allocated(Pb_vec)) allocate(Pb_vec(dimL_orb)) 

            i=1
            pind=0
            last_ii=1
            proj_index=1 
            counter=1
            proj_nonzero=1
            do while (i.le.dimL_orb/2-nprocs+1)
              Pb_vec=0d0  
              proj_vals_inter=0d0
              proj_rows_inter=0
              proj_cols_inter=0
              pos=mod(i+myid,ND)


              if(pos.eq.0) pos=ND
              do k=1, 2*Morb
                 if(i+myid.le.k*ND) then
                  orb_pos=k
                  exit
                 end if
              end do
              if(orb_pos.gt.Morb) orb_pos=orb_pos-Morb
              call get_proj_col(PSI,Pb_vec,i+myid,pos)
              Pb_vec=Conjg(Pb_vec) 

              do k=1, nprocs
                pind(k)=0
              end do 

              if(myid.eq.0) proj_index(i)=counter

              do k=1+(orb_pos-1)*ND, orb_pos*ND
                if(myid.eq.0) then
                  if(abs(Pb_vec(k)).gt.tol1_proj) then
                   proj_vals_coo(counter)=Pb_vec(k)
                   proj_rows_coo(counter)=i+myid
                   proj_cols_coo(counter)=k
!                   write(256,'(2I10,10x,2F20.10)') 
!     &                              proj_rows_coo(counter), 
!     &                              proj_cols_coo(counter), 
!     &                         Real(proj_vals_coo(counter)), 
!     &                         Dimag(proj_vals_coo(counter))
                   pind(1)=pind(1)+1
                   counter=counter+1
                  end if
                else
                  if(abs(Pb_vec(k)).gt.tol1_proj) then
                   pind(myid+1)=pind(myid+1)+1
                   proj_vals_inter(pind(myid+1))=Pb_vec(k)
                   proj_rows_inter(pind(myid+1))=i+myid
                   proj_cols_inter(pind(myid+1))=k
                  end if
                end if
              end do 

             
            do k=2, nprocs
              if(myid.eq.0) CALL MPI_RECV(pind(k),1,MPI_INTEGER,k-1,2,
     .                       MPI_COMM_WORLD,s_status,ierr) 
              if(myid.eq.k-1) CALL MPI_SEND(pind(k),1,MPI_INTEGER,0,2,
     .                       MPI_COMM_WORLD,ierr)
  

             if(pind(k).ne.0) then
               if(myid.eq.0) CALL MPI_RECV(proj_vals_inter(1:pind(k)),
     .                          pind(k),MPI_DOUBLE_COMPLEX,k-1,0,
     .                          MPI_COMM_WORLD,s_status,ierr) 
               if(myid.eq.k-1) CALL MPI_SEND(proj_vals_inter(1:pind(k)
     .                         ),pind(k),MPI_DOUBLE_COMPLEX,0,0,
     .                         MPI_COMM_WORLD,ierr)
               if(myid.eq.0) CALL MPI_RECV(proj_cols_inter(1:pind(k)),
     .                       pind(k),MPI_INTEGER,k-1,1,
     .                       MPI_COMM_WORLD,s_status,ierr) 
               if(myid.eq.k-1) CALL MPI_SEND(proj_cols_inter(1:pind(k)
     .                         ),pind(k),MPI_INTEGER,0,1,
     .                         MPI_COMM_WORLD,ierr)
               if(myid.eq.0) CALL MPI_RECV(proj_rows_inter(1:pind(k)),
     .                       pind(k),MPI_INTEGER,k-1,3,
     .                       MPI_COMM_WORLD,s_status,ierr) 
               if(myid.eq.k-1) CALL MPI_SEND(proj_rows_inter(1:pind(k)
     .                         ),pind(k),MPI_INTEGER,0,3,
     .                         MPI_COMM_WORLD,ierr)

                if(myid.eq.0) then
                  proj_index(i+k-1)=counter 
                  do ii=1, pind(k)
                   proj_vals_coo(counter)=proj_vals_inter(ii)
                   proj_cols_coo(counter)=proj_cols_inter(ii)
                   proj_rows_coo(counter)=proj_rows_inter(ii)
!                 write(256,'(2I10,10x,2F20.10)') proj_rows_coo(counter),
!     &                              proj_cols_coo(counter), 
!     &                              Real(proj_vals_coo(counter)), 
!     &                              Dimag(proj_vals_coo(counter))
                   counter=counter+1
                  end do
                  proj_vals_inter=0d0
                  proj_cols_inter=0
                  proj_rows_inter=0 
                end if


             end if 
            end do

              
              i=i+nprocs
            end do


           if(i.gt.dimL_orb/2-nprocs+1.and.i.le.dimL_orb/2.and.
     .             myid.eq.0) then

            do while(i.le.dimL_orb/2) 

              proj_index(i)=counter

              Pb_vec=0d0  
              pos=mod(i,ND)
              if(pos.eq.0) pos=ND
              do k=1, 2*Morb
                 if(i.le.k*ND) then
                  orb_pos=k
                  exit
                 end if
              end do
              if(orb_pos.gt.Morb) orb_pos=orb_pos-Morb

              call get_proj_col(PSI,Pb_vec,i,pos)
              Pb_vec=Conjg(Pb_vec) 

              do k=1+(orb_pos-1)*ND, orb_pos*ND
                  if(abs(Pb_vec(k)).gt.tol1_proj) then
                   pind(1)=pind(1)+1
                   proj_vals_coo(counter)=Pb_vec(k)
                   proj_rows_coo(counter)=i
                   proj_cols_coo(counter)=k
!                   write(256,'(2I10,10x,2F20.10)') 
!     &                              proj_rows_coo(counter), 
!     &                              proj_cols_coo(counter), 
!     &                         Real(proj_vals_coo(counter)), 
!     &                         Dimag(proj_vals_coo(counter))
                   counter=counter+1
                  end if
              end do
               

              i=i+1
            end do 

           end if


           if(myid.eq.0) then
              proj_nonzero=counter-1 
              proj_index(dimL_orb/2+1)=proj_nonzero+1 
           end if


4000        continue


              if(.not.allocated(send_counts))
     .          allocate(send_counts(nprocs))
              if(.not.allocated(displ))
     .          allocate(displ(nprocs))
              if(.not.allocated(proj_nonzeros))
     .          allocate(proj_nonzeros(nprocs))

              call MPI_BCAST(proj_index,dimL_orb/2+1,MPI_INTEGER,0,
     .                           MPI_COMM_WORLD,ierr) 
              call MPI_BCAST(proj_nonzero,1,MPI_INTEGER,0,
     .                           MPI_COMM_WORLD,ierr) 


              rows_proj=(dimL_orb/2)/nprocs
              if(myid.eq.0) rows_proj=rows_proj+mod((dimL_orb/2),nprocs)

              if(myid.eq.0) then  
                if(.not.allocated(projdist_vals_coo)) 
     .            allocate(projdist_vals_coo(proj_index(rows_proj+1)-1))
                if(.not.allocated(projdist_cols_coo)) 
     .            allocate(projdist_cols_coo(proj_index(rows_proj+1)-1))
                if(.not.allocated(projdist_rows_coo)) 
     .            allocate(projdist_rows_coo(proj_index(rows_proj+1)-1))

                send_counts(1)=
     .                proj_index(rows_proj+1)-1
                displ(1)=0 


              end if

              rows_proj=(dimL_orb/2)/nprocs

               do k=2, nprocs
                if(myid.eq.k-1) then
               from=proj_index(mod(dimL_orb/2,nprocs)+(k-1)*rows_proj+1)
               till=proj_index(mod(dimL_orb/2,nprocs)+k*rows_proj+1)-1 
                if(.not.allocated(projdist_vals_coo)) 
     .            allocate(projdist_vals_coo(till-from+1))
                if(.not.allocated(projdist_cols_coo)) 
     .            allocate(projdist_cols_coo(till-from+1))
                if(.not.allocated(projdist_rows_coo)) 
     .            allocate(projdist_rows_coo(till-from+1))

                 send_counts(k)=till-from+1
                 displ(k)=from-1
                end if
               end do

            do k=1, nprocs
              call MPI_BCAST(displ(k),1,MPI_INTEGER,k-1,
     .                           MPI_COMM_WORLD,ierr) 
              call MPI_BCAST(send_counts(k),1,MPI_INTEGER,k-1,
     .                           MPI_COMM_WORLD,ierr) 
            end do
            proj_nonzeros=send_counts  


          call MPI_SCATTERV(proj_vals_coo,send_counts,displ,
     .               MPI_DOUBLE_COMPLEX,projdist_vals_coo,
     .               send_counts(myid+1),
     .               MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)

          call MPI_SCATTERV(proj_rows_coo,send_counts,displ,
     .               MPI_INTEGER,projdist_rows_coo,
     .               send_counts(myid+1),
     .               MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

          call MPI_SCATTERV(proj_cols_coo,send_counts,displ,
     .               MPI_INTEGER,projdist_cols_coo,
     .               send_counts(myid+1),
     .               MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

           do k=2, nprocs
            if(myid.eq.k-1) then
             projdist_rows_coo=projdist_rows_coo-
     .          (mod(dimL_orb/2,nprocs)+(k-1)*rows_proj)
            end if
           end do 

           call MPI_BARRIER(MPI_COMM_WORLD,ierr)




           if(allocated(vals_store)) deallocate(vals_store)
           if(allocated(rows_store)) deallocate(rows_store)
           if(allocated(cols_store)) deallocate(cols_store)
           if(allocated(proj_vals_inter)) deallocate(proj_vals_inter)
           if(allocated(proj_rows_inter)) deallocate(proj_rows_inter)
           if(allocated(proj_cols_inter)) deallocate(proj_cols_inter)

            if(myid.eq.0) close(256)
            write(*,*) "Projector prec done!" 
       call CPU_TIME(t2)
          END IF         



!         h2_PSI and conjugate
          if(h2_prec) then 


          if(.not.allocated(h2_PSI)) allocate(h2_PSI(NDX*NDY*NDZ,Morb))
          if(.not.allocated(h2_PSI_conjg)) 
     .         allocate(h2_PSI_conjg(NDX*NDY*NDZ,Morb))
          call CPU_TIME(start)
          CALL get_h2_PSI(PSI,h2_PSI)
          CALL get_h2_PSI(Conjg(PSI),h2_PSI_conjg)
          call CPU_TIME(finish)
          write(*,*) "get_h2_PSI done",finish-start
          end if


!         WSL matrix

          if (WSL_prec) then

          if(xlambda_0.ne.0d0) then
                if(.not.allocated(WSL_mat)) 
     .            allocate(WSL_mat(NDX*NDY*NDZ,Morb,Morb))
          call CPU_TIME(start)
            do i=1, Morb
              do j=1, Morb
                CALL  Get_WSL_omp(WSL_mat(:,i,j),psi(:,i),psi(:,j)) 
              end do
            end do 
          end if 
          call CPU_TIME(finish)
          write(*,*) "get_WSL_omp done",finish-start

          end if    


!         KSL matrices array

          if(KSL_prec) then       

          IF(xlambda_0.ne.0d0) then

               if(DIM_MCTDHB.eq.1) KSL_maxnonzero=7e5
               if(DIM_MCTDHB.eq.2) then
                 if(NDX*NDY.le.64*64) KSL_maxnonzero=1.0e7 
                 if(NDX*NDY.gt.64*64.and.NDX*NDY.le.256*256)
     .                   KSL_maxnonzero=6.25e7 
               end if
 
            tol1_KSL=tol1 

            if(myid.eq.0) then
                if(.not.allocated(KSL_mat_vals_csc).and.myid.eq.0) 
     .            allocate(KSL_mat_vals_csc(KSL_maxnonzero,Morb,Morb))
                if(.not.allocated(KSL_mat_ind_csc).and.myid.eq.0) 
     .            allocate(KSL_mat_ind_csc(NDX*NDY*NDZ+1,Morb,Morb))
                if(.not.allocated(KSL_mat_rows_csc).and.myid.eq.0) 
     .            allocate(KSL_mat_rows_csc(KSL_maxnonzero,Morb,Morb))

                if(.not.allocated(KSL_conjgmat_vals_csc).and.myid.eq.0) 
     .         allocate(KSL_conjgmat_vals_csc(KSL_maxnonzero,Morb,Morb))
                if(.not.allocated(KSL_conjgmat_ind_csc).and.myid.eq.0) 
     .         allocate(KSL_conjgmat_ind_csc(NDX*NDY*NDZ+1,Morb,Morb))
                if(.not.allocated(KSL_conjgmat_rows_csc).and.myid.eq.0) 
     .         allocate(KSL_conjgmat_rows_csc(KSL_maxnonzero,Morb,Morb))


               KSL_mat_vals_csc=0d0
               KSL_mat_ind_csc=1
               KSL_mat_rows_csc=0
               KSL_conjgmat_vals_csc=0d0
               KSL_conjgmat_ind_csc=1
               KSL_conjgmat_rows_csc=0
            end if


            if(.not.allocated(KSL_hilf)) allocate(KSL_hilf(ND))
            if(.not.allocated(KSL_hilf_p0)) allocate(KSL_hilf_p0(ND))
            if(.not.allocated(KSL_conjg_hilf)) 
     .                     allocate(KSL_conjg_hilf(ND))
            if(.not.allocated(KSL_conjg_hilf_p0)) 
     .                     allocate(KSL_conjg_hilf_p0(ND))

             ll=1; kk=1; ii=1; jj=1; kkmax=1; kkmax_conjg=1;
             last_ii=1; last_kk=1; last_ll=1; last_jj=1;
             ind1=1
             ind1_conjg=1

!          if(myid.eq.0) then
!           INQUIRE(FILE="KSL_mat.dat", EXIST=exists) 

!           if(exists.eqv..FALSE.) then   !File does not exist so far

!            open(unit=256, file='KSL_mat.dat', 
!     .        status='replace', action='readwrite')
!            open(unit=257, file='KSL_conjg_mat.dat', 
!     .        status='replace', action='readwrite')


!           else 

!             open(unit=256, file='KSL_mat.dat', 
!     .        status='old', action='readwrite')
!             open(unit=257, file='KSL_conjg_mat.dat', 
!     .        status='old', action='readwrite')
           
!                if(.not.allocated(KSL_mat_vals_coo)) 
!     .            allocate(KSL_mat_vals_coo(KSL_maxnonzero))
!                if(.not.allocated(KSL_mat_cols_coo)) 
!     .            allocate(KSL_mat_cols_coo(KSL_maxnonzero))
!                if(.not.allocated(KSL_mat_rows_coo)) 
!     .            allocate(KSL_mat_rows_coo(KSL_maxnonzero))

!                if(.not.allocated(KSL_mat_vals_csr_inter))
!     .            allocate(KSL_mat_vals_csr_inter(KSL_maxnonzero))
!                if(.not.allocated(KSL_mat_cols_csr_inter))
!     .            allocate(KSL_mat_cols_csr_inter(KSL_maxnonzero))
!                if(.not.allocated(KSL_mat_ind_csr_inter)) 
!     .            allocate(KSL_mat_ind_csr_inter(ND+1))

!                if(.not.allocated(KSL_conjgmat_vals_coo)) 
!     .            allocate(KSL_conjgmat_vals_coo(KSL_maxnonzero))
!                if(.not.allocated(KSL_conjgmat_cols_coo)) 
!     .            allocate(KSL_conjgmat_cols_coo(KSL_maxnonzero))
!                if(.not.allocated(KSL_conjgmat_rows_coo)) 
!     .            allocate(KSL_conjgmat_rows_coo(KSL_maxnonzero))

!                if(.not.allocated(KSL_conjgmat_vals_csr_inter)) 
!     .            allocate(KSL_conjgmat_vals_csr_inter(KSL_maxnonzero))
!                if(.not.allocated(KSL_conjgmat_cols_csr_inter)) 
!     .            allocate(KSL_conjgmat_cols_csr_inter(KSL_maxnonzero))
!                if(.not.allocated(KSL_conjgmat_ind_csr_inter)) 
!     .            allocate(KSL_conjgmat_ind_csr_inter(ND+1))

!             ind1=1
!             do 
!               read(256,*,iostat=st) ll,kk,ii,jj,a,b
!               if(st.ne.0) exit

!                 if(jj.ne.last_jj) then

!                   mjob(1)=2
!                   mjob(2)=1
!                   mjob(3)=1
!                   mjob(4)=2
!                   mjob(5)=KSL_maxnonzero
!                   mjob(6)=0
!           call mkl_zcsrcoo(mjob, ND, KSL_mat_vals_csr_inter, 
!     .             KSL_mat_cols_csr_inter,
!     .             KSL_mat_ind_csr_inter, ind1-1, KSL_mat_vals_coo, 
!     .             KSL_mat_rows_coo, KSL_mat_cols_coo, info) 


!                   mjob(1)=0
!                   mjob(2)=1
!                   mjob(3)=1
!                   mjob(4)=2
!                   mjob(5)=KSL_maxnonzero
!                   mjob(6)=1
!           call mkl_zcsrcsc(mjob, ND, KSL_mat_vals_csr_inter, 
!     .             KSL_mat_cols_csr_inter, KSL_mat_ind_csr_inter,
!     .             KSL_mat_vals_csc(:,last_ii,last_jj), 
!     .             KSL_mat_rows_csc(:,last_ii,last_jj),
!     .             KSL_mat_ind_csc(:,last_ii,last_jj), info)


!                   kkmax=1
!                   ind1=1
!                   KSL_mat_vals_coo=0d0
!                   KSL_mat_rows_coo=0        
!                   KSL_mat_cols_coo=0
!                   KSL_mat_vals_csr_inter=0d0
!                   KSL_mat_cols_csr_inter=0
!                   KSL_mat_ind_csr_inter=1         

!                   KSL_mat_vals_coo(ind1)=Cmplx(a,b,kind=8)
!                   KSL_mat_rows_coo(ind1)=ll
!                   KSL_mat_cols_coo(ind1)=kk
!                   ind1=ind1+1
!                   last_ll=ll 
!                   last_kk=kk 
!                   last_ii=ii 
!                   last_jj=jj 
!                   if(kk.gt.kkmax) kkmax=kk 
           
!                 else
!                    KSL_mat_vals_coo(ind1)=Cmplx(a,b,kind=8)
!                    KSL_mat_rows_coo(ind1)=ll
!                    KSL_mat_cols_coo(ind1)=kk
!                    ind1=ind1+1
!                    last_ll=ll 
!                    last_kk=kk 
!                    last_ii=ii 
!                    last_jj=jj 
!                    if(kk.gt.kkmax) kkmax=kk 
!                 end if

!             end do  !Take last k,i,j as starting point to continue

!             counter=ind1-1
!             mjob(1)=2
!             mjob(2)=1
!             mjob(3)=1
!             mjob(4)=2
!             mjob(5)=KSL_maxnonzero
!             mjob(6)=0
!           call mkl_zcsrcoo(mjob, ND, KSL_mat_vals_csr_inter, 
!     .             KSL_mat_cols_csr_inter,
!     .             KSL_mat_ind_csr_inter, counter, KSL_mat_vals_coo, 
!     .             KSL_mat_rows_coo, KSL_mat_cols_coo, info) 


!             mjob(1)=0
!             mjob(2)=1
!             mjob(3)=1
!             mjob(4)=2
!             mjob(5)=KSL_maxnonzero
!             mjob(6)=1
!           call mkl_zcsrcsc(mjob, ND, KSL_mat_vals_csr_inter, 
!     .             KSL_mat_cols_csr_inter, KSL_mat_ind_csr_inter,
!     .             KSL_mat_vals_csc(:,last_ii,last_jj), 
!     .             KSL_mat_rows_csc(:,last_ii,last_jj),
!     .             KSL_mat_ind_csc(:,last_ii,last_jj), info)



!             ind1_conjg=1
!             last_ii=1; last_kk=1; last_ll=1; last_jj=1;
!             do 
!               read(257,*,iostat=st) ll,kk,ii,jj,a,b
!               if(st.ne.0) exit

!                 if(jj.ne.last_jj) then

!                   mjob(1)=2
!                   mjob(2)=1
!                   mjob(3)=1
!                   mjob(4)=2
!                   mjob(5)=KSL_maxnonzero
!                   mjob(6)=0
!           call mkl_zcsrcoo(mjob, ND, KSL_conjgmat_vals_csr_inter, 
!     .             KSL_conjgmat_cols_csr_inter,
!     .             KSL_conjgmat_ind_csr_inter, ind1_conjg-1,
!     .             KSL_conjgmat_vals_coo, 
!     .             KSL_conjgmat_rows_coo, KSL_conjgmat_cols_coo, info) 


!                   mjob(1)=0
!                   mjob(2)=1
!                   mjob(3)=1
!                   mjob(4)=2
!                   mjob(5)=KSL_maxnonzero
!                   mjob(6)=1
!           call mkl_zcsrcsc(mjob, ND, KSL_conjgmat_vals_csr_inter, 
!     .             KSL_conjgmat_cols_csr_inter,
!     .             KSL_conjgmat_ind_csr_inter,
!     .             KSL_conjgmat_vals_csc(:,last_ii,last_jj), 
!     .             KSL_conjgmat_rows_csc(:,last_ii,last_jj),
!     .             KSL_conjgmat_ind_csc(:,last_ii,last_jj), info)

!                   ind1_conjg=1
!                   kkmax_conjg=1
!                   KSL_conjgmat_vals_coo=0d0
!                   KSL_conjgmat_rows_coo=0        
!                   KSL_conjgmat_cols_coo=0
!                   KSL_conjgmat_vals_csr_inter=0d0
!                   KSL_conjgmat_cols_csr_inter=0
!                   KSL_conjgmat_ind_csr_inter=1         

!                   KSL_conjgmat_vals_coo(ind1_conjg)=Cmplx(a,b,kind=8)
!                   KSL_conjgmat_rows_coo(ind1_conjg)=ll
!                   KSL_conjgmat_cols_coo(ind1_conjg)=kk
!                   ind1_conjg=ind1_conjg+1
!                   last_ll=ll 
!                   last_kk=kk 
!                   last_ii=ii 
!                   last_jj=jj 
!                   if(kk.gt.kkmax_conjg) kkmax_conjg=kk 
           
!                 else
!                    KSL_conjgmat_vals_coo(ind1_conjg)=Cmplx(a,b,kind=8)
!                    KSL_conjgmat_rows_coo(ind1_conjg)=ll
!                    KSL_conjgmat_cols_coo(ind1_conjg)=kk
!                    ind1_conjg=ind1_conjg+1
!                    last_ll=ll 
!                    last_kk=kk 
!                    last_ii=ii 
!                    last_jj=jj 
!                    if(kk.gt.kkmax_conjg) kkmax_conjg=kk 
!                 end if
!             end do  !Take last k,i,j as starting point to continue



!             counter=ind1_conjg-1
!             mjob(1)=2
!             mjob(2)=1
!             mjob(3)=1
!             mjob(4)=2
!             mjob(5)=KSL_maxnonzero
!             mjob(6)=0
!           call mkl_zcsrcoo(mjob, ND, KSL_conjgmat_vals_csr_inter, 
!     .             KSL_conjgmat_cols_csr_inter,
!     .             KSL_conjgmat_ind_csr_inter, counter, 
!     .             KSL_conjgmat_vals_coo, 
!     .             KSL_conjgmat_rows_coo, KSL_conjgmat_cols_coo, info) 


!             mjob(1)=0
!             mjob(2)=1
!             mjob(3)=1
!             mjob(4)=2
!             mjob(5)=KSL_maxnonzero
!             mjob(6)=1
!           call mkl_zcsrcsc(mjob, ND, KSL_conjgmat_vals_csr_inter, 
!     .             KSL_conjgmat_cols_csr_inter, 
!     .             KSL_conjgmat_ind_csr_inter,
!     .             KSL_conjgmat_vals_csc(:,last_ii,last_jj), 
!     .             KSL_conjgmat_rows_csc(:,last_ii,last_jj),
!     .             KSL_conjgmat_ind_csc(:,last_ii,last_jj), info)



!         if(allocated(KSL_mat_vals_coo))
!     .   deallocate(KSL_mat_vals_coo)   
!         if(allocated(KSL_conjgmat_vals_coo))
!     .   deallocate(KSL_conjgmat_vals_coo)   
!         if(allocated(KSL_mat_rows_coo))
!     .   deallocate(KSL_mat_rows_coo)   
!         if(allocated(KSL_conjgmat_rows_coo))
!     .   deallocate(KSL_conjgmat_rows_coo)   
!         if(allocated(KSL_mat_cols_coo))
!     .   deallocate(KSL_mat_cols_coo)   
!         if(allocated(KSL_conjgmat_cols_coo))
!     .   deallocate(KSL_conjgmat_cols_coo)   
!         if(allocated(KSL_mat_vals_csr_inter))
!     .   deallocate(KSL_mat_vals_csr_inter)   
!         if(allocated(KSL_conjgmat_vals_csr_inter))
!     .   deallocate(KSL_conjgmat_vals_csr_inter)   
!         if(allocated(KSL_mat_cols_csr_inter))
!     .   deallocate(KSL_mat_cols_csr_inter)   
!         if(allocated(KSL_conjgmat_cols_csr_inter))
!     .   deallocate(KSL_conjgmat_cols_csr_inter)   
!         if(allocated(KSL_mat_ind_csr_inter))
!     .   deallocate(KSL_mat_ind_csr_inter)   
!         if(allocated(KSL_conjgmat_ind_csr_inter))
!     .   deallocate(KSL_conjgmat_ind_csr_inter)   

!           end if
!          end if
         

          if(myid.eq.0) kk=min(kkmax,kkmax_conjg)+1

          CALL MPI_BCAST(exists,1,MPI_LOGICAL,0,
     .             MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(ll,1,MPI_INTEGER,0,
     .             MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(kk,1,MPI_INTEGER,0,
     .             MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(ii,1,MPI_INTEGER,0,
     .             MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(jj,1,MPI_INTEGER,0,
     .             MPI_COMM_WORLD,ierr)
          
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)


          call CPU_TIME(start)


          if(ii.eq.Morb.and.jj.eq.Morb.and.Morb.gt.1) goto 567

!            goto 490

              IF(myid.eq.0) then

                if(.not.allocated(KSL_temp)) 
     .            allocate(KSL_temp(NDX*NDY*NDZ,NDX*NDY*NDZ))
                if(.not.allocated(KSL_temp_conjg)) 
     .            allocate(KSL_temp_conjg(NDX*NDY*NDZ,NDX*NDY*NDZ))

               do i=1, Morb 
                 
                 KSL_temp=0d0
                 KSL_temp_conjg=0d0
                 CALL  Get_KSL(KSL_temp
     .                         ,psi(:,i))
                 if(orb_real) then
                   KSL_temp_conjg=KSL_temp
                 else
                   CALL  Get_KSL(KSL_temp_conjg,
     .                    CONJG(psi(:,i))) 
                 end if

                 do j=1, Morb


                   k=1
                   do while (k.le.NDX*NDY*NDZ)
                     KSL_mat_ind_csc(k,i,j)=ind1
                     KSL_conjgmat_ind_csc(k,i,j)=ind1_conjg
                     p=1 
                     do while (p.le.NDX*NDY*NDZ)
                  
                      if(abs(KSL_temp(k,p)*Psi(k,j)).gt.tol1_KSL) then
                       KSL_mat_vals_csc(ind1,i,j)=KSL_temp(k,p)*Psi(k,j)
                       KSL_mat_rows_csc(ind1,i,j)=p
!                        write(256,'(4I10,10x,2F20.10)') p,k,i,j,
!     .                     Real(KSL_temp(k,p)*Psi(k,j))
!     .                     ,Dimag(KSL_temp(k,p)*Psi(k,j))
                       ind1=ind1+1
                      end if 

                      if(abs(KSL_temp_conjg(k,p)*Psi(k,j)).gt.
     .                                             tol1_KSL) then
                       KSL_conjgmat_vals_csc(ind1_conjg,i,j)=
     .                                 KSL_temp_conjg(k,p)*Psi(k,j)
                       KSL_conjgmat_rows_csc(ind1_conjg,i,j)=p
!                        write(257,'(4I10,10x,2F20.10)') p,k,i,j,
!     .                     Real(KSL_temp_conjg(k,p)*Psi(k,j))
!     .                     ,Dimag(KSL_temp_conjg(k,p)*Psi(k,j))
                       ind1_conjg=ind1_conjg+1
                      end if 

                      p=p+1
                     end do

               


                     k=k+1
                   end do 
                  if(myid.eq.0) then 
                    KSL_mat_ind_csc(NDX*NDY*NDZ+1,i,j)=ind1
                    KSL_conjgmat_ind_csc(NDX*NDY*NDZ+1,i,j)=ind1_conjg
                  end if  
                  ind1=1
                  ind1_conjg=1
                 end do



               end do
                if(allocated(KSL_temp)) deallocate(KSL_temp)
                if(allocated(KSL_temp_conjg)) deallocate(KSL_temp_conjg)
!          if(myid.eq.0) close(256)
!          if(myid.eq.0) close(257)
              END IF


          call MPI_BARRIER(MPI_COMM_WORLD,ierr)


567       continue
          call CPU_TIME(finish)
          write(*,*) "get_KSL_mat done",finish-start
 


          END IF
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)


          if(allocated(KSL_hilf))
     .    deallocate(KSL_hilf)
          if(allocated(KSL_hilf_p0))
     .    deallocate(KSL_hilf_p0)
          if(allocated(KSL_conjg_hilf))
     .    deallocate(KSL_conjg_hilf)
          if(allocated(KSL_conjg_hilf_p0))
     .    deallocate(KSL_conjg_hilf_p0)


          end if 



!    END of KSL mat construction


!========  left_Loo construction
       if(left_Loo_prec) then 
 
     
          write(*,*) "STR 3.9"
         if(.not.allocated(KSL_dist))  then
          allocate(KSL_dist(ND,Morb,Morb),stat=ierr)
          szgb= ND*Morb*Morb*16*2.0/1024/1024/1024
          write(*,*)"Size KLS_dist Gigabytes:", szgb
          if(ierr /= 0) then
          write(*,*)"allocation error for KSL_dist"
          stop "Err alloc KLS_dist"
          else
          if(ierr == 0)write(*,*)"allocation ok for KSL_dist"
          endif
         endif 

           if(.not.allocated(KSL_conjg_dist))  then
              allocate(KSL_conjg_dist(ND,Morb,Morb), stat=ierr)
          szgb= ND*Morb*Morb*16*2.0/1024/1024/1024
          write(*,*)"Size KLS_conj_dist Gigabytes:", szgb
          if(ierr /= 0) then
          write(*,*)"allocation error for KSL_conj_dist"
          stop "Err alloc KLS_conj_dist"
          else
          if(ierr == 0)write(*,*)"allocation ok for KSL_conj_dist"
          endif
           endif 

           call CPU_TIME(start)
           call construct_left_of_Loo_upfront_sparse(PSI) !CALLED BY ALL PEs
           call CPU_TIME(finish)



           write(*,*) "construct left_Loo done",finish-start
           call MPI_BARRIER(MPI_COMM_WORLD,ierr)


       end if 


             
          write(*,*) "STR 4 ok"


!======== V_PrdCIJ Coefficient vectors needed for L_co/L_oc
!======== not in parallel because fast anyway, all PEs execute the same code 

!      if(Nconf*Morb*Morb*Morb*Morb.lt.7.0e8) then       !THIS IS BELOW 10 GB 

       if(V_CIJ_prec) then

         if(.not.allocated(V_PrdCIJ))  then
        allocate(V_PrdCIJ(Nconf,Morb,Morb,Morb,Morb), stat=ierr) 
          szgb= Nconf*Morb**4.0/1024/1024/1024
          write(*,*)"Size V_PrdCIJ Gigabytes:", szgb
          if(ierr /= 0) then
          write(*,*)"allocation error for V_PrdCIJ"
          stop "Err alloc V_PrdCIJ"
          else
          if(ierr == 0)write(*,*)"allocation ok for V_PrdCIJ"
          endif
           endif 
         if(.not.allocated(V_PrdCIJ2))  then
       allocate(V_PrdCIJ2(Nconf,Morb,Morb), stat=ierr) 
          szgb= Nconf*Morb**2/1024/1024/1024
          write(*,*)"Size V_PrdCIJ2 Gigabytes:", szgb
          if(ierr /= 0) then
          write(*,*)"allocation error for V_PrdCIJ2"
          stop "Err alloc V_PrdCIJ2"
          else
          if(ierr == 0)write(*,*)"allocation ok for V_PrdCIJ2"
          endif
           endif 
              write(6,*) "STR 4.1"

           call CPU_TIME(start)

         if(myid.eq.0) INQUIRE(FILE="V_PrdCIJ.dat", EXIST=exists) 
            CALL MPI_BCAST(exists,1,MPI_LOGICAL,0,
     .             MPI_COMM_WORLD,ierr)
!           INQUIRE(FILE="V_PrdCIJ2.dat", EXIST=exists) 

       if(myid.eq.0) then
        if(exists.eqv..FALSE.) then  !files do not exist       

           open(unit=256, file='V_PrdCIJ.dat', 
     .        status='replace', action='readwrite')
           open(unit=257, file='V_PrdCIJ2.dat', 
     .        status='replace', action='readwrite')


        else

           open(unit=256, file='V_PrdCIJ.dat', 
     .        status='old', action='readwrite')
           open(unit=257, file='V_PrdCIJ2.dat', 
     .        status='old', action='readwrite')


        end if 
       end if

           V_PrdCIJ=0d0  
           V_PrdCIJ2=0d0  

              write(6,*) "STR 4.2"
         do q=1, Morb
           do k=1, Morb
                 Vhelp=VIN;Vhelp2=0d0
!               CALL Produce_Cij(Vhelp,V_PrdCIJ2(:,q,k),q,k,0,0,Nc)
               CALL Produce_Cij(Vhelp,Vhelp2,q,k,0,0,Nc)
               V_PrdCIJ2(:,q,k)=Vhelp2
              if(myid.eq.0) then
               do s=1, NC
                  if(abs(V_PrdCIJ2(s,q,k)).gt.tol1) then 
                    write(257,'(3I10,10x,2F20.10)') q,k,s
     .                ,Real(V_PrdCIJ2(s,q,k))
     .                ,Dimag(V_PrdCIJ2(s,q,k))
                  end if
               end do           
!               write(*,*) V_PrdCIJ2(:,1,1)
!               stop
              end if      
              call MPI_BARRIER(MPI_COMM_WORLD,ierr)

             do s=1, Morb
               do ll=1, Morb 
                 Vhelp=VIN;Vhelp2=0d0
                CALL Produce_Cij(Vhelp,Vhelp2,q,k,s,ll,Nc)
                V_PrdCIJ(:,q,k,s,ll)=Vhelp2

              if(myid.eq.0) then
                  do ii=1, NC
                     if(abs(V_PrdCIJ(ii,q,k,s,ll)).gt.tol1) then 
                       write(256,'(5I10,10x,2F20.10)') q,k,s,ll,
     .                   ii,Real(V_PrdCIJ(ii,q,k,s,ll))
     .                  ,Dimag(V_PrdCIJ(ii,q,k,s,ll))
                     end if
                  end do                 
              end if
              call MPI_BARRIER(MPI_COMM_WORLD,ierr)
              end do 
             end do
           end do 
         end do

        if(myid.eq.0) close(256)
        if(myid.eq.0) close(257)


        call CPU_TIME(finish)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      write(*,*) "Produce_Cij done",finish-start
       end if


!======  left Loc construction
      if(left_Loc_prec.and.Morb.gt.1) then

!       Loc_maxnonzero=maxnonzero 


              write(6,*) "STR 4.3"
           call CPU_TIME(start)
           call construct_Loc_matrix_full_sparse(PSI,VIN,Nconf) 
           call CPU_TIME(finish)
              write(6,*) "STR 4.4"



        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        write(*,*) "construct left_Loc done",finish-start

      end if 

        write(*,*) "PRECONSTRUCTIONS DONE!!!"
434     continue


!================================================================================================
!===============  END OF PRECONSTRUCTION PART ===================================================
!================================================================================================





 
!===============================================================================



            
!================================================================================================
!===============  TESTING SECTION ===================================================

!      if(myid.eq.0) then
!           open(unit=256, file='orb_col.dat', 
!     .        status='replace', action='readwrite')

!        do i=1, dimL_orb/2
!             L_orb_col=0d0 
!             call construct_orb_col(PSI,L_orb_col,i)
        
!      call construct_full_row_orb(PSI,L_orb_rows,512)
!             do j=1, dimL_orb 
!              if(abs(L_orb_rows(1,j)).ge.tol1) then
!                write(256,'(I10,10x,2F20.10)') j,
!     .                   Real(L_orb_rows(1,j)),
!     .                   Dimag(L_orb_rows(1,j))
!              end if
!             end do
!        end do
!        close(256)
!       end if 
!           call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!           stop


!=============== END OF TESTING SECTION ===================================================
!================================================================================================




!=============== MATRIX CONSTRUCTION IF SMALL ENOUGH ===================================================
!=============  SWITCH: IF DIM <= 20000 --> Save the LR-matrix, only 1 construction of each row


!        ROWS of LR-matrix per PE
         if(myid.lt.nprocs-1) then 
           rows=nloc
         else
           rows=nloc+mod(nx,nprocs)
         end if    


!      if(dimL.le.maxdim_prec.and.restart.eqv..FALSE.) then 
 
          if(preconstr.eqv..FALSE.) then 
             if(.not.allocated(h2_PSI)) 
     .           allocate(h2_PSI(NDX*NDY*NDZ,Morb)) 
             if(.not.allocated(h2_PSI_conjg)) 
     .           allocate(h2_PSI_conjg(NDX*NDY*NDZ,Morb)) 
          end if

 
!        if(.not.allocated(LR_mat)) then
!          allocate(LR_mat(nx,nx),stat=ierr) 
!          if(ierr /= 0) then
!              write(*,*) myid,ierr,"allocation error in LR_MAT"
!              stop 
!          end if
!        end if  


!======= CONSTRUCTING UPPER PART OF LR-MATRIX  ===============================================


        call CPU_TIME(start_gen)

!      LR_mat=0d0
      skip_upper_part=.FALSE.
      start_from_upper=1
      start_from_lower=1
      element_count=0
      if(myid.eq.0) then 

       if(binary.eqv..FALSE.) 
     . INQUIRE(FILE="LR_matrix_save.dat", EXIST=exists) 
       if(binary.eqv..TRUE.) 
     . INQUIRE(FILE="LR_bin", EXIST=exists) 

      if(exists.eqv..TRUE.) then

      write(6,*)"CANNOT Continue STEP.1:  Matrix construction"
      write(6,*)"File(s) LR_matrix_save.dat or LR_bin already EXIST!"
      write(6,*)"To reconstruct the matrix remove them:"
          write(6,*)"rm -f LR_matrix_save.dat"
          write(6,*)"rm -f LR_bin"
      write(6,*)"For STEP.2 Arnolri diagonalization change:"
      write(6,*)"    task=2 in lr_arnoldi.in "
       stop("Needs your action...stop")
       if(binary.eqv..FALSE.) open(unit=256, file='LR_matrix_save.dat', 
     .   form='formatted', status='unknown')
     !!.        status='unknown', action='readwrite')
     !!.        status='old', action='readwrite')

       if(binary.eqv..TRUE.)Open (unit=260, File='LR_bin',
     .  Status='unknown', 
     .  Form='unformatted', Access='sequential')

          last_row=1
          do 
            if(binary.eqv..FALSE.) read(256,*,iostat=st) i,j,a,b
            if(binary.eqv..TRUE.) Read(260,iostat=st) i,j,tempp
            if(st.ne.0) exit 
            if(i.gt.last_row) last_row=i 
            
          end do  

!          last_row=i
          if(last_row.gt.dimL_orb) then
             skip_upper_part=.TRUE.
             start_from_lower=mod(i-dimL_orb+1,NC)
             if(start_from_lower.eq.0) start_from_lower=NC
          end if 
          if(last_row.lt.dimL_orb) then
             skip_upper_part=.FALSE.
             start_from_upper=mod(i+1,ND)
             if(start_from_upper.eq.0) start_from_upper=ND
          end if 
          if(last_row.eq.dimL_orb) then
             skip_upper_part=.TRUE.
             start_from_lower=1
          end if 

      else

       if(binary.eqv..FALSE.) open(unit=256, file='LR_matrix_save.dat', 
     .        status='replace', action='readwrite')

       if(binary.eqv..TRUE.)Open (unit=260, File='LR_bin',
     .  Status='unknown', 
     .  Form='unformatted', Access='sequential')


      end if 
      end if

         CALL MPI_BCAST(exists,1,MPI_LOGICAL,0,
     .             MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(last_row,1,MPI_INTEGER,0,
     .             MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(skip_upper_part,1,MPI_LOGICAL,0,
     .             MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(start_from_upper,1,MPI_INTEGER,0,
     .             MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(start_from_lower,1,MPI_INTEGER,0,
     .             MPI_COMM_WORLD,ierr)



!       goto 654
     
      if(skip_upper_part.eqv..FALSE.) then
      i=start_from_upper
      if(.not.allocated(left_Loo_hilf2)) 
     .     allocate(left_Loo_hilf2(Morb,dimL_orb))
      if(.not.allocated(left_Loo_hilf2_p0)) 
     .     allocate(left_Loo_hilf2_p0(Morb,dimL_orb))
      if(.not.allocated(left_Loc_hilf2)) 
     .     allocate(left_Loc_hilf2(Morb,2*Nconf))
      if(.not.allocated(left_Loc_hilf2_p0)) 
     .     allocate(left_Loc_hilf2_p0(Morb,2*Nconf))

      if(.not.allocated(Loo_u)) 
     .     allocate(Loo_u(dimL_orb/2,Morb*nprocs))
      if(.not.allocated(Loo_v)) 
     .     allocate(Loo_v(dimL_orb/2,Morb*nprocs))

            if(myid.eq.0) then
             if(.not.allocated(c1_mat)) 
     &        allocate(c1_mat(rows_proj+mod(dimL_orb/2,nprocs),
     &                       Morb*nprocs))
            else 
             if(.not.allocated(c1_mat)) 
     &        allocate(c1_mat(rows_proj,nprocs*Morb))
            end if            
            if(myid.eq.0) then
             if(.not.allocated(c2_mat)) 
     &        allocate(c2_mat(rows_proj+mod(dimL_orb/2,nprocs),
     &                       nprocs*Morb))
            else 
             if(.not.allocated(c2_mat)) 
     &        allocate(c2_mat(rows_proj,nprocs*Morb))
            end if            
 

      do while(i.le.ND-nprocs+1)
!      do i=start_from_upper, ND

!     Orbital part 
      call CPU_TIME(start)

!     distribute left_Loo parts to different PEs
      if(left_Loo_prec.eqv..TRUE.) then
         left_Loo_hilf2=0d0

         if(.not.allocated(Loo_row))
     .     allocate(Loo_row(dimL_orb))

!              if(myid.eq.0) left_Loo_hilf2(j,:)=
!     .                        left_Loo(i+(j-1)*ND,:)

        if(myid.eq.0) then 
         do j=1, Morb
            Loo_row=0d0
            m=Loo_ind_csr(i+(j-1)*ND)
            if(m.ne.Loo_ind_csr(i+(j-1)*ND+1)) then
            pp=i+(j-1)*ND+1
            do m=Loo_ind_csr(i+(j-1)*ND), Loo_ind_csr(pp)-1 
!                if(abs(Loc_vals_csr(k)).gt.tol1) then
                 Loo_row(Loo_col_csr(m))=Loo_vals_csr(m)
!                end if 
            end do 
            end if

              left_Loo_hilf2(j,:)=
     .                        Loo_row(:)
         end do
        end if
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

         do k=2, nprocs
           left_Loo_hilf2_p0=0d0
           do j=1, Morb
              if(myid.eq.0) then 
               Loo_row=0d0
               ll=Loo_ind_csr(i+k-1+(j-1)*ND)
               if(ll.ne.Loo_ind_csr(i+k-1+(j-1)*ND+1)) then
               pp=i+k-1+(j-1)*ND+1
               do ll=Loo_ind_csr(i+k-1+(j-1)*ND),Loo_ind_csr(pp)-1 
!                   if(abs(Loc_vals_csr(ll)).gt.tol1) then
                    Loo_row(Loo_col_csr(ll))=Loo_vals_csr(ll)
!                   end if
               end do 
               end if
                  left_Loo_hilf2_p0(j,:)=
     .                        Loo_row(:)
              end if 
!              if(myid.eq.0) left_Loo_hilf2_p0(j,:)=
!     .                        left_Loo(i+k-1+(j-1)*ND,:)

              if(myid.eq.k-1) CALL MPI_RECV(left_Loo_hilf2_p0(j,:),
     .         dimL_orb,MPI_DOUBLE_COMPLEX,0,0,
     .         MPI_COMM_WORLD,s_status,ierr) 
              if(myid.eq.0) CALL MPI_SEND(left_Loo_hilf2_p0(j,:),
     .         dimL_orb,MPI_DOUBLE_COMPLEX,k-1,0,MPI_COMM_WORLD,ierr)
!              CALL MPI_BCAST(left_Loo_hilf2_p0(j,:),dimL_orb,
!     .                   MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr) 
              if(myid.eq.k-1) left_Loo_hilf2(j,:)=
     .                        left_Loo_hilf2_p0(j,:)
           end do
           call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         end do
      end if

!     DISTRIBUTE KSL entries to different PEs
             IF(left_Loo_prec.eqv..FALSE..and.KSL_prec.eqv..TRUE.) then

               if(.not.allocated(KSL_hilf2_p0)) 
     .               allocate(KSL_hilf2_p0(ND))
               if(.not.allocated(KSL_conjg_hilf2_p0)) 
     .               allocate(KSL_conjg_hilf2_p0(ND))


              KSL_dist=0d0
              KSL_conjg_dist=0d0

             IF(myid.eq.0) then
              pos=mod(i,ND)
              if(pos.eq.0) pos=ND
              do mm=1, Morb
               do ll=1, Morb
                  k=KSL_mat_ind_csc(pos,mm,ll)
                  if(k.ne.KSL_mat_ind_csc(pos+1,mm,ll)) then
                  pp=pos+1
                  do k=KSL_mat_ind_csc(pos,mm,ll), 
     .                          KSL_mat_ind_csc(pp,mm,ll)-1 
                    if(abs(KSL_mat_vals_csc(k,mm,ll)).gt.tol1) then
                      KSL_dist(KSL_mat_rows_csc(k,mm,ll),mm,ll)=
     .                           KSL_mat_vals_csc(k,mm,ll)
                    end if
                  end do 
                  end if
 
                  k=KSL_conjgmat_ind_csc(pos,mm,ll)
                  if(k.ne.KSL_conjgmat_ind_csc(pos+1,mm,ll)) then
                  pp=pos+1
                  do k=KSL_conjgmat_ind_csc(pos,mm,ll), 
     .                             KSL_conjgmat_ind_csc(pp,mm,ll)-1 
                    if(abs(KSL_conjgmat_vals_csc(k,mm,ll)).gt.tol1) then
                   KSL_conjg_dist(KSL_conjgmat_rows_csc(k,mm,ll),mm,ll)=
     .                           KSL_conjgmat_vals_csc(k,mm,ll)
                    end if
                  end do 
                  end if
!                   KSL_dist(:,mm,ll)=KSL_mat(:,pos,mm,ll)  
!                   KSL_conjg_dist(:,mm,ll)=
!     .                  KSL_conjg_mat(:,pos,mm,ll) 
               end do
              end do 
             END IF

             do m=2, nprocs
              call MPI_BARRIER(MPI_COMM_WORLD,ierr)
              pos2=mod(i+m-1,ND)
              if(pos2.eq.0) pos2=ND
              do mm=1, Morb
                do ll=1, Morb
!                  if(myid.eq.0) KSL_hilf2_p0(:)= KSL_mat(:,pos2,mm,ll)
!                  if(myid.eq.0) KSL_conjg_hilf2_p0(:)= 
!     .                             KSL_conjg_mat(:,pos2,mm,ll)
              KSL_hilf2_p0=0d0
              KSL_conjg_hilf2_p0=0d0
                 if(myid.eq.0) then
                  k=KSL_mat_ind_csc(pos2,mm,ll)
                  if(k.ne.KSL_mat_ind_csc(pos2+1,mm,ll)) then
                  pp=pos2+1
                  do k=KSL_mat_ind_csc(pos2,mm,ll), 
     .                        KSL_mat_ind_csc(pp,mm,ll)-1 
                    if(abs(KSL_mat_vals_csc(k,mm,ll)).gt.tol1) then
                     KSL_hilf2_p0(KSL_mat_rows_csc(k,mm,ll))=
     .                           KSL_mat_vals_csc(k,mm,ll)
                    end if 
                  end do 
                  end if 

                  k=KSL_conjgmat_ind_csc(pos2,mm,ll)
                  if(k.ne.KSL_conjgmat_ind_csc(pos2+1,mm,ll)) then
                  pp=pos2+1
                  do k=KSL_conjgmat_ind_csc(pos2,mm,ll), 
     .                             KSL_conjgmat_ind_csc(pp,mm,ll)-1 
                    if(abs(KSL_conjgmat_vals_csc(k,mm,ll)).gt.tol1) then
                     KSL_conjg_hilf2_p0(KSL_conjgmat_rows_csc(k,mm,ll))=
     .                           KSL_conjgmat_vals_csc(k,mm,ll)
                    end if
                  end do 
                  end if
                 end if

                  if(myid.eq.m-1) CALL MPI_RECV(KSL_hilf2_p0(:),ND,
     .             MPI_DOUBLE_COMPLEX,0,0,MPI_COMM_WORLD,s_status,ierr) 
                  if(myid.eq.0) CALL MPI_SEND(KSL_hilf2_p0(:),ND,
     .             MPI_DOUBLE_COMPLEX,m-1,0,MPI_COMM_WORLD,ierr)
                  if(myid.eq.m-1) CALL MPI_RECV(KSL_conjg_hilf2_p0(:),ND
     .             ,MPI_DOUBLE_COMPLEX,0,1,MPI_COMM_WORLD,s_status,ierr)
                  if(myid.eq.0) CALL MPI_SEND(KSL_conjg_hilf2_p0(:),ND,
     .             MPI_DOUBLE_COMPLEX,m-1,1,MPI_COMM_WORLD,ierr)
                  if(myid.eq.m-1) KSL_dist(:,mm,ll)=KSL_hilf2_p0(:)  
                  if(myid.eq.m-1) KSL_conjg_dist(:,mm,ll)=
     .                  KSL_conjg_hilf2_p0(:) 
                end do
              end do
             end do

             END IF  
    
      L_orb_rows=0d0
      call construct_full_row_orb(PSI,L_orb_rows,i+myid,i)
      call CPU_TIME(finish)
      write(*,*) "orbital row",finish-start


        if(myid.eq.0) then

          L_orb_rows_2=0d0; 
            L_orb_rows_2=L_orb_rows

          do j=1, Morb
           L_row_top(j,1:dimL_orb)=L_orb_rows_2(j,:)

           L_row_top(Morb+j,1:dimL_orb/2)=-Conjg(L_orb_rows_2(j,
     .                                dimL_orb/2+1:dimL_orb))
          L_row_top(Morb+j,dimL_orb/2+1:dimL_orb)=-Conjg(L_orb_rows_2(j,
     .                      1:dimL_orb/2 ))


           if(myid.eq.0) then
            do q=1, dimL_orb
              if(abs(L_row_top(j,q)).ge.tol1) then
               if(binary.eqv..FALSE.) write(256,'(2I10,10x,2F20.10)') 
     .                   i+(j-1)*ND,q,
     .                   Real(L_row_top(j,q)),
     .                   Dimag(L_row_top(j,q))
               if(binary.eqv..TRUE.) 
     .           write(260) i+(j-1)*ND,q,L_row_top(j,q)
                element_count=element_count+1
              end if
              if(abs(L_row_top(j+Morb,q)).ge.tol1) then
                if(binary.eqv..FALSE.) 
     .           write(256,'(2I10,10x,2F20.10)') i+(Morb+j-1)*ND,q,
     .                   Real(L_row_top(j+Morb,q)),
     .                   Dimag(L_row_top(j+Morb,q))
                if(binary.eqv..TRUE.) 
     .           write(260) i+(Morb+j-1)*ND,q,L_row_top(j+Morb,q)
                element_count=element_count+1
              end if
            end do 
           end if 
          end do     

        end if
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         


        do k=2, numprocs  

       L_orb_rows_2=0d0; 

          if(myid.eq.k-1) then
            L_orb_rows_2=L_orb_rows
          end if         

          if(myid.eq.0) CALL MPI_RECV(L_orb_rows_2,
     .         Morb*dimL_orb,MPI_DOUBLE_COMPLEX,k-1,0,
     .         MPI_COMM_WORLD,s_status,ierr) 
          if(myid.eq.k-1) CALL MPI_SEND(L_orb_rows_2,Morb*dimL_orb,
     .         MPI_DOUBLE_COMPLEX,0,0,MPI_COMM_WORLD,ierr)
!       call MPI_BCAST(L_orb_rows_2,Morb*dimL_orb,MPI_DOUBLE_COMPLEX,k-1,
!     .           MPI_COMM_WORLD,ierr)
!       call MPI_BCAST(OCmat_rows_2,Morb*2*Nconf,MPI_DOUBLE_COMPLEX,k-1,
!     .           MPI_COMM_WORLD,ierr)


          if(myid.eq.0) then
           L_row_top=0d0 
           do j=1, Morb
           L_row_top(j,1:dimL_orb)=L_orb_rows_2(j,:)

           L_row_top(Morb+j,1:dimL_orb/2)=-Conjg(L_orb_rows_2(j,
     .                                dimL_orb/2+1:dimL_orb))
          L_row_top(Morb+j,dimL_orb/2+1:dimL_orb)=-Conjg(L_orb_rows_2(j,
     .                      1:dimL_orb/2 ))


            do q=1, dimL_orb
              if(abs(L_row_top(j,q)).ge.tol1) then
               if(binary.eqv..FALSE.)
     .            write(256,'(2I10,10x,2F20.10)') i+(k-1)+(j-1)*ND,q,
     .                   Real(L_row_top(j,q)),
     .                   Dimag(L_row_top(j,q))
               if(binary.eqv..TRUE.) 
     .            write(260) i+(k-1)+(j-1)*ND,q,
     .                   L_row_top(j,q)
                element_count=element_count+1
              end if
              if(abs(L_row_top(j+Morb,q)).ge.tol1) then
               if(binary.eqv..FALSE.)
     .          write(256,'(2I10,10x,2F20.10)') i+(k-1)+(Morb+j-1)*ND,q,
     .                   Real(L_row_top(j+Morb,q)),
     .                   Dimag(L_row_top(j+Morb,q))
               if(binary.eqv..TRUE.)
     .          write(260) i+(k-1)+(Morb+j-1)*ND,q,
     .                   L_row_top(j+Morb,q)
                element_count=element_count+1
              end if
            end do 

           end do     
          end if 
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        end do 


        write(*,*) "row",i+myid," constructed" 

 
       i=i+nprocs
      end do

!    REMAINING rows only by root process
!      if(i.gt.ND-nprocs+1.and.myid.eq.0.and.i.le.ND) then
      if(i.gt.ND-nprocs+1.and.i.le.ND) then
        do while(i.le.ND)

!        Orbital  
         call CPU_TIME(start)


!     store needed left_Loo parts in variable
      IF (myid.eq.0) THEN
      if(left_Loo_prec.eqv..TRUE.) then
         left_Loo_hilf2=0d0
!         do k=1, nprocs
           do j=1, Morb

            Loo_row=0d0
            m=Loo_ind_csr(i+(j-1)*ND)
            if(m.ne.Loo_ind_csr(i+(j-1)*ND+1)) then
            pp=i+(j-1)*ND+1
            do m=Loo_ind_csr(i+(j-1)*ND), Loo_ind_csr(pp)-1 
!                if(abs(Loc_vals_csr(k)).gt.tol1) then
                 Loo_row(Loo_col_csr(m))=Loo_vals_csr(m)
!                end if 
            end do 
            end if

              left_Loo_hilf2(j,:)=
     .                        Loo_row(:)

!              left_Loo_hilf2(j,:)=
!     .                        left_Loo(i+(j-1)*ND,:)
           end do
!         end do
      end if


!        PUT KSL-mat entries in KSL_dist if preconstructed
           if(left_Loo_prec.eqv..FALSE..and.KSL_prec.eqv..TRUE.) then
              pos=mod(i,ND)
              if(pos.eq.0) pos=ND
                 do mm=1, Morb
                   do ll=1, Morb
                    k=KSL_mat_ind_csc(pos,mm,ll)
                    if(k.ne.KSL_mat_ind_csc(pos+1,mm,ll)) then
                    pp=pos+1
                    do k=KSL_mat_ind_csc(pos,mm,ll), 
     .                          KSL_mat_ind_csc(pp,mm,ll)-1 
                      if(abs(KSL_mat_vals_csc(k,mm,ll)).gt.tol1) then
                       KSL_dist(KSL_mat_rows_csc(k,mm,ll),mm,ll)=
     .                           KSL_mat_vals_csc(k,mm,ll)
                      end if
                    end do 
                    end if

                    k=KSL_conjgmat_ind_csc(pos,mm,ll)
                    if(k.ne.KSL_conjgmat_ind_csc(pos+1,mm,ll)) then
                    pp=pos+1
                    do k=KSL_conjgmat_ind_csc(pos,mm,ll), 
     .                             KSL_conjgmat_ind_csc(pp,mm,ll)-1 
                    if(abs(KSL_conjgmat_vals_csc(k,mm,ll)).gt.tol1) then
                   KSL_conjg_dist(KSL_conjgmat_rows_csc(k,mm,ll),mm,ll)=
     .                           KSL_conjgmat_vals_csc(k,mm,ll)
                    end if
                    end do 
                    end if
!                       KSL_dist(:,mm,ll)=  
!     .                        KSL_mat(:,i,mm,ll)
!                       KSL_conjg_dist(:,mm,ll)= 
!     .                        KSL_conjg_mat(:,i,mm,ll)
                   end do
                 end do 
           end if
         END IF
         
         L_orb_rows=0d0
         call construct_full_row_orb(PSI,L_orb_rows,i,i)
         call CPU_TIME(finish)
         if(myid.eq.0) write(*,*) "orbital row",finish-start


         IF(myid.eq.0) THEN
          do j=1, Morb
           L_row_top(j,1:dimL_orb)=L_orb_rows(j,:)

           L_row_top(Morb+j,1:dimL_orb/2)=-Conjg(L_orb_rows(j,
     .                                dimL_orb/2+1:dimL_orb))
           L_row_top(Morb+j,dimL_orb/2+1:dimL_orb)=-Conjg(L_orb_rows(j,
     .                      1:dimL_orb/2 ))


            do q=1, Nconf+dimL_orb
              if(abs(L_row_top(j,q)).ge.tol1) then
               if(binary.eqv..FALSE.)
     .             write(256,'(2I10,10x,2F20.10)') i+(j-1)*ND,q,
     .                   Real(L_row_top(j,q)),
     .                   Dimag(L_row_top(j,q))
               if(binary.eqv..TRUE.) 
     .           write(260) i+(j-1)*ND,q,L_row_top(j,q)
                element_count=element_count+1
              end if
              if(abs(L_row_top(j+Morb,q)).ge.tol1) then
                if(binary.eqv..FALSE.) 
     .             write(256,'(2I10,10x,2F20.10)') i+(Morb+j-1)*ND,q,
     .                   Real(L_row_top(j+Morb,q)),
     .                   Dimag(L_row_top(j+Morb,q))
                if(binary.eqv..TRUE.) 
     .           write(260) i+(Morb+j-1)*ND,q,L_row_top(j+Morb,q)
                element_count=element_count+1
              end if
            end do 

          end do     
         END IF
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)

          if(myid.eq.0) write(*,*) "row",i," constructed" 
          i=i+1 
        end do
      end if
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

       
!      Now write Loc matrix to file
       IF(myid.eq.0) THEN
          do q=1, Loc_nonzeroes
               if(binary.eqv..FALSE.)
     .             write(256,'(2I10,10x,2F20.10)') Loc_rows_hilf(q),
     .                   Loc_cols_hilf(q)+dimL_orb,  
     .                   Real(Loc_vals_hilf(q)),
     .                   Dimag(Loc_vals_hilf(q))
               if(binary.eqv..TRUE.) 
     .           write(260) Loc_rows_hilf(q),dimL_orb+Loc_cols_hilf(q),
     .                      Loc_vals_hilf(q)
            
            if(Loc_cols_hilf(q).le.NC) then
               if(binary.eqv..FALSE.)
     .             write(256,'(2I10,10x,2F20.10)')
     .                   Loc_rows_hilf(q)+dimL_orb/2,
     .                   Loc_cols_hilf(q)+dimL_orb+NC,  
     .                   -1.0*Real(Loc_vals_hilf(q)),
     .                   Dimag(Loc_vals_hilf(q))
               if(binary.eqv..TRUE.) 
     .           write(260) Loc_rows_hilf(q)+dimL_orb/2,
     .                      dimL_orb+Loc_cols_hilf(q)+NC,
     .                      -1.0*Conjg(Loc_vals_hilf(q))
            else
               if(binary.eqv..FALSE.)
     .             write(256,'(2I10,10x,2F20.10)')
     .                   Loc_rows_hilf(q)+dimL_orb/2,
     .                   Loc_cols_hilf(q)+dimL_orb-NC,  
     .                   -1.0*Real(Loc_vals_hilf(q)),
     .                   Dimag(Loc_vals_hilf(q))
               if(binary.eqv..TRUE.) 
     .           write(260) Loc_rows_hilf(q)+dimL_orb/2,
     .                      dimL_orb+Loc_cols_hilf(q)-NC,
     .                      -1.0*Conjg(Loc_vals_hilf(q))
            end if
          end do
       END IF   
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)


        call CPU_TIME(finish_gen)
        write(*,*) "Upper part constructed in ",
     .                       finish_gen-start_gen," seconds"
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)


      end if  !end if upper part already constructed
!      stop


654     continue
!======= CONSTRUCTING LOWER PART OF LR-MATRIX  ===============================================
      if(Morb.gt.1) then

        call CPU_TIME(start_gen)
     
        CI_Production_1b=.FALSE.
        CI_Production_2b=.FALSE.
        CI_Prd=.FALSE.

        allocate(CI_row(2*Nconf,nprocs)) 

        bare=.FALSE.
        i=start_from_lower
!        do i=start_from_lower, Nconf 
        do while(i.le.Nconf-nprocs+1)
          call construct_CO_matrix_row(PSI,VIN,
     .       bare,COmat_row,
     .       COmat_bare_row,Nconf,i+myid)


          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
          do k=1, nprocs
            call construct_CI_row(CI_row(:,k),Nconf,i+k-1)
!           BCAST unnötig da nur auf p=0 benötigt  
!            call MPI_BCAST(CI_row(:,k),2*Nconf,MPI_DOUBLE_COMPLEX,0,     
!     .           MPI_COMM_WORLD,ierr)
          end do
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)

          do k=2, nprocs
           if(myid.eq.k-1) CALL MPI_RECV(CI_row(:,k),
     .         2*Nconf,MPI_DOUBLE_COMPLEX,0,0,
     .         MPI_COMM_WORLD,s_status,ierr) 
           if(myid.eq.0) CALL MPI_SEND(CI_row(:,k),2*Nconf,
     .         MPI_DOUBLE_COMPLEX,k-1,0,MPI_COMM_WORLD,ierr)
          end do     


          L_row_bottom(1,1:dimL_orb)=COmat_row
          L_row_bottom(1,dimL_orb+1:dimL_orb+2*Nconf)=CI_row(:,myid+1)

          L_row_bottom(2,1:dimL_orb/2)=-Conjg(COmat_row(dimL_orb/2+1
     .                                   :dimL_orb))
          L_row_bottom(2,dimL_orb/2+1:dimL_orb)=
     .                       -Conjg(COmat_row(1:dimL_orb/2))
          L_row_bottom(2,dimL_orb+1:dimL_orb+Nconf)=
     .                       -Conjg(CI_row(NC+1:2*NC,myid+1))
          L_row_bottom(2,dimL_orb+1+Nconf:dimL_orb+2*NC)=
     .                       -Conjg(CI_row(1:NC,myid+1))


         if(myid.eq.0) then 
             do j=1, 2*Nconf+dimL_orb

               if(abs(L_row_bottom(1,j)).ge.tol1) then
                 if(binary.eqv..FALSE.)
     .              write(256,'(2I10,10x,2F22.16)') i+dimL_orb,j,
     .                     Dreal(L_row_bottom(1,j)),
     .                     Dimag(L_row_bottom(1,j))
                 if(binary.eqv..TRUE.)
     .              write(260) i+dimL_orb,j,
     .                     L_row_bottom(1,j)
                element_count=element_count+1
               end if
               if(abs(L_row_bottom(2,j)).ge.tol1) then
                 if(binary.eqv..FALSE.)
     .             write(256,'(2I10,10x,2F20.10)') i+dimL_orb+NC,j,
     .                     Real(L_row_bottom(2,j)),
     .                     Dimag(L_row_bottom(2,j))
                 if(binary.eqv..TRUE.)
     .              write(260) i+dimL_orb+NC,j,
     .                     L_row_bottom(2,j)
                element_count=element_count+1
               end if

             end do 
         end if


         do k=2, nprocs

          L_row_bottom_2=0d0 
          if(myid.eq.k-1) then
            L_row_bottom_2=L_row_bottom
          end if

          if(myid.eq.0) CALL MPI_RECV(L_row_bottom_2,
     .         2*dimL,MPI_DOUBLE_COMPLEX,k-1,0,
     .         MPI_COMM_WORLD,s_status,ierr) 
          if(myid.eq.k-1) CALL MPI_SEND(L_row_bottom_2,2*dimL,
     .         MPI_DOUBLE_COMPLEX,0,0,MPI_COMM_WORLD,ierr)
          

!        if(myid.eq.0) then
!          L_row_bottom(1,1:dimL_orb)=COmat_row_2
!          L_row_bottom(1,dimL_orb+1:dimL_orb+2*Nconf)=CI_row(:,myid+1)

!          L_row_bottom(2,1:dimL_orb/2)=-Conjg(COmat_row_2(dimL_orb/2+1
!     .                                   :dimL_orb))
!          L_row_bottom(2,dimL_orb/2+1:dimL_orb)=
!     .                       -Conjg(COmat_row_2(1:dimL_orb/2))
!          L_row_bottom(2,dimL_orb+1:dimL_orb+Nconf)=
!     .                       -Conjg(CI_row(NC+1:2*NC,myid+1))
!          L_row_bottom(2,dimL_orb+1+Nconf:dimL_orb+2*NC)=
!     .                       -Conjg(CI_row(1:NC,myid+1))
!        end if      



!         L_row_bottom_2=0d0
!         if(myid.eq.k-1) then 
!            L_row_bottom_2=L_row_bottom
!         end if  
!         call MPI_BCAST(L_row_bottom_2,2*dimL,MPI_DOUBLE_COMPLEX,k-1,
!     .           MPI_COMM_WORLD,ierr)
       
           if(myid.eq.0) then
             do j=1, 2*Nconf+dimL_orb

               if(abs(L_row_bottom_2(1,j)).ge.tol1) then
                if(binary.eqv..FALSE.) 
     .             write(256,'(2I10,10x,2F20.10)') i+(k-1)+dimL_orb,j,
     .                     Real(L_row_bottom_2(1,j)),
     .                     Dimag(L_row_bottom_2(1,j))
                if(binary.eqv..TRUE.) 
     .             write(260) i+(k-1)+dimL_orb,j,
     .                     L_row_bottom_2(1,j)
                element_count=element_count+1
               end if
               if(abs(L_row_bottom_2(2,j)).ge.tol1) then
                if(binary.eqv..FALSE.) 
     .            write(256,'(2I10,10x,2F20.10)') i+(k-1)+dimL_orb+NC,j,
     .                     Real(L_row_bottom_2(2,j)),
     .                     Dimag(L_row_bottom_2(2,j))
                if(binary.eqv..TRUE.) 
     .             write(260) i+(k-1)+dimL_orb+NC,j,
     .                     L_row_bottom_2(2,j)
                element_count=element_count+1
               end if

             end do 
           end if
         end do


         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         write(*,*) "lower part row",i+myid," constructed" 

         i=i+nprocs
        end do 


        if(i.gt.NC-nprocs+1.and.i.le.NC) then

          do while(i.le.NC)


          call construct_CO_matrix_row(PSI,VIN,
     .       bare,COmat_row,
     .       COmat_bare_row,Nconf,i)


          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            call construct_CI_row(CI_row(:,1),Nconf,i)
!            call MPI_BCAST(CI_row(:,1),2*Nconf,MPI_DOUBLE_COMPLEX,0,
!     .           MPI_COMM_WORLD,ierr)
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)

         if(myid.eq.0) then
          L_row_bottom(1,1:dimL_orb)=COmat_row
          L_row_bottom(1,dimL_orb+1:dimL_orb+2*Nconf)=CI_row(:,1)

          L_row_bottom(2,1:dimL_orb/2)=-Conjg(COmat_row(dimL_orb/2+1
     .                                   :dimL_orb))
          L_row_bottom(2,dimL_orb/2+1:dimL_orb)=
     .                       -Conjg(COmat_row(1:dimL_orb/2))
          L_row_bottom(2,dimL_orb+1:dimL_orb+Nconf)=
     .                       -Conjg(CI_row(NC+1:2*NC,1))
          L_row_bottom(2,dimL_orb+1+Nconf:dimL_orb+2*NC)=
     .                       -Conjg(CI_row(1:NC,1))

       
             do j=1, 2*Nconf+dimL_orb

               if(abs(L_row_bottom(1,j)).ge.tol1) then
                 if(binary.eqv..FALSE.)
     .              write(256,'(2I10,10x,2F20.10)') i+dimL_orb,j,
     .                     Real(L_row_bottom(1,j)),
     .                     Dimag(L_row_bottom(1,j))
                 if(binary.eqv..TRUE.)
     .              write(260) i+dimL_orb,j,
     .                     L_row_bottom(1,j)
                element_count=element_count+1
               end if
               if(abs(L_row_bottom(2,j)).ge.tol1) then
                 if(binary.eqv..FALSE.)
     .              write(256,'(2I10,10x,2F20.10)') i+dimL_orb+NC,j,
     .                     Real(L_row_bottom(2,j)),
     .                     Dimag(L_row_bottom(2,j))
                 if(binary.eqv..TRUE.)
     .              write(260) i+dimL_orb+NC,j,
     .                     L_row_bottom(2,j)
                element_count=element_count+1
               end if

             end do 
         end if


            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            if(myid.eq.0) write(*,*) "lower part row",i," constructed" 
            i=i+1
          end do  

        end if 

      end if  !! if Morb gt 1 

        call CPU_TIME(finish_gen)
        write(*,*) "Lower part constructed in ",
     .             finish_gen-start_gen," seconds"
     
     
      if(binary.eqv..FALSE..and.myid.eq.0) close(256)
      if(binary.eqv..TRUE..and.myid.eq.0) close(260)
      if(myid.eq.0) close(257)


 
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)  

!      end if

      call CPU_TIME(finish_tot)
     
      if(myid.eq.0) then 
        write(*,*) 
        write(*,*) "LR-MATRIX constructed in",finish_tot-start_tot,
     .                                             "seconds" 
       if(binary.eqv..FALSE.) 
     .    write(*,*) "It is stored in 'LR_matrix_save.dat'."
       if(binary.eqv..TRUE.) 
     .    write(*,*) "It is stored in 'LR_bin'."
      end if
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)  


!================================================================================================
!=============== END OF MATRIX CONSTRUCTION =====================================================
!================================================================================================


9000    continue

         if(myid.eq.0) then
             write(6,*)" Cleaning some intermediate files:"
               file2='Tkin_2D.dat'
           call remove(file2)
               file2='left_Loc.dat'
           call remove(file2)
               file2='left_Loc_inter.dat'
           call remove(file2)
               file2='left_Loo.dat'
           call remove(file2)
               file2='left_Loo_inter.dat'
           call remove(file2)
               file2='V_PrdCIJ.dat'
           call remove(file2)
               file2='V_PrdCIJ2.dat'
           call remove(file2)
               file2='KSL_mat.dat'
           call remove(file2)
               file2='KSL_conjg_mat.dat'
           call remove(file2)
               file2='proj.dat'
           call remove(file2)
         end if
c      
      end subroutine LR_ARNOLDI
