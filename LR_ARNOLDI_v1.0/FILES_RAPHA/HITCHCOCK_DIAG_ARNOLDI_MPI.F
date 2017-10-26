      subroutine DIAG_ARNOLDI_MPI(PSI,VIN) 
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
                  
      integer*8           ldv, int1, int2, maxdim_prec, ldz
      parameter         ( 
     &                     maxdim_prec=50000)
      integer             file_length
      parameter           (file_length=1461582)                   
c
c     %--------------%
c     | Local Arrays |
c     %--------------%
c
      integer           iparam(11), ipntr(14), mjob(6)
      logical           select(maxncv)
      logical,allocatable ::        sel(:)
      logical           wantt,wantz,tempp_log,is_real,cu_analysis 
      logical             ::   squared
      Complex*16, allocatable ::
     &                  ax(:), d(:), 
     &                  v(:,:), workd(:), 
     &                  workev(:), resid(:), 
     &                  workl(:),
     &                  vec(:) 
      Complex*16, allocatable ::  matrix(:,:)
      Complex*16, allocatable ::  matrix2(:,:)
      Complex*16, allocatable ::  v_final(:,:)
      Complex*16, allocatable ::  v_inter(:,:)
      Complex*16, allocatable ::  mat_val(:)
      Complex*16, allocatable ::  mat_val_up(:)
      Complex*16, allocatable ::  mat_val_down(:)
      Complex*16, allocatable ::  mat_val2(:)
      Complex*16, allocatable ::  vec_new(:)
      Complex*16, allocatable ::  vec_hilf(:)
      Complex*16, allocatable ::  vec2(:)
      Complex*16, allocatable ::  vec2_new(:)
      Complex*16, allocatable ::  resf(:)
      Complex*16, allocatable ::  cu_vec(:)
      Complex*16, allocatable ::  cv_vec(:)
      Complex*16, allocatable ::  cu(:,:)
      Complex*16, allocatable ::  cv(:,:)
      Complex*16, allocatable ::  cu_temp(:)
      Complex*16, allocatable ::  cv_temp(:)
      Complex*16, allocatable ::  rij_eigen(:)
      Complex*16, allocatable ::  u_vec(:,:)
      Complex*16, allocatable ::  v_vec(:,:)
      Complex*16, allocatable ::  CCu(:,:,:)
      Complex*16, allocatable ::  CCv(:,:,:)
      Complex*16, allocatable ::  t(:)
      Complex*16, allocatable ::  rij_temp1(:,:)
      Complex*16, allocatable ::  rij_temp2(:,:)
      Complex*16, allocatable ::  rij_mat(:,:)
      Complex*16, allocatable ::  work_eigen(:)
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
      Complex*16, allocatable ::  dens_resp_orb(:,:)
      Complex*16, allocatable ::  dens_resp_CI(:,:)
      Complex*16, allocatable ::  Rho(:,:)
      Complex*16, allocatable ::  RhoSQR(:,:)
      Complex*16, allocatable ::  RhoSQR1(:,:)
      Complex*16, allocatable ::  RhoSQR2(:,:)
      Complex*16, allocatable ::  Rho_InvSqr(:,:)
      Complex*16, allocatable ::  Comp_Rho_InvSqr(:,:)
      Complex*16, DIMENSION(NDX*NDY*NDZ,Morb), INTENT(IN) ::  PSI
      Complex*16, DIMENSION(Nconf), INTENT(IN) ::  VIN
      Complex*16, allocatable ::  Vhelp(:)
      Complex*16, allocatable ::  VIN_help(:)
      Double precision, allocatable ::  rwork_hess(:)
      Double precision, allocatable ::  rwork_eigen(:)
      Complex*16, allocatable   ::  mat_val_csr(:)
      Complex*16, allocatable   ::  mat_val_csr_up(:)
      Complex*16, allocatable   ::  mat_val_csr_down(:)

      COMPLEX*16, allocatable :: ut(:,:,:), vt(:,:,:)
      COMPLEX*16, allocatable :: u2(:,:,:), v2(:,:,:)
      COMPLEX*16, allocatable :: tr(:,:)
      COMPLEX*16, allocatable :: weight_orb(:,:), weight_CI(:,:)
      COMPLEX*16, allocatable :: CCu_r(:,:,:)
      COMPLEX*16, allocatable :: CCv_r(:,:,:)
      COMPLEX*16, allocatable :: orb_int(:,:,:)

      integer*4, allocatable    ::  mat_csr1(:)
      integer*4, allocatable    ::  mat_csr1_up(:)
      integer*4, allocatable    ::  mat_csr1_down(:)
      integer*4, allocatable    ::  mat_csr2(:)
      integer*4, allocatable    ::  mat_csr2_up(:)
      integer*4, allocatable    ::  mat_csr2_down(:)
      integer*4, allocatable    ::  mat_coo1(:)
      integer*4, allocatable    ::  mat_coo2(:)
      integer*4, allocatable    ::  mat_coo1_up(:)
      integer*4, allocatable    ::  mat_coo2_up(:)
      integer*4, allocatable    ::  mat_coo1_down(:)
      integer*4, allocatable    ::  mat_coo2_down(:)
      integer, allocatable    ::  ifaill(:)
      integer, allocatable    ::  ifailr(:)
      logical,allocatable     ::  sel_temp(:)
      Double precision 
     &                  rwork(maxncv), rd(maxncv,3)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      character         bmat*1, which*2, howmny*1
      character(50)     file1,file2
      character(20)     path,pathORB,pathCI,pathType,pathRes,pathRij
      parameter         (file1='LR_MATRIX_Original.dat')
      integer           myid, numprocs, ido, n, nx, lworkl, info, j,
     &                  ierr, nconv, maxitr, ishfts, mode, i, k, l,
     &                  status_open, rows, sum, iter, iter2,infoo,
     &                  conv_count,lwork, ll, s ,q, ndigit, logfil,
     &                  mcaupd, counter, ii, kk, jj, ss,
     &                  last_row, start_from_upper, start_from_lower,
     &                  start_from_i, start_from_j, start_counter,
     &                  nonzero, lll, qq, counter_up,
     &                  counter_down, nonzero_up, nonzero_down,
     &                  nloc_up, nloc_down, rows_up, rows_down, ind, 
     &                  s_status(MPI_STATUS_SIZE), summe, sum1, sum2, st
      integer*8         NC 
      Complex*16
     &                  sigma,tempp,alpha,beta
      Real*8
     &                  tol, a, b, ticks, start, finish,time,
     &                  t2, t1, a2, b2, a3, b3, sum_u, sum_v
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
!      logical :: eof


      if(Morb.gt.1) dimL = 2*Morb*(NDX*NDY*NDZ) + 2*Nconf
      if(Morb.eq.1) dimL = 2*Morb*(NDX*NDY*NDZ) 
      dimL_orb = 2*Morb*NDX*NDY*NDZ 
      NC = Nconf
      ND = NDX*NDY*NDZ 

!      ldv=maxn
      ldv=dimL


!       allocate(ax(maxn), d(maxncv),  
!     &          v(ldv,maxncv), workd(3*maxn), 
!     &          workev(3*maxncv), resid(maxn), 
!     &          workl(3*maxncv*maxncv+5*maxncv),
!     &          vec(maxn)) 

       allocate(ax(maxn), !d(maxncv),  
     &          workd(3*dimL), !v(ldv,ncv),
     &          workev(2*ncv), resid(maxn), 
     &          workl(3*ncv*ncv+5*ncv),
     &          vec(maxn)) 

!      if(.not.allocated(PSI)) allocate(PSI(NDX*NDY*NDZ,Morb)) 
!      if(.not.allocated(VIN)) allocate(VIN(Nconf))
      allocate(Vhelp(Nconf),VIN_help(Nconf)) 

c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c 

c      call CPU_TIME(start)

c      call MPI_INIT( ierr )
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
!      nev   = 250
!      ncv   = 500
!      min_conv = 120 
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
      tol    = 1.0e-10 
!      tol1   = 1.0e-12 
      tol2   = 1.0e-10 
!      upperlimit = 0.5e1
      ido    = 0

      restart=.FALSE.
      if(restart.eqv..FALSE.) then
        info   = 0  !random initial resid vector
      else

        info   = 1  !user-defined initial resid vector
        open(unit=216, file='residual_vec.dat', 
     .      status='old', action='READ', iostat=status_open)
          do j=1, dimL
             read(216,*) a,b
             resid(j)=cmplx(a,b,kind=8)
          end do
        close(216)
       

!          open(unit=218, file='arnoldi_basis.dat', 
!     .      status='old', action='READ', iostat=status_open)
!           do while (.not.eof(218)) 
!             read(218,*) i,j,a,b
!             v(j,i)=cmplx(a,b,kind=8)
!           end do
!          close(218)



!====  Cu/Cv analysis section  =========
       if(myid.eq.0) then
      
       cu_analysis=.FALSE. 
!       CI_Production_1b=.FALSE.
         
       if(cu_analysis.eqv..TRUE.) then

         LWORK=2*Morb-1
         allocate(cu_temp(NC),cv_temp(NC))
         allocate(rij_temp1(Morb,Morb),rij_temp2(Morb,Morb))
         allocate(rij_mat(Morb,Morb),rij_eigen(Morb))
         allocate(rwork_eigen(3*Morb-2),work_eigen(Lwork))

         pathRij='Eigenv'//path_sep//'rij'//path_sep
         call mk_dir(pathRij)


         open(unit=148, file=trim(adjustl(pathRij))//
     .         'exc_occupation.dat',
     .         status='replace',action='READWrite',
     .         iostat=status_open)


         open(unit=216, file='Eigenv/typ/type.dat', 
     .      status='old', action='READ', iostat=status_open)
          i=0
          do while(i.le.200)
!             read(216,*,iostat=st) i,a
!             if(st.ne.0) exit
!             if(abs(abs(a)-1.0).le.1e-2) then
              i=i+1
              if(i.le.9) write(file2,'(A26,I1,A4)')
     .                  'Eigenv/CI/eigenvectors_CI0',i,'.dat'
              if(i.gt.9) write(file2,'(A25,I2,A4)') 
     .                  'Eigenv/CI/eigenvectors_CI',i,'.dat'
              if(i.gt.99.and.i.lt.1000) write(file2,'(A25,I3,A4)') 
     .                  'Eigenv/CI/eigenvectors_CI',i,'.dat'

              open(unit=217, file=file2, 
     .          status='old', action='READ', iostat=status_open)
                 j=1
                 do 
                   read(217,*,iostat=st) a2,b2,a3,b3
                   if(st.ne.0) exit 
                   cu_temp(j)=Cmplx(a2,b2,kind=8) 
                   cv_temp(j)=Cmplx(a3,b3,kind=8) 
                  
                   j=j+1
                 end do 
              close(217) 
             
!              call HPSI(cu_temp,Vhelp)
!              write(*,*) ZRIJ
!              stop 
             
              DO k=1,Morb
                DO l=1,Morb

                  VIN_help=cu_temp; Vhelp=0.0d0; 
                  IF(Morb.gt.1) THEN
                   CALL Produce_Cij(VIN_help,Vhelp,k,l,0,0,Nconf)
                  END IF
                  rij_temp1(k,l)=dot_product(cu_temp(:),Vhelp(:))
!                  rij_temp1(k,l)=dot_product(VIN(:),Vhelp(:))

                  VIN_help=cv_temp; Vhelp=0.0d0; 
                  IF(Morb.gt.1) THEN
                   CALL Produce_Cij(VIN_help,Vhelp,k,l,0,0,Nconf)
                  END IF
                  rij_temp2(k,l)=dot_product(cv_temp(:),Vhelp(:))
!                  rij_temp2(k,l)=dot_product(VIN(:),Vhelp(:))

                END DO
              END DO


          rij_mat=rij_temp1
          rij_eigen=0d0
          call zheev('N','U',Morb,rij_mat,Morb,rij_eigen,
     .         work_eigen,LWORK,rwork_eigen,infoo)
     
          write(*,*) i,Real(rij_eigen(1))
          do k=1, Morb
            write(148,*) i,k,Abs(rij_eigen(k))
          end do

              if(i.le.9) write(file2,'(A8,I1,A4)')
     .                  'newRIJ_0',i,'.dat'
              if(i.gt.9) write(file2,'(A7,I2,A4)') 
     .                  'newRIJ_',i,'.dat'
              if(i.gt.99.and.i.lt.1000) write(file2,'(A7,I3,A4)') 
     .                  'newRIJ_',i,'.dat'

              open(unit=145, file=trim(adjustl(pathRij))//file2,
     .         status='replace',action='READWrite',
     .         iostat=status_open)
               do k=1, Morb
                do l=1, Morb  
                  write(145,'(2I5,4E20.10)') k,l,Real(rij_temp1(k,l)),
     .                Dimag(rij_temp1(k,l)),Real(rij_temp2(k,l)),
     .                Dimag(rij_temp2(k,l)) 
                end do
               end do  
              close(145)

!             end if
          end do
         close(216)
         close(148)

         goto 8000

       end if
!====  End CU/Cv analysis
       end if
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)

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

      if(myid.eq.0) then
      write(*,*) "%---------------------------------------------------%"
      write(*,*) "LR-simulation started with ", nprocs,
     .                    " MPI process(es)"
      write(*,*) "Dimension of the LR-matrix is ", nx
      write(*,*) "You ask for",min_conv,"eigenvectors of ", which
      write(*,*) "up to",upperlimit
      write(*,*) "with residual tolerance",tol2
      write(*,*) "# of Krylov vectors is",ncv
      write(*,*) "# of wanted eigenvalues per Arnoldi iteration is",nev
      write(*,*) "%---------------------------------------------------%"
      end if

!==========================================================================================

!===     READ-IN of LR-MATRIX      
!          if(.not.allocated(LR_mat)) allocate(LR_mat(nx,nx),stat=ierr) 
!          if(ierr /= 0) then
!              write(*,*) myid,ierr,"allocation error in LR_MAT"
!              stop 
!          end if

      if(binary.eqv..TRUE.) then
        Open (unit=260, File='LR_bin', Status='old', 
     .  Form='unformatted', Access='sequential')
      else
        open(unit=217, file='LR_matrix_save.dat', 
     .      status='old', action='READ', iostat=status_open)
      end if
!      counter=1
!      do  
!      Read(260,iostat=st) i,j,tempp
!      if(st.ne.0) exit
!      write(*,*) i,j,tempp
!      counter=counter+1
!      end do
!      write(*,*) st, counter
!      stop


!       if(Morb.eq.1) then

        allocate(mat_val(maxnonzero/nprocs),mat_coo1(maxnonzero/nprocs)
     .           ,mat_coo2(maxnonzero/nprocs))
!        allocate(mat_val_csr(maxnonzero/nprocs),mat_csr1(nx+1)
!     .           ,mat_csr2(maxnonzero/nprocs))


          nloc = (nx / nprocs)

           counter=1  
           do  
             if(binary.eqv..FALSE.) read(217,*,iostat=st) i,j,a,b
             if(binary.eqv..TRUE.) Read(260,iostat=st) i,j,tempp
             if(st.ne.0) exit
             if(i.gt.myid*nloc.and.i.le.(myid+1)*nloc
!     .          .and.sqrt(a**2+b**2).gt.tol1) then
     .          ) then
               if(binary.eqv..FALSE.) mat_val(counter)=Cmplx(a,b,kind=8)
               if(binary.eqv..TRUE.) mat_val(counter)=tempp
                mat_coo1(counter)=i-myid*nloc
                mat_coo2(counter)=j  
                counter=counter+1
             end if
             if(myid.eq.nprocs-1.and.i.gt.(myid+1)*nloc
!     .         .and.sqrt(a**2+b**2).gt.tol1) then
     .         ) then
               if(binary.eqv..FALSE.) mat_val(counter)=Cmplx(a,b,kind=8)
               if(binary.eqv..TRUE.) mat_val(counter)=tempp
                mat_coo1(counter)=i-myid*nloc
                mat_coo2(counter)=j  
                counter=counter+1
             end if
           end do
          if(binary.eqv..FALSE.) close(217)
          if(binary.eqv..TRUE.) close(260)
          nonzero=counter-1
          counter=nonzero

        allocate(mat_val_csr(nonzero),mat_csr1(nx+1)
     .           ,mat_csr2(nonzero))

          write(*,*) myid, counter

          call MPI_ALLREDUCE(counter,summe,1,MPI_INTEGER,MPI_SUM,
     .                     MPI_COMM_WORLD,ierr)
!          call MPI_REDUCE(counter,summe,1,MPI_INTEGER,MPI_SUM,0,
!     .                     MPI_COMM_WORLD,ierr)
          if(myid.eq.0) write(*,*) "AllSum:",summe


!       else


!        allocate(mat_val_up(maxnonzero/nprocs),
!     .            mat_coo1_up(maxnonzero/nprocs)
!     .           ,mat_coo2_up(maxnonzero/nprocs))
!        allocate(mat_val_down(maxnonzero/nprocs),
!     .            mat_coo1_down(maxnonzero/nprocs)
!     .           ,mat_coo2_down(maxnonzero/nprocs))


!          nloc_up = dimL_orb/nprocs
!          nloc_down = 2*NC/nprocs


!           counter_up=1  
!           counter_down=1  
!           do  

!             if(binary.eqv..FALSE.) read(217,*,iostat=st) i,j,a,b
!             if(binary.eqv..TRUE.) Read(260,iostat=st) i,j,tempp
!             if(st.ne.0) exit

!             if(i.gt.myid*nloc_up.and.i.le.(myid+1)*nloc_up
!     .           ) then
!               if(binary.eqv..FALSE.) 
!     .               mat_val_up(counter_up)=Cmplx(a,b,kind=8)
!               if(binary.eqv..TRUE.) mat_val_up(counter_up)=tempp
!                mat_coo1_up(counter_up)=i-myid*nloc_up
!                mat_coo2_up(counter_up)=j  
!                counter_up=counter_up+1
!             end if

!             if(myid.eq.nprocs-1.and.i.gt.(myid+1)*nloc_up.and.
!     .          i.le.dimL_orb) then
!               if(binary.eqv..FALSE.) 
!     .            mat_val_up(counter_up)=Cmplx(a,b,kind=8)
!               if(binary.eqv..TRUE.) mat_val_up(counter_up)=tempp
!                  mat_coo1_up(counter_up)=i-myid*nloc_up
!                  mat_coo2_up(counter_up)=j  
!                  counter_up=counter_up+1
!             end if                


!             if(i.gt.dimL_orb+myid*nloc_down.and.
!     .          i.le.dimL_orb+(myid+1)*nloc_down
!     .          ) then
!               if(binary.eqv..FALSE.) 
!     .              mat_val_down(counter_down)=Cmplx(a,b,kind=8)
!               if(binary.eqv..TRUE.) mat_val_down(counter_down)=tempp
!                mat_coo1_down(counter_down)=
!     .                     i-(dimL_orb+myid*nloc_down)
!                mat_coo2_down(counter_down)=j  
!                counter_down=counter_down+1
!             end if


!             if(myid.eq.nprocs-1.and.
!     .               i.gt.dimL_orb+(myid+1)*nloc_down
!     .              ) then
!               if(binary.eqv..FALSE.) 
!     .             mat_val_down(counter_down)=Cmplx(a,b,kind=8)
!               if(binary.eqv..TRUE.) mat_val_down(counter_down)=tempp
!                 mat_coo1_down(counter_down)=i-(dimL_orb+myid*nloc_down)
!                 mat_coo2_down(counter_down)=j  
!                 counter_down=counter_down+1
!             end if                


!           end do
!          if(binary.eqv..FALSE.) close(217)
!          if(binary.eqv..TRUE.) close(260)
!          nonzero_up=counter_up-1
!          counter_up=nonzero_up
!          nonzero_down=counter_down-1
!          counter_down=nonzero_down


!        allocate(mat_val_csr_up(nonzero_up),
!     .            mat_csr1_up(nx+1)
!     .           ,mat_csr2_up(nonzero_up))
!        allocate(mat_val_csr_down(nonzero_down),mat_csr1_down(nx+1)
!     .           ,mat_csr2_down(nonzero_down))
      
!          write(*,*) myid,nonzero_up,nonzero_down
!          call MPI_REDUCE(counter_up,sum1,1,MPI_INTEGER,MPI_SUM,0,
!     .                     MPI_COMM_WORLD,ierr)
!          call MPI_REDUCE(counter_down,sum2,1,MPI_INTEGER,MPI_SUM,0,
!     .                     MPI_COMM_WORLD,ierr)
!          if(myid.eq.0) write(*,*) "Sum:",sum1+sum2
          
!       end if 
 

      call CPU_TIME(start)



!=============== MATRIX CONSTRUCTION IF SMALL ENOUGH ===================================================
!=============  SWITCH: IF DIM <= 20000 --> Save the LR-matrix, only 1 construction of each row


!        ROWS of LR-matrix per PE
       
!        if(Morb.eq.1) then       

          if(myid.lt.nprocs-1) then 
            rows=nloc
          else
            rows=nloc+mod(nx,nprocs)
          end if    

          allocate(v(rows,ncv))

!        else
 
!          if(myid.lt.nprocs-1) then 
!            rows_up=nloc_up
!          else
!            rows_up=nloc_up+mod(dimL_orb,nprocs)
!          end if    


!          if(myid.lt.nprocs-1) then 
!            rows_down=nloc_down
!          else
!            rows_down=nloc_down+mod(2*NC,nprocs)
!          end if    

!        end if

!       write(*,*) myid, rows
!       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!       stop



!=====  Calculate rows of matrix for each MPI process, last takes the rest
c         write(*,*) "MyID:", myid,"Rows:",rows

!         allocate(matrix(rows,nx+1))
!          do i=1, rows
!            matrix(i,:)=LR_mat(myid*nloc+i,:)
!          end do 


!        allocate(mat_csr1(maxnonzero),mat_val2(maxnonzero)
!     .           ,mat_csr2(rows+1))
        allocate(vec_new(nx))
        allocate(resf(ncv),vec2(nx),vec2_new(nx))
        allocate(final_evctrs(nx,ncv),evals(ncv)) 
        allocate(evals_sel(ncv))


!      Convert to sparse csr format

!      if(Morb.eq.1) then

!====  For storing matrix in sparse format
        mjob(1)=1
        mjob(2)=1
        mjob(3)=1
        mjob(4)=2
        mjob(5)=nonzero
        mjob(6)=0
      
        call mkl_zcsrcoo(mjob, rows, mat_val_csr, mat_csr2, mat_csr1,
     .             nonzero, mat_val, 
     .             mat_coo1, mat_coo2, info) 


        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!      else

!====  For storing matrix in sparse format
!        mjob(1)=1
!        mjob(2)=1
!        mjob(3)=1
!        mjob(4)=2
!        mjob(5)=nonzero_up
!        mjob(6)=0

!        call mkl_zcsrcoo(mjob, rows_up, mat_val_csr_up, mat_csr2_up, 
!     .             mat_csr1_up, nonzero_up, mat_val_up(1:nonzero_up), 
!     .             mat_coo1_up(1:nonzero_up), mat_coo2_up(1:nonzero_up)
!     .             , info) 


!        mjob(5)=nonzero_down

!        call mkl_zcsrcoo(mjob, rows_down, mat_val_csr_down,  
!     .             mat_csr2_down, mat_csr1_down, nonzero_down,  
!     .             mat_val_down(1:nonzero_down), 
!     .             mat_coo1_down(1:nonzero_down),
!     .             mat_coo2_down(1:nonzero_down), info) 


!      end if




c
c     %-------------------------------------------%
c     | M A I N   L O O P (Reverse communication) | 
c     %-------------------------------------------%
c

      conv_count=0  
     

      iter=0
      iter2=0 
      info=0


 10   continue
c
c        %---------------------------------------------%
c        | Repeatedly call the routine ZNAUPD and take |
c        | actions indicated by parameter IDO until    |
c        | either convergence is indicated or maxitr   |
c        | has been exceeded.                          |
c        %---------------------------------------------%
c




         call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
         call pznaupd ( MPI_COMM_WORLD, ido, bmat, rows, which, 
!         call pznaupd ( MPI_COMM_WORLD, ido, bmat, nx, which, 
!     &        nev, tol2, resid, ncv, v, ldv, iparam, ipntr, 
     &        nev, tol, resid, ncv, v, rows, iparam, ipntr, 
     &        workd, workl, lworkl, rwork,info )

!        if(myid.eq.0) then 
!         call znaupd ( ido, bmat, nx, which, 
!     .        nev, tol, resid, ncv, v, ldv, iparam, ipntr, 
!     .        workd, workl, lworkl, rwork,info )
!        end if
         call MPI_BARRIER(MPI_COMM_WORLD,ierr) 

!          do i=1, nprocs
!           if(myid.eq.i-1) write(*,*) rows,ipntr(1),ipntr(2)
!           if(myid.eq.i-1) write(*,*) rows,ipntr(8),workd(ipntr(8))
!           call MPI_BARRIER(MPI_COMM_WORLD,ierr) 
!          end do
!          stop 
    
!          call MPI_BCAST(ipntr,14,
!     .          MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

!          call MPI_BCAST(workd(ipntr(1)),nx,
!     .          MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)

!          call MPI_BCAST(ido,1,
!     .          MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

!         do i=1, ncv    
!          write(*,*) workd(ipntr(1)+i-1)
!         end do

         if (ido .eq. -1 .or. ido .eq. 1) then
c
c
c           %-------------------------------------------%
c           | Perform matrix vector multiplication      |
c           |                y <--- OP*x                |
c           | The user should supply his/her own        |
c           | matrix vector multiplication routine here |
c           | that takes workd(ipntr(1)) as the input   |
c           | vector, and return the matrix vector      |
c           | product to workd(ipntr(2)).               | 
c           %-------------------------------------------%
c

          iter2=iter2+1
 
c=========  OWN CONVERGENCE TEST
          done=0 
          if(myid.eq.0.and.mod(iter2,ncv).eq.0.and.iter.ge.0) then
             conv_count=0
             do i=1, ncv
c               if( (abs(workl(ipntr(8)+i-1))
c     .           -tol*abs(workl(ipntr(6)+i-1))).lt.0.0 ) then
               if( 
cARPACK crit     &             abs(workl(ipntr(8)+i-1)).lt.tol2*
cARPACK crit     &               abs(workl(ipntr(6)+i-1))
     &             abs(workl(ipntr(8)+i-1)).lt.tol2
     &        .and.abs(Real(workl(ipntr(6)+i-1))).lt.upperlimit 
     &        .and.Real(workl(ipntr(6)+i-1)).gt.0d0 
     &           ) conv_count=conv_count+1
!               end if
             end do 
!             if(conv_count.ge.nev.or.conv_count.ge.min_conv) done=1
             if(conv_count.ge.min_conv) done=1
          end if
          call MPI_BCAST(done,1,
     .          MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
          if(done.eq.1) goto 9000


       call CPU_TIME(t1)

c====   Combine all parts of v-vector from PEs
         vec_new=0.0

         do i=1, nprocs
          if(i.lt.nprocs) then
            if(myid.eq.i-1) then
              do k=1, nloc
               vec_new((i-1)*nloc+k)=workd(ipntr(1)+k-1) 
              end do 
            end if
            call MPI_BCAST(vec_new((i-1)*nloc+1:i*nloc),nloc,
     .                MPI_DOUBLE_COMPLEX,i-1,MPI_COMM_WORLD,ierr)
          else
            if(myid.eq.i-1) then
              do k=1, nloc+mod(nx,nprocs)
                vec_new((i-1)*nloc+k)=
     .                        workd(ipntr(1)+k-1) 
              end do
            end if
            call MPI_BCAST(vec_new((i-1)*nloc+1:i*nloc+mod(nx,nprocs)),
     .                nloc+mod(nx,nprocs),
     .                MPI_DOUBLE_COMPLEX,i-1,MPI_COMM_WORLD,ierr)
          end if
         end do 
!        write(*,*) "Myid:", myid, "Vec_new done"
!        write(*,*) nx,nprocs,nloc
!        call MPI_BARRIER(MPI_COMM_WORLD,ierr)


c====   Matrix*vector multiplication using sparse mkl routine


!      if(Morb.eq.1) then

!        call mkl_zcoogemv('N', rows, mat_val, mat_coo1,
!     .                   mat_coo2, nonzero, workd(ipntr(1)),
!     .                   workd(ipntr(2)+myid*nloc))


!         do i=1, nx    
!          vec_new(i)=workd(ipntr(1)+i-1)
!          write(*,*) vec_new(i)
!         end do

         call mkl_zcsrgemv('N', rows, mat_val_csr, mat_csr1,
     .                   mat_csr2, vec_new,
     .                   workd(ipntr(2)))


!      else

!        call mkl_zcoogemv('N', rows_up, mat_val_up, mat_coo1_up,
!     .                   mat_coo2_up, nonzero_up, workd(ipntr(1)),
!     .                   workd(ipntr(2)+myid*nloc_up))

!        call mkl_zcsrgemv('N', rows_up, mat_val_csr_up, mat_csr1_up,
!     .                   mat_csr2_up, workd(ipntr(1)),
!     .                   workd(ipntr(2)+myid*nloc_up))

!        call mkl_zcoogemv('N', rows_down, mat_val_down, mat_coo1_down,
!     .                   mat_coo2_down, nonzero_down, workd(ipntr(1)),
!     .                   workd(ipntr(2)+dimL_orb+myid*nloc_down))

!      call mkl_zcsrgemv('N', rows_down, mat_val_csr_down, mat_csr1_down,
!     .                   mat_csr2_down, workd(ipntr(1)),
!     .                   workd(ipntr(2)+dimL_orb+myid*nloc_down))

!      end if 

      
c====   squared or not
        squared=.FALSE.
        IF(squared.eqv..TRUE.) then
         vec_new=0d0 

         do i=1, nprocs
          if(i.lt.nprocs) then
            if(myid.eq.i-1) then
              vec_new((i-1)*nloc+1:i*nloc)=workd(ipntr(2):
     .                                       ipntr(2)+nloc) 
            end if
            call MPI_BCAST(vec_new((i-1)*nloc+1:i*nloc),nloc,
     .                MPI_DOUBLE_COMPLEX,i-1,MPI_COMM_WORLD,ierr)
          else
            if(myid.eq.i-1) then
              vec_new((i-1)*nloc+1:i*nloc+mod(nx,nprocs))=
     .                     workd(ipntr(2):ipntr(2)+nloc+mod(nx,nprocs)) 
            end if
            call MPI_BCAST(vec_new((i-1)*nloc+1:i*nloc+mod(nx,nprocs)),
     .                nloc+mod(nx,nprocs),
     .                MPI_DOUBLE_COMPLEX,i-1,MPI_COMM_WORLD,ierr)
          end if
         end do 

         call mkl_zcsrgemv('N', rows, mat_val_csr, mat_csr1,
     .                   mat_csr2, vec_new,
     .                   workd(ipntr(2)))

        END IF   

       call CPU_TIME(t2)
!       write(*,*) "Mult-Time:", t2-t1,"seconds"


       
       call CPU_TIME(t1)
     
!       write(*,*) "BCAST Time:", t2-t1,"seconds"
       


c====    Terminal output every NCV iterations from root process 0
        if(myid.eq.0.and.ishfts.eq.1.and.mod(iter2,ncv).eq.0) then
             iter=iter+1
          write(*,*)  "# root",
     .               "Ritz Values (Real)"
     .            ,"Ritz Values (Imag)"
     .            ,"Error"
c     .            ,"Converged for neg value"
          do i=1, ncv
            write(*,'(I6,10x,6(E20.10E3,10x))')
     .      ncv-i+1, 
     .      Real(workl(ipntr(6)+i-1)),
     .      DImag(workl(ipntr(6)+i-1)),
     .      abs(workl(ipntr(8)+i-1))
c     .      abs(workl(ipntr(8)+i-1))-tol2*abs(workl(ipntr(6)+i-1))
          end do

          call CPU_TIME(time)
     
          write(*,*)
          write(*,*) "Arnoldi iteration:", iter
          write(*,*) "ZNAUPD info:", info   
          write(*,*) "Converged values:", conv_count
          write(*,*) "ido:", ido
          write(*,*) "Time elapsed:", time-start,"seconds"
          write(*,*) "======================================="
 
!       Write current residual vector to file
!        open(unit=214, file="residual_vec.dat", 
!     .      status='replace', action='readwrite')
!          do j=1, dimL
!              write(214,'(2F20.10)') Real(resid(j))
!     .         ,Dimag(resid(j))
!          end do
!        close(214)


!       Write current arnoldi basis to file
!        open(unit=218, file="arnoldi_basis.dat", 
!     .      status='replace', action='readwrite')
!        do i=1, ncv
!          do j=1, dimL
!              write(218,'(2I10,10x,2F20.10)') i,j,Real(v(j,i))
!     .         ,Dimag(v(j,i))
!          end do
!        end do
!        close(218)

   
        end if
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)


c
c           %-----------------------------------------%
c           | L O O P   B A C K to call ZNAUPD again. |
c           %-----------------------------------------%
c
            go to 10
         end if



 9000 continue
c
c           %-----------------------------------------%
c           | Convergence reached. Continue with Residuals and
c           | Eigenvectors
c           %-----------------------------------------%
c

c=========  OWN CONVERGENCE TEST
          if(myid.eq.0.and.mod(iter2,ncv).eq.0.and.iter.ge.0) then
             conv_count=0
             done=0 
             do i=1, ncv
c               if( (abs(workl(ipntr(8)+i-1))
c     .           -tol*abs(workl(ipntr(6)+i-1))).lt.0.0 ) then
               if( 
     &             abs(workl(ipntr(8)+i-1)).le.tol2
!     &               abs(workl(ipntr(6)+i-1))
     &        .and.abs(Real(workl(ipntr(6)+i-1))).lt.upperlimit 
     &        .and.Real(workl(ipntr(6)+i-1)).gt.0d0 
     &           ) conv_count=conv_count+1
!               end if
             end do 
!             if(conv_count.ge.nev.or.conv_count.ge.min_conv) done=1
             if(conv_count.ge.min_conv) done=1
          end if
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

c====    Last normal Terminal output before final output and eigenvector calc
        if(myid.eq.0.and.ishfts.eq.1.and.mod(iter2,ncv).eq.0) then
             iter=iter+1
 
          write(*,*) "# root","             ",
     .               "Ritz Values","                  "
     .            ,"Error Ritz Values","              "
     .            ,"Converged for neg value"
          do i=1, ncv
            write(*,'(I6,10x,6(E20.10E3,10x))')
     .      ncv-i+1, 
     .      Real(workl(ipntr(6)+i-1)),
     .      Dimag(workl(ipntr(6)+i-1)),
!     .      abs(workl(ipntr(8)+i-1))-tol2*abs(workl(ipntr(6)+i-1))
     .      abs(workl(ipntr(8)+i-1))
          end do

          call CPU_TIME(time)
     
          write(*,*)
          write(*,*) "Arnoldi iteration:", iter
          write(*,*) "ZNAUPD info:", info   
          write(*,*) "Converged values:", conv_count
          write(*,*) "ido:", ido
          write(*,*) "Time elapsed:", time-start,"seconds"
          write(*,*) "======================================="

!       Write current residual vector to file
!        open(unit=214, file="residual_vec.dat", 
!     .      status='replace', action='readwrite')
!          do j=1, dimL
!              write(214,'(2F20.10)') Real(resid(j))
!     .         ,Dimag(resid(j))
!          end do
!        close(214)

!       Write current arnoldi basis to file
!        open(unit=218, file="arnoldi_basis.dat", 
!     .      status='replace', action='readwrite')
!        do i=1, ncv
!          do j=1, dimL
!              write(218,'(2I10,10x,2F20.10)') i,j,Real(v(j,i))
!     .         ,Dimag(v(j,i))
!          end do
!        end do
!        close(218)
   
        end if
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)


!      PARPACK convergence criterion
c          nconv  = 0
c          do 25 i = 1, ncv
c             rtemp = max( eps23, dlapy2 ( dble (ritz(np+i)),
c     &                                  dimag (ritz(np+i)) ) )
c             if ( dlapy2 (dble (bounds(np+i)),dimag (bounds(np+i)))
c     &                 .le. tol*rtemp ) then
c                nconv = nconv + 1
c             end if
c      25    continue


      
        goto 18
 
         if(myid.eq.0) iparam(5)=conv_count 
          call MPI_BCAST(iparam(5),1,
     .          MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         rvec=.TRUE.
         howmny='A' 
       allocate(d(nev+1))
       allocate(z(nx,nev))
       z=0d0
       ldz=nx

!===   Create array of selected eigenvalues
       allocate(sel(ncv))
       if(myid.eq.0) write(*,*) "HERE"
       do i=1,ncv
         if( abs(workl(ipntr(8)+i-1))
     .           .lt.tol2 
     .        .and.abs(Real(workl(ipntr(6)+i-1))).lt.upperlimit 
     .        .and.Real(workl(ipntr(6)+i-1)).gt.0d0
     .     ) then    
                   sel(i)=.TRUE.
         else
                   sel(i)=.FALSE.
         end if
       end do 

       if(myid.eq.0) then
         call pzneupd
     &         ( MPI_COMM_WORLD , rvec  , howmny, sel, d    ,
     &           z    , ldz   , sigma , workev, bmat ,
     &           rows    , which , nev   , tol2   , resid,
     &           ncv  , v     , ldv   , iparam, ipntr,
     &           workd, workl , lworkl, rwork , info )
       end if  

         allocate(vec_hilf(rows))
         call mkl_zcsrgemv('N', rows, mat_val_csr, mat_csr1,
     .                   mat_csr2, v(:,1),
     .                   vec_hilf)

           vec_new=0d0
         do i=1, nprocs
          if(i.lt.nprocs) then
            if(myid.eq.i-1) then
              vec_new((i-1)*nloc+1:i*nloc)=vec_hilf
            end if
            call MPI_BCAST(vec_new((i-1)*nloc+1:i*nloc),nloc,
     .                MPI_DOUBLE_COMPLEX,i-1,MPI_COMM_WORLD,ierr)
          else
              vec_new((i-1)*nloc+1:i*nloc+mod(nx,nprocs))=vec_hilf
            call MPI_BCAST(vec_new((i-1)*nloc+1:i*nloc+mod(nx,nprocs)),
     .                nloc+mod(nx,nprocs),
     .                MPI_DOUBLE_COMPLEX,i-1,MPI_COMM_WORLD,ierr)
          end if
         end do 
       if(myid.eq.0) then
         do i=1, nev+1
           write(*,*) d(i)
         end do
           write(*,*) sqrt(dot_product(vec_new,vec_new))
       end if 
       if(myid.eq.0) write(*,*) "LINE 1124", info
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       stop 


18        continue
c=========   OWN CODE TOP CALCULATE EIGENVECTORS AND RESIDUALS OF
c=========        WANTED!!!     RITZ VALUES


c==== Store final arnoldi basis on root process
          if(myid.eq.0) allocate(v_final(nx,ncv))

          if(.not.allocated(v_inter)) 
     .          allocate(v_inter(nloc+mod(nx,nprocs),ncv)) 

          if(myid.eq.0) then
            do i=1, ncv
              v_final(1:nloc,i)=v(:,i)
            end do
          end if

!          DO k=1, nprocs
!             if(myid.eq.(k-1)) write(*,*) "myid:", myid, v(1,1)
!          END DO
!          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!          stop 

          DO k=1, nprocs-1
           if(myid.eq.0) v_inter=0d0
           IF(k.ne.nprocs-1) then     
            if(myid.eq.0) then
               do i=1, ncv
                 CALL MPI_RECV(v_inter(1:nloc,i),nloc,
     .               MPI_DOUBLE_COMPLEX,k,2,
     .               MPI_COMM_WORLD,s_status,ierr) 
                 do j=1, nloc
                  v_final(k*nloc+j,i)=v_inter(j,i)
                 end do
               end do 
            end if  
            if(myid.eq.k) then
                do i=1, ncv
                 do j=1, nloc
                  v_inter(j,i)=v(j,i)
                 end do  
                 CALL MPI_SEND(v_inter(1:nloc,i),nloc,
     .               MPI_DOUBLE_COMPLEX,0,2,
     .               MPI_COMM_WORLD,s_status,ierr) 
                end do
            end if
           ELSE
            if(myid.eq.0) then
               do i=1, ncv
                 CALL MPI_RECV(v_inter(1:nloc+mod(nx,nprocs),i),
     .               nloc+mod(nx,nprocs),MPI_DOUBLE_COMPLEX,k,1,
     .               MPI_COMM_WORLD,s_status,ierr) 
                 do j=1, nloc+mod(nx,nprocs)
                  v_final(k*nloc+j,i)=
     .                     v_inter(j,i)
                 end do          
               end do 
            end if  
            if(myid.eq.k) then
                do i=1, ncv
                 do j=1, nloc+mod(nx,nprocs)
                  v_inter(j,i)=
     .                        v(j,i)
                 end do 
                 CALL MPI_SEND(v_inter(1:nloc+mod(nx,nprocs),i),
     .               nloc+mod(nx,nprocs),MPI_DOUBLE_COMPLEX,0,1,
     .               MPI_COMM_WORLD,s_status,ierr) 
                end do
            end if
           END IF  
          END DO    

          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
          if(allocated(v_inter)) 
     .          deallocate(v_inter) 




      if(myid.eq.0) then


      write(*,*) "%-----------------------------------------------%"
      write(*,*) conv_count,"Ritz values converged"
      write(*,*) "starting now eigenvector/residual calculation..."
      write(*,*) "%-----------------------------------------------%"

!===   Create array of selected eigenvalues
       allocate(sel(ncv))
       do i=1,ncv
         if( 
c     &             abs(workl(ipntr(8)+i-1)).lt.tol2*
c     &               abs(workl(ipntr(6)+i-1))
     .    abs(workl(ipntr(8)+i-1))
     .           .lt.tol2 
     .        .and.abs(Real(workl(ipntr(6)+i-1))).lt.upperlimit 
     .        .and.Real(workl(ipntr(6)+i-1)).gt.0d0
     .     ) then    
                   sel(i)=.TRUE.
         else
                   sel(i)=.FALSE.
         end if
       end do 


!===    allocate arrays for eigenvector calc
        allocate(hess(ncv,ncv),vl(ncv,ncv),vr(ncv,ncv),
     .             work_hess(ncv*ncv),rwork_hess(ncv),
     .             ifaill(ncv),ifailr(ncv),hess2(ncv,ncv))
   
!===      Define hess and evals
         do i=1, ncv
           evals(i)=workl(ipntr(6)+i-1) 
           do j=1, ncv
             hess(j,i)=workl(ipntr(5)+j-1+(i-1)*ncv)
           end do
         end do
         hess2=hess





        

!=====    Eigenvectors of Hessenberg matrix

!=====   sort eigenvalues
         allocate(eigenw_temp(ncv),sel_temp(ncv))
         eigenw_temp=evals
         sel_temp=sel
         do i = 1, ncv
                do j = 2, ncv
                 if(Real(eigenw_temp(j-1)).gt.Real(eigenw_temp(j))) then
                      tempp_log=sel_temp(j-1)
                      tempp = eigenw_temp(j-1);
                      eigenw_temp(j-1) = eigenw_temp(j);
                      eigenw_temp(j) = tempp;
                      sel_temp(j-1) = sel_temp(j);
                      sel_temp(j) = tempp_log;
                 end if
                end do
         end do
         evals=eigenw_temp

        
         k=1
         do i=1, ncv
           if(sel_temp(i)) then 
            evals_sel(k)=eigenw_temp(i) 
            k=k+1
           end if
         end do
         call zhsein('R','N','N',sel_temp,ncv,hess,ncv,eigenw_temp,
     .               vl,ncv,vr,ncv,ncv,conv_count,work_hess,
     .               rwork_hess,ifaill,ifailr,infoo)
c         call zhsein('R','Q','N',select,ncv,hess,ncv,evals,
c     .               vl,ncv,vr,ncv,ncv,ncv,work_hess,
c     .               rwork_hess,ifaill,ifailr,infoo)

!         write(*,*) "Infoo of zhsein:",infoo
!         if(infoo.ne.0) write(*,*) ifailr


!=====   Normalize eigenvectors of Hessenberg matrix
          do i=1, conv_count
            vr(:,i)=vr(:,i)/Sqrt(dot_product(vr(:,i),vr(:,i)))
          end do


!====   Transform eigenvectors to full space

          call zgemm('N','N',nx,ncv,ncv,alpha,v_final,nx,vr,ncv,
     .               beta,final_evctrs,nx)
           
!=====   Normalize eigenvectors of Hessenberg matrix
          do i=1, conv_count
            final_evctrs(:,i)=final_evctrs(:,i)
     .     /Sqrt(dot_product(final_evctrs(:,i),final_evctrs(:,i)))
          end do


   
!====  Output Hessenberg matrix
c           write(file2,'(A10,I1,A4)') 'hessenberg',myid,'.dat'
c       open(unit=14, file=file2,
c     . status='replace',action='READWrite',
c     .iostat=status_open)

c           do i=1, ncv !nx
c             write(14,'(100(2F20.10,10x))') 
c     .          (hess(i,j), 
c     .            j=1,ncv,1)
c           end do 
c       close(14)




      end if
!      end of proc 0 calculation so far..       
       

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)


!=====   Calculation of Residuals on all PEs
        call MPI_BCAST(evals_sel,ncv,
     .          MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_BCAST(final_evctrs,nx*ncv,
     .          MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        call MPI_BCAST(conv_count,1,
     .          MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  



        do j=1,  conv_count

             vec2(:)=final_evctrs(:,j)

!         call mkl_zcsrgemv('N', rows, mat_val, mat_csr2,
!     .                   mat_csr1, vec2,
!     .                   vec2_new(1+myid*nloc))



!      if(Morb.eq.1) then

!         call mkl_zcoogemv('N', rows, mat_val, mat_coo1,
!     .                   mat_coo2, nonzero, vec2,
!     .                   vec2_new(1+myid*nloc))

         vec2_new=0d0
         call mkl_zcsrgemv('N', rows, mat_val_csr, mat_csr1,
     .                   mat_csr2, vec2,
     .                   vec2_new(1+myid*nloc))



          do i=1, nprocs
            if(i.lt.nprocs) then
              call MPI_BARRIER(MPI_COMM_WORLD,ierr)
              call MPI_BCAST(vec2_new(1+(i-1)*nloc),nloc,
     .          MPI_DOUBLE_COMPLEX,i-1,MPI_COMM_WORLD,ierr)
              call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            else
              call MPI_BARRIER(MPI_COMM_WORLD,ierr)
              call MPI_BCAST(vec2_new(1+(i-1)*nloc),nloc+mod(nx,nprocs),
     .          MPI_DOUBLE_COMPLEX,i-1,MPI_COMM_WORLD,ierr)
              call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            end if
          end do

!      else

!        call mkl_zcoogemv('N', rows_up, mat_val_up, mat_coo1_up,
!     .                   mat_coo2_up, nonzero_up, vec2,
!     .                   vec2_new(1+myid*nloc_up))


!        call mkl_zcoogemv('N', rows_down, mat_val_down, mat_coo1_down,
!     .                   mat_coo2_down, nonzero_down, vec2,
!     .                   vec2_new(dimL_orb+1+myid*nloc_down))


!          do i=1, nprocs
!            if(i.lt.nprocs) then
!              call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!              call MPI_BCAST(vec2_new(1+(i-1)*nloc_up),nloc_up,
!     .          MPI_DOUBLE_COMPLEX,i-1,MPI_COMM_WORLD,ierr)
!              call MPI_BCAST(vec2_new(dimL_orb+1+(i-1)*nloc_down),
!     .          nloc_down,MPI_DOUBLE_COMPLEX,i-1,MPI_COMM_WORLD,ierr)
!              call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!            else
!              call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!              call MPI_BCAST(vec2_new(1+(i-1)*nloc_up),
!     .          nloc_up+mod(dimL_orb,nprocs),
!     .          MPI_DOUBLE_COMPLEX,i-1,MPI_COMM_WORLD,ierr)
!              call MPI_BCAST(vec2_new(dimL_orb+1+(i-1)*nloc_down),
!     .          nloc_down+mod(2*NC,nprocs),
!     .          MPI_DOUBLE_COMPLEX,i-1,MPI_COMM_WORLD,ierr)
!              call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!            end if
!          end do


!      end if 


c         call zgemv('N',rows,nx,alpha,matrix,rows,
c     .         vec2,1,beta,
c     .         vec2_new(1+myid*nloc),1)




         if(myid.eq.0) then
          vec2=0.0
          vec2(:)=final_evctrs(:,j)
          call zaxpy(nx, -evals_sel(j), vec2, 1,vec2_new,1)
          resf(j)=sqrt(dot_product(vec2_new,vec2_new))
c          call cblas_zdotc_sub(nx,vec2_new,1,vec2_new,1,resf(j))
c          write(*,*) "Res",j,Real(resf(j))
         end if 
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        end do
 
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      call MPI_FINALIZE()





!=====================================================================================================
!==========   OUTPUT eigenvectors + values / Type / Density Response / Response weights  =============
!=====================================================================================================

!====  OUTPUT and preliminary calculations only done by process 0


!==========
!====  Output eigenvectors (split up in orbital and CI parts)  
      if(myid.eq.0) then

           path='Eigenv'
           call mk_dir(path) 
           pathORB='Eigenv'//path_sep//'orb'//path_sep
           call mk_dir(pathORB)
           pathCI='Eigenv'//path_sep//'CI'//path_sep
           call mk_dir(pathCI)
           pathType='Eigenv'//path_sep//'typ'//path_sep
           call mk_dir(pathType)
           pathRes='Eigenv'//path_sep//'res'//path_sep
           call mk_dir(pathRes)
!           pathRij='Eigenv'//path_sep//'rij'//path_sep
!           call mk_dir(pathRij)

      DIM: select case (dim_mctdhb)

       case(1)

           do i=1, conv_count
            if(i.le.9) write(file2,'(A17,I1,A4)')
     .                  'eigenvectors_orb0',i,'.dat'
            if(i.gt.9.and.i.lt.100) write(file2,'(A16,I2,A4)') 
     .                  'eigenvectors_orb',i,'.dat'
            if(i.gt.99.and.i.lt.1000) write(file2,'(A16,I3,A4)') 
     .                  'eigenvectors_orb',i,'.dat'
            open(unit=123, file=trim(adjustl(pathORB))//file2,
     .       status='replace',action='READWrite',
     .       iostat=status_open)

          write(123,'(100(A15))') "x","    ","Re(u_i)","   ","Im(u_i)",
     .            "   ","Re(v_i)","   ","Im(v_i)","   ", "i=1..Morb"

            counter=1
            do k=1, NDX
             write(123,'(F20.10,10x,1000(F20.10,10x,F20.10,10x))') 
     .             Ort_x(k),
     .            (Real(final_evctrs((ll-1)*ND+counter,i)), 
     .             Dimag(final_evctrs((ll-1)*ND+counter,i)),
     .             Real(final_evctrs((ll-1+Morb)*ND+counter,i)), 
     .             Dimag(final_evctrs((ll-1+Morb)*ND+counter,i)),
     .             ll=1,Morb,1)
               counter=counter+1
            end do
            close(123)
!              write(123,*) 
           end do



          do j=1, conv_count
           if(j.le.9) write(file2,'(A16,I1,A4)') 
     .            'eigenvectors_CI0',j,'.dat'
           if(j.gt.9.and.j.lt.100) write(file2,'(A15,I2,A4)') 
     .            'eigenvectors_CI',j,'.dat'
           if(j.gt.99.and.j.lt.1000) write(file2,'(A15,I3,A4)') 
     .                  'eigenvectors_CI',j,'.dat'
           open(unit=133, file=trim(adjustl(pathCI))//file2,
     .       status='replace',action='READWrite',
     .       iostat=status_open)


          write(133,'(100(A15))') "Real(C_u)","   ","Imag(C_u)","   ",
     .                 "Real(C_v)","   ","Imag(C_v)"

           do i=1, Nconf
               write(133,'(1000(F20.10,10x,F20.10,10x))')
     .             Real(final_evctrs(dimL_orb+i,j)), 
     .             Dimag(final_evctrs(dimL_orb+i,j)),
     .             Real(final_evctrs(dimL_orb+NC+i,j)), 
     .             Dimag(final_evctrs(dimL_orb+NC+i,j))

           end do
           close(133)
          end do
 


       case(2)


           do i=1, conv_count
            if(i.le.9) write(file2,'(A17,I1,A4)')
     .                  'eigenvectors_orb0',i,'.dat'
            if(i.gt.9) write(file2,'(A16,I2,A4)') 
     .                  'eigenvectors_orb',i,'.dat'
            open(unit=123, file=trim(adjustl(pathORB))//file2,
     .       status='replace',action='READWrite',
     .       iostat=status_open)

          write(123,'(100(A15))') "x","    ","y","    ","Real(u_i)",
     .                 "   ","Imag(u_i)","   ",
     .                 "Real(v_i)","   ","Imag(v_i)","   ", "i=1..Morb"

           counter=1
            do j=1, NDY
             do k=1, NDX  
               write(123,'(F20.10,10x,F20.10,10x,
     .                 500(E20.10,10x,E20.10,10x))') 
     .             Ort_x(k),Ort_y(j),
     .            (Real(final_evctrs((ll-1)*ND+counter,i)), 
     .             Dimag(final_evctrs((ll-1)*ND+counter,i)),
     .             Real(final_evctrs((ll-1+Morb)*ND+counter,i)), 
     .             Dimag(final_evctrs((ll-1+Morb)*ND+counter,i)),
     .             ll=1,Morb,1)
               counter=counter+1
             end do
             write(123,*) 
            end do
            close(123) 
           end do



          do j=1, conv_count
           if(j.le.9) write(file2,'(A16,I1,A4)') 
     .            'eigenvectors_CI0',j,'.dat'
           if(j.gt.9) write(file2,'(A15,I2,A4)') 
     .            'eigenvectors_CI',j,'.dat'
           open(unit=133, file=trim(adjustl(pathCI))//file2,
     .       status='replace',action='READWrite',
     .       iostat=status_open)


          write(133,'(100(A15))') "Real(C_u)","   ","Imag(C_u)","   ",
     .                 "Real(C_v)","   ","Imag(C_v)"

           do i=1, Nconf
               write(133,'(1000(F20.10,10x,F20.10,10x))')
     .             Real(final_evctrs(dimL_orb+i,j)), 
     .             Dimag(final_evctrs(dimL_orb+i,j)),
     .             Real(final_evctrs(dimL_orb+NC+i,j)), 
     .             Dimag(final_evctrs(dimL_orb+NC+i,j))

           end do
           close(133) 
          end do

 
       end select DIM




!=====================
!====== Type CI or Orb

        if(Morb.gt.1) then  
          allocate(cu_vec(NC),cv_vec(NC))       
          allocate(u_vec(Morb,ND),v_vec(Morb,ND))  
          allocate(t(conv_count))       
          write(file2,'(A8)') "type.dat"
            open(unit=133, file=trim(adjustl(pathType))//"type.dat",
     .       status='replace',action='READWrite',
     .       iostat=status_open)

          write(133,'(100(A15))') "#","    ","type","   ","sum_u","   ",
     .                 "sum_v","   ","Cu cont.","   ", "Cv cont.",
     .                 "   ","u_i(1...Morb)","   ","v_i(1...Morb)" 
          do j=1, conv_count
            cu_vec=final_evctrs(dimL_orb+1:dimL_orb+NC,j)
            cv_vec=final_evctrs(dimL_orb+NC+1:dimL_orb+2*NC,j)
            sum_u=0d0;sum_v=0d0; 
            do i=1, Morb
              u_vec(i,:)=final_evctrs((i-1)*ND+1:i*ND,j)
              v_vec(i,:)=final_evctrs(dimL_orb/2+(i-1)*ND+1:
     .                                 dimL_orb/2+i*ND,j)
              sum_u=sum_u+Real(dot_product(u_vec(i,:),u_vec(i,:)))
              sum_v=sum_v+Real(dot_product(v_vec(i,:),v_vec(i,:)))
            end do
            t(j)=dot_product(cu_vec,cu_vec)-
     .           dot_product(cv_vec,cv_vec)

            write(133,'(I4,50F20.10)') j,Real(t(j)),sum_u,sum_v,
     .             Real(dot_product(cu_vec,cu_vec)),
     .             Real(dot_product(cv_vec,cv_vec)),
     .            (Real(dot_product(u_vec(k,:),u_vec(k,:))),k=1,Morb,1),
     .            (Real(dot_product(v_vec(k,:),v_vec(k,:))),k=1,Morb,1)
          end do

          close(133)

        end if



!=======================
!===== DENSITY RESPONSES


!      READ ZRIJ
!       if(myid.eq.0) then  
          open(unit=218, file='ZRIJ.out', 
     .      status='old', action='READ', iostat=status_open)
            do   
              read(218,*,iostat=st) i,a,b 
              if(st.ne.0) exit
              ZRIJ(i)=Cmplx(a,b,kind=8)
            end do
          close(218)
 
        allocate(dens_resp_orb(conv_count,dimL_orb),
     .           Rho_InvSqr(Morb,Morb),Comp_Rho_InvSqr(Morb,Morb),
     .           RhoSQR1(Morb,Morb),RhoSQR2(Morb,Morb),Rho(Morb,Morb))
!        if(.not.allocated(PSI)) allocate( PSI(NDX*NDY*NDZ,Morb))


!     Check if ZRIJ is real
        is_real=.TRUE. 
        do i=1, Morb*Morb
           if(Dimag(ZRIJ(i)).gt.1e-8) is_real=.FALSE.
        end do
 
!     Calculate sqrt of ZRIJ
      Rho_InvSqr=0d0; Comp_Rho_InvSqr=0d0; Rho=0d0;
      call Get_Full_Rij(Rho)
      call squarerootRInv_arnoldi(Rho_InvSqr)
      Comp_Rho_InvSqr=Cmplx(Rho_InvSqr)  
      RhoSQR1=matmul(Rho,Rho_InvSqr)
      RhoSQR2=matmul(Comp_Rho_InvSqr,Rho)  


      dens_resp_orb=0d0


!     == ORBITAL PART ==
!     CASE OF REAL ZRIJ
!        if(is_real.eqv..TRUE.) then 
         do k=1, conv_count
           do i=1, Morb
            do j=1, Morb 
              do ll=1, NDX*NDY*NDZ
!         write(*,*) PSI
               

                 dens_resp_orb(k,ll)=dens_resp_orb(k,ll)+
     .            RhoSQR1(i,j)*
     .            Conjg(PSI(ll,i))*
     .            (final_evctrs((j-1)*ND+ll,k)+
     .            Conjg(final_evctrs((j-1+Morb)*ND+ll,k)))+
     .            (final_evctrs((i-1)*ND+ll,k)+
     .            Conjg(final_evctrs((i-1+Morb)*ND+ll,k)))*
     .            RhoSQR2(i,j)*PSI(ll,j)

              end do                  
            end do
           end do
          
         end do
!        end if


!    == CI PART ==
       allocate(cu(NC,conv_count),cv(NC,conv_count))  
       allocate(CCu(Morb,Morb,conv_count),CCv(Morb,Morb,conv_count))  
       allocate(dens_resp_CI(dimL_orb,conv_count))

       Cu=0d0;Cv=0d0;CCu=0d0;CCv=0d0;dens_resp_CI=0d0;  

       DO k=1, conv_count
          Cu(:,k)=final_evctrs(dimL_orb+1:dimL_orb+NC,k)
          Cv(:,k)=final_evctrs(dimL_orb+NC+1:dimL_orb+2*NC,k)

          DO i=1,Morb
             DO j=1,Morb

                VIN_help=VIN; Vhelp=0.0d0; 
                IF(Morb.gt.1) THEN
                   CALL Produce_Cij(VIN_help,Vhelp,i,j,0,0,Nconf)
                END IF
                CCv(i,j,k)=dot_product(Conjg(Cv(:,k)),Vhelp(:))

                VIN_help=VIN; Vhelp=0.0d0; 
                IF(Morb.gt.1) THEN
                   CALL Produce_Cij(VIN_help,Vhelp,j,i,0,0,Nconf)
                END IF
                CCu(j,i,k)=dot_product(Vhelp(:),Cu(:,k))
!                CCu(j,i,k)=SUM(Cu(:,k)*CONJG(Vhelp(:)))
!                CCv(i,j,k)=SUM(Cv(:,k)*Vhelp(:))

             END DO
          END DO
       END DO

!      Multiply with orbitals

       dens_resp_CI=0.0d0  
       DO k=1, conv_count
          DO i=1,Morb
             DO j=1,Morb

                dens_resp_CI(:,k)=dens_resp_CI(:,k)+  
     .            CONJG(PSI(:,i))*PSI(:,j)*(CCu(i,j,k)+CCv(i,j,k))

             END DO
          END DO
       END DO
!   ====


!     == OUTPUT ==
        do i=1, conv_count
            if(i.le.9) write(file2,'(A10,I1,A4)')
     .                  'densresp_0',i,'.dat'
            if(i.gt.9) write(file2,'(A9,I2,A4)') 
     .                  'densresp_',i,'.dat'
            if(i.gt.99.and.i.lt.1000) write(file2,'(A9,I3,A4)') 
     .                  'densresp_',i,'.dat'
            open(unit=144, file=trim(adjustl(pathRes))//file2,
     .       status='replace',action='READWrite',
     .       iostat=status_open)



          if(dim_mctdhb.eq.1) then


          write(144,'(100(A15))') "x","   ","Re(DR_orb)","   "
     .                 ,"Im(DR_orb)","   ",
     .                 "Re(DR_CI)","   ","Im(DR_CI)",
     .                 "   ", "Re(DR_tot)","    ",
     .                 "Im(DR_tot)"

             do k=1, NDX  
               write(144,'(F20.10,10x,
     .                 E20.10,10x,E20.10,10x,E20.10,10x
     .                 E20.10,10x,E20.10,10x,E20.10 )') 
     .             Ort_x(k),
     .             Real(dens_resp_orb(i,k)), 
     .             Dimag(dens_resp_orb(i,k)),
     .             Real(dens_resp_CI(k,i)), 
     .             Dimag(dens_resp_CI(k,i)),
     .             Real(dens_resp_orb(i,k))+Real(dens_resp_CI(k,i)), 
     .             Dimag(dens_resp_orb(i,k))+Dimag(dens_resp_CI(k,i))
             end do
             write(144,*) 
          end if 

          if(dim_mctdhb.eq.2) then


          write(144,'(100(A15))') "x","   ","y","   ","Re(DR_orb)","   "
     .                 ,"Im(DR_orb)","   ",
     .                 "Re(DR_CI)","   ","Im(DR_CI)",
     .                 "   ", "Re(DR_tot)","    ",
     .                 "Im(DR_tot)"

            do j=1, NDY
             do k=1, NDX  
               write(144,'(F20.10,10x,F20.10,10x,
     .                 E20.10,10x,E20.10,10x,E20.10,10x,E20.10,10x,  
     .                 E20.10,10x,E20.10)') 

     .             Ort_x(k),Ort_y(j),

     .             Real(dens_resp_orb(i,(j-1)*NDX+k)), 
     .             Dimag(dens_resp_orb(i,(j-1)*NDX+k)),
     .             Real(dens_resp_CI((j-1)*NDX+k,i)), 
     .             Dimag(dens_resp_CI((j-1)*NDX+k,i)),
     .             Real(dens_resp_orb(i,(j-1)*NDX+k))+
     .                   Real(dens_resp_CI((j-1)*NDX+k,i)), 
     .             Dimag(dens_resp_orb(i,(j-1)*NDX+k))+
     .                   Dimag(dens_resp_CI((j-1)*NDX+k,i))

             end do
             write(144,*) 
            end do
          end if 

            close(144) 
        end do

!=============


!      goto 8001 
!======================================================
!=====  RESPONSE  WEIGHTS  (taken from the old LR code)

         allocate(RhoSQR(Morb,Morb))
         allocate(tr(ND,2))
         allocate(weight_orb(conv_count,2),weight_CI(conv_count,2))
         allocate(ut(ND,conv_count,Morb),vt(ND,conv_count,Morb))
         allocate(u2(ND,conv_count,Morb),v2(ND,conv_count,Morb))
         allocate(CCu_r(Morb,Morb,conv_count))
         allocate(CCv_r(Morb,Morb,conv_count))
         allocate(orb_int(Morb,Morb,2))


!        ORBITAL!!!!!!!!!!!!!!!!!!!!!!

!        first get square roots of densities
         IF(Morb.gt.1) THEN
            CALL squarerootR_arnoldi(RhoSQR)
         ELSE 
            RhoSQR(1,1)=SQRT(REAL(AllZRIJ(1,1)))
         END IF

!        multiply square roots with amplitudes

          do j=1, conv_count
            do i=1, Morb
              u2(:,j,i)=final_evctrs((i-1)*ND+1:i*ND,j)
              v2(:,j,i)=final_evctrs(dimL_orb/2+(i-1)*ND+1:
     .                                 dimL_orb/2+i*ND,j)
            end do 
          end do

         ut=0.0d0; vt=0.0d0

         DO k=1,conv_count
            DO j=1,Morb
               DO i=1,Morb

                  ut(:,k,j)=ut(:,k,j)+RhoSQR(j,i)*Conjg(u2(:,k,i))
                  vt(:,k,j)=vt(:,k,j)+RhoSQR(j,i)*v2(:,k,i)
                
               END DO
            END DO
         END DO

!        weights for different external driving: f=x,x^2
             tr=(0d0,0d0)   !tr=f^+=f^-
          if(dim_mctdhb.eq.1) then 
            Do i=1,ND
            tr(i,1)=ort_x(i)
            tr(i,2)=ort_x(i)**2
            ENDDO
          end if
          if(dim_mctdhb.eq.2) then 
            ind=1  
            Do i=1,NDX
              Do j=1, NDY
                tr(ind,1)=sqrt(ort_x(i)**2+ort_y(j)**2)
                tr(ind,2)=ort_x(i)**2+ort_y(j)**2
                ind=ind+1  
              ENDDO
            ENDDO
          end if

         weight_orb=0.0d0  
         DO k=1,conv_count
            DO i=1,Morb

               weight_orb(k,1)=weight_orb(k,1) +  
     .   SUM(tr(:,1) *(PSI(:,i)*ut(:,k,i)+   
     .   CONJG(PSI(:,i))*CONJG(vt(:,k,i))))
              weight_orb(k,2)=weight_orb(k,2) +  
     .   SUM(tr(:,2)*(PSI(:,i)*ut(:,k,i)+  
     .   CONJG(PSI(:,i))*CONJG(vt(:,k,i))))
            END DO


         END DO

!      CI!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!      Produce Cij

       DO k=1,conv_count
          DO i=1,Morb
             DO j=1,Morb

                VIN_help=VIN; Vhelp=0.0d0; 
                IF(Morb.gt.1) THEN
                    CALL Produce_Cij(VIN_help,Vhelp,i,j,0,0,Nconf)
                END IF
                CCu_r(i,j,k)=SUM(CONJG(Cu(:,k))*Vhelp(:))
                CCv_r(j,i,k)=SUM(CONJG(Cv(:,k)*Vhelp(:)))
             END DO
          END DO
       END DO

!      Orbital integrals for different external driving: f=1,x,x^2

       DO i=1,Morb
          DO j=1,Morb

             orb_int(i,j,1)=SUM(CONJG(PSI(:,i))*PSI(:,j)*tr(:,1))
             orb_int(i,j,2)=SUM(CONJG(PSI(:,i))*PSI(:,j)*tr(:,2))

          END DO
       END DO

!      Weights for different external driving: f=1,x,x^2

         weight_CI=0.0d0  
         DO k=1,conv_count
            DO i=1,Morb
               DO j=1,Morb

                  weight_CI(k,1)=weight_CI(k,1) +  
     .               orb_int(i,j,1)*(CCu_r(i,j,k)+CCv_r(i,j,k))
                  weight_CI(k,2)=weight_CI(k,2) +  
     .               orb_int(i,j,2)*(CCu_r(i,j,k)+CCv_r(i,j,k))
!            IF(k==265) WRITE(196,*) ABS(wo(k)), i,j,  
!     .                            CCu_r(i,j,k),CCv_r(i,j,k)
               END DO
            END DO
         END DO

       path='Eigenv'
            open(unit=128, file=trim(adjustl(path))//'/'//
     .                 'response_weights.dat',
     .       status='replace',action='READWrite',
     .       iostat=status_open)

          write(128,'(100(A15))') 
     .      "# root","   ",
     .               "w_i (Re)","   ",
     .            "w_i (Im)","   ",
     .            "w_i (abs)","   ","Op=1,2"
          write(128,*) '  '

          do i=1, conv_count
            write(128,'(I6,10x,6(E20.10E3,10x))')
     .      i, 
     .      (Real(weight_orb(i,j)+weight_CI(i,j)),
     .      Dimag(weight_orb(i,j)+weight_CI(i,j)),
     .      Abs(weight_orb(i,j)+weight_CI(i,j)),j=1,2,1)
          end do

            close(128)

!8001  continue



!================================
!====  Output final Arnoldi basis
!           write(file2,'(A7,I1,A4)') 'a_basis',myid,'.dat'
!       open(unit=14, file=file2,
!     . status='replace',action='READWrite',
!     .iostat=status_open)

!           do i=1, nx
!             write(14,'(100(F20.10,10x,F20.10,10x))') 
!     .          (Real(v_final(i,j)), 
!     .             Dimag(v_final(i,j)),j=1,ncv,1)
!           end do 
!       close(14)





c=====   FINAL OUTPUT EIGENVALUES (SCREEN AND FILE)


      write(*,*) " "
      write(*,*) "Eigenvectors and residuals calculated."
      print *, '  '
      print *, '  '

       path='Eigenv'
            open(unit=123, file=trim(adjustl(path))//'/'//
     .                 'eigenvalues.dat',
     .       status='replace',action='READWrite',
     .       iostat=status_open)

          write(123,'(A6,10x,A18,10x,A18,10x,A9)') 
     .      "# root",
     .               "Ritz Values (Real)"
     .            ,"Ritz Values (Imag)"
     .            ,"Residuals"
          write(123,*) '  '
          do i=1, conv_count
            write(123,'(I6,10x,6(E20.10E3,10x))')
     .      i, 
     .      Real(evals_sel(i)),
     .      Dimag(evals_sel(i)),
c     .      Sqrt(abs(Real(evals_sel(i)))/2.0),
c     .      Sqrt(abs(Dimag(evals_sel(i)))/2.0),
     .      Real(resf(i))
          end do

           close(123)


          write(*,'(A6,10x,A18,10x,A18,10x,A9)') 
     .      "# root",
     .               "Ritz Values (Real)"
     .            ,"Ritz Values (Imag)"
     .            ,"Residuals"
      print *, '  '
          do i=1, conv_count
            write(*,'(I6,10x,6(E20.10E3,10x))')
     .      i, 
     .      Real(evals_sel(i)),
     .      Dimag(evals_sel(i)),
c     .      Sqrt(abs(Real(evals_sel(i)))/2.0),
c     .      Sqrt(abs(Dimag(evals_sel(i)))/2.0),
     .      Real(resf(i))
          end do


         print *, '  '
         print *, '  '

         print *, ' Dimension of Krylov subspace ', ncv
         print *, ' Tolerance ', tol2
         print *, ' Mat*vec operations', iter2
         print *, ' Arnoldi restarts', iter
         print *, ' Converged Ritz values', conv_count
         print *, ' Eigenvectors stored in subfolder "Eigenv/"'
         call CPU_TIME(finish)
         print *, ' Computation time ', finish-start, "seconds"
      
      end if
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

8000  continue  
c      
      end subroutine DIAG_ARNOLDI_MPI
