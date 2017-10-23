         MODULE LR_ARNOLDI_MOD
         IMPLICIT NONE
         SAVE 

        
! OLD tails (shift to main program)
        REAL*8  ::  time_psi_MAX, time_cic_MAX 
        REAL*8  ::  T_From     = 0.d0
         INTEGER,      public :: task
         logical,      public :: binary
         INTEGER,      public :: ncv
         INTEGER,      public :: nev
         INTEGER,      public :: min_conv
         INTEGER,      public :: maxncv
         INTEGER,      public :: maxnev
         INTEGER,      public :: maxn
         REAL*8,    public  ::  tol1     
         REAL*8,    public  ::  tol2     
         REAL*8,    public  ::  upperlimit
         INTEGER*8,    public  ::  maxnonzero

        NAMELIST /LR/T_From,task,binary,ncv,nev,min_conv,tol1,tol2, &
                 upperlimit,maxn,maxncv,maxnev,maxnonzero




!        INTEGER, public ::  MYID,numprocs  

        END MODULE LR_ARNOLDI_MOD



      MODULE  LR_RAPHA

         COMPLEX*16, ALLOCATABLE ::  WSL_mat(:,:,:)
         COMPLEX*16, ALLOCATABLE ::  KSL_mat(:,:,:,:)
         COMPLEX*16, ALLOCATABLE ::  KSL_conjg_mat(:,:,:,:)
         LOGICAL :: orb_real
         Complex*16, ALLOCATABLE :: V_PrdCIJ(:,:,:,:,:)
         Complex*16, ALLOCATABLE :: V_PrdCIJ2(:,:,:)
         COMPLEX*16, ALLOCATABLE :: h2_PSI(:,:)
         COMPLEX*16, ALLOCATABLE :: h2_PSI_conjg(:,:)
         Real*8 :: time_LR
         Complex*16, allocatable :: PSI_LR(:,:)
         Complex*16, allocatable :: InvZRIJ_LR(:,:)
         COMPLEX*16, allocatable :: VIN_LR(:)
         Real*8, allocatable :: Prefactors_1b_Arnoldi(:,:)
         Real*8, allocatable :: Prefactors_2b_Arnoldi(:,:)
         Integer, allocatable :: Ind_CI_1b_Arnoldi(:,:)
         Integer, allocatable :: Ind_CI_2b_Arnoldi(:,:)
         Integer  :: mysize,overallsize
         logical :: preconstr, left_Loo_prec, left_Loc_prec, Tkin_prec
         logical :: h2_prec, WSL_prec, KSL_prec, V_CIJ_prec, proj_prec
         COMPLEX*16, allocatable :: Tkin_2D(:,:)
         COMPLEX*16, allocatable :: Tkin_1D(:,:)
         COMPLEX*16, allocatable :: left_Loo(:,:)
         COMPLEX*16, allocatable :: left_Loc(:,:)
         COMPLEX*16, allocatable :: Loc_vals_csr(:)
         INTEGER*4, allocatable :: Loc_col_csr(:)
         INTEGER*4, allocatable :: Loc_ind_csr(:)
         COMPLEX*16, allocatable :: Loc_vals_coo(:)
         INTEGER*4, allocatable :: Loc_col_coo(:)
         INTEGER*4, allocatable :: Loc_row_coo(:)
         COMPLEX*16, allocatable :: proj_vals_coo(:)
         INTEGER*4, allocatable :: proj_cols_coo(:)
         INTEGER*4, allocatable :: proj_rows_coo(:)
         COMPLEX*16, allocatable :: proj_vals_inter(:)
         INTEGER*4, allocatable :: proj_cols_inter(:)
         INTEGER*4, allocatable :: proj_rows_inter(:)
         COMPLEX*16, allocatable :: projdist_vals_coo(:)
         INTEGER*4, allocatable :: projdist_cols_coo(:)
         INTEGER*4, allocatable :: projdist_rows_coo(:)
         COMPLEX*16, allocatable :: left_Loo_hilf2(:,:) 
         COMPLEX*16, allocatable :: left_Loo_hilf2_p0(:,:) 
         COMPLEX*16, allocatable :: left_Loc_hilf2(:,:) 
         COMPLEX*16, allocatable :: left_Loc_hilf2_p0(:,:) 
         COMPLEX*16, allocatable :: KSL_dist(:,:,:) 
         COMPLEX*16, allocatable :: KSL_conjg_dist(:,:,:) 
         Complex*16, allocatable ::  KSL_hilf2_p0(:)
         Complex*16, allocatable ::  KSL_conjg_hilf2_p0(:)
         Complex*16, allocatable ::  KSL_mat_vals_csc(:,:,:)
         Complex*16, allocatable ::  KSL_conjgmat_vals_csc(:,:,:)
         INTEGER*4, allocatable ::  KSL_mat_ind_csc(:,:,:)
         INTEGER*4, allocatable ::  KSL_conjgmat_ind_csc(:,:,:)
         INTEGER*4, allocatable ::  KSL_mat_rows_csc(:,:,:)
         INTEGER*4, allocatable ::  KSL_conjgmat_rows_csc(:,:,:)
         COMPLEX*16, allocatable :: Tkin_col_p0(:)
         COMPLEX*16, allocatable :: Loo_vals_csr(:)
         INTEGER*4, allocatable :: Loo_col_csr(:)
         INTEGER*4, allocatable :: Loo_ind_csr(:)
         COMPLEX*16, allocatable :: Loo_vals_hilf(:)
         INTEGER*4, allocatable :: Loo_cols_hilf(:)
         INTEGER*4, allocatable :: Loo_rows_hilf(:)
         integer*4, allocatable    ::  proj_nonzeros(:)

         COMPLEX*16, allocatable :: LR_mat(:,:) 
         COMPLEX*16, allocatable :: Loo_u(:,:) 
         COMPLEX*16, allocatable :: Loo_v(:,:) 
          Complex*16, allocatable :: c1_mat(:,:), c2_mat(:,:)
          INTEGER*4, allocatable :: Loc_rows_hilf(:),Loc_cols_hilf(:)
          COMPLEX*16, allocatable :: Loc_vals_hilf(:)

         INTEGER*4 :: proj_nonzero,Loc_nonzeroes

      END MODULE


