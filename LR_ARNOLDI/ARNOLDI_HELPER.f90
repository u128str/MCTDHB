MODULE ARNOLDI_MOD
! USE LINEAR_RESPONSE
! USE Prop_MB
 USE LR_ARNOLDI_MOD
 USE W_INTERPARTICLE
 USE PASS_ARG
 USE SHARED_DIMS
 USE rR_hW
 USE LR_RAPHA
 USE CI_Prod 

      IMPLICIT NONE
         COMPLEX*16, ALLOCATABLE, public :: H_ij(:,:), W_ijkl(:,:,:,:)
         INTEGER, public  :: dimL_orb, dimL, ND, numeig
         INTEGER, DIMENSION(100), public  :: c1_i,c1_j
         INTEGER, DIMENSION(10000), public  :: c2_i,c2_j,c2_k,c2_l

   




      CONTAINS


       subroutine get_Tkin_2D(T_2D)
 
          USE CI_ALL
          USE Parallel_CI 
          USE W_INTERPARTICLE 
          USE DVR_ALL
          USE rR_hW

          IMPLICIT NONE

          COMPLEX*16, DIMENSION(NDX*NDY*NDZ,NDX*NDY*NDZ), INTENT(OUT) :: T_2D
          INTEGER  :: i1,i2,j1,j2,q1,q2,I,J

          Do j1=1,NDY
          Do j2=1,NDY
          IF(j1.eq.j2) then
          q1=1
          else
          q1=0
          endIF
             Do i1=1,NDX
             Do i2=1,NDX
             IF(i1.eq.i2) then
             q2=1
             else
             q2=0
             endIF
             I=NDX*(j1-1)+i1
             J=NDX*(j2-1)+i2
             T_2D(I,J)=q1*Op_X(i1,i2)+q2*Op_Y(j1,j2)
!         IF(I.eq.J) ham(I,I)=ham(I,I)+Vext(I)
!      write(6,'(i2,i2,a10,4i2)') I,J,"i1i2_j1j2",q1*i1,q1*i2,q2*j1,q2*j2
            endDo
            endDo
          endDo
          endDo
!          Tkin_col=ham(:,pos) 

       end subroutine get_Tkin_2D


       subroutine get_KSL_2D(T_2D)
 
          USE CI_ALL
          USE Parallel_CI 
          USE W_INTERPARTICLE 
          USE DVR_ALL
          USE rR_hW

          IMPLICIT NONE

          COMPLEX*16, DIMENSION(NDX*NDY*NDZ,NDX*NDY*NDZ), INTENT(OUT) :: T_2D
          INTEGER  :: i1,i2,j1,j2,q1,q2,I,J

          Do j1=1,NDY
          Do j2=1,NDY
          IF(j1.eq.j2) then
          q1=1
          else
          q1=0
          endIF
             Do i1=1,NDX
             Do i2=1,NDX
             IF(i1.eq.i2) then
             q2=1
             else
             q2=0
             endIF
             I=NDX*(j1-1)+i1
             J=NDX*(j2-1)+i2
             T_2D(I,J)=q1*Op_X(i1,i2)+q2*Op_Y(j1,j2)
!         IF(I.eq.J) ham(I,I)=ham(I,I)+Vext(I)
!      write(6,'(i2,i2,a10,4i2)') I,J,"i1i2_j1j2",q1*i1,q1*i2,q2*j1,q2*j2
            endDo
            endDo
          endDo
          endDo
!          Tkin_col=ham(:,pos) 

       end subroutine get_KSL_2D


       subroutine HPSI_LR(VIN,VOUT,Nc)
       use PASS_ARG
       USE CI_SUBR
       USE CI_ALL
        USE CI_prod 
       implicit NONE
      INTERFACE 
      SUBROUTINE GetCIJKL_All_body_Par(MYID,VIN)
       COMPLEX*16 ::  VIN(:)
       INTEGER    :: MYID
      END SUBROUTINE GetCIJKL_All_body_Par
      END INTERFACE 
!=================== MPI 
       INCLUDE 'mpif.h'
       INTEGER ::  ierr,MYID,numprocs,Nc
!==========================================================
       integer i,j,k,l,P,nn

       integer  maxsil,NL,Iter_Total,Nterms
       INTEGER ::  n,FromN,TillN
!=================== F90 F95

        COMPLEX*16, INTENT(IN) :: VIN(Nc)
        COMPLEX*16, INTENT(OUT) :: VOUT(Nc)
!====================== For SIL
        REAL*8 :: time, Error_SIL, E_state

        COMPLEX*16 :: zrho=ZERO,Z,Z1
        COMPLEX*16 :: Escale
!==============================================
      real*4 start,finish,exec_time ,finish1, iter_time       
      REAL*8 :: start_mpi,start_mpi_all
      REAL*8 :: mpi_time_bc,mpi_time_cp,mpi_time_rd

      LOGICAL  CNV,SIL,DNS

!===========================================================
            start_mpi_all=MPI_WTIME(ierr)
         call MPI_COMM_RANK(MPI_COMM_WORLD,MYID,ierr)
         call MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,ierr)
!===========================================================

             start_mpi=MPI_WTIME(ierr)
!=========================================================
               SIL=.TRUE.
      call MPI_BCAST(SIL,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
!=========================================================


              VOUT   =ZERO
              ZRIJKL =Zero
              ZRIJKL1=Zero
              ZRIJ   =Zero 
              ZRIJ1  =Zero
!==================================================================
!      Broadcast nicht notwendig... 
!      call MPI_BCAST(VIN,Nconf,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)

            start_mpi=MPI_WTIME(ierr)
      !       write(*,*) myid,"bevor"
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)  
      CALL GetCIJKL_All_body_Par(MYID,VIN)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)  
       !      write(*,*) myid,"danach"
      ! write(*,*) myid, size(VIN)

            start_mpi=MPI_WTIME(ierr)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)  
      CALL MPI_REDUCE(VIN,VOUT,Nconf,MPI_DOUBLE_COMPLEX,MPI_SUM,0, &
    &                                           MPI_COMM_WORLD,ierr)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)  

!      nn=Rdim*(Rdim+1)/2
!              ZRIJKL1 =ZRIJKL
!              ZRIJKL =Zero
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)  
!      CALL MPI_REDUCE(ZRIJKL1,ZRIJKL,nn,MPI_DOUBLE_COMPLEX,MPI_SUM,0, &
!    &                                           MPI_COMM_WORLD,ierr)
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)  
!              ZRIJ1 =ZRIJ
!              ZRIJ =Zero
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)  
!      CALL MPI_REDUCE(ZRIJ1,ZRIJ,Rdim,MPI_DOUBLE_COMPLEX,MPI_SUM,0, &
!    &                                           MPI_COMM_WORLD,ierr)
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)  

      end subroutine HPSI_LR





!===CONSTRUCT CI MATRIX ROW==============================================================
       SUBROUTINE construct_CI_row(CI_row,Nc,row)

          USE CI_ALL
          USE Parallel_CI 
          IMPLICIT NONE

          INCLUDE 'mpif.h'
          INTEGER ::  ierr,MYID,numprocs
          COMPLEX*16 :: Vhelp1(Nc), Vhelp2(Nc)
          COMPLEX*16,DIMENSION(2*Nc),INTENT(OUT) :: CI_row
          INTEGER :: row, pos, i, Nc, exists


         call MPI_COMM_RANK(MPI_COMM_WORLD,MYID,ierr)
         call MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,ierr)

          pos=mod(row,Nc)
          if(pos.eq.0) pos=Nc

             Vhelp1=0.0d0
             Vhelp1(pos)=1.0d0
             
!             CALL HPSI(Vhelp1,Vhelp2)
             CALL HPSI_LR(Vhelp1,Vhelp2,Nc)
             Vhelp2=Conjg(Vhelp2)
!             write(*,*) myid, Vhelp2
             Vhelp2(pos)=Vhelp2(pos)-Energy

!             write(*,*) myid, Energy

             CI_row=0d0
             if(row.le.Nc) then 
               CI_row(1:Nc)=Vhelp2
             else
               CI_row(Nc+1:2*Nc)=-Conjg(Vhelp2)
             end if

           

       END SUBROUTINE construct_CI_row




!=========================================================================================================


!===CONSTRUCT ORB MATRIX FULL ROW ORBITAL==============================================================
       SUBROUTINE construct_full_row_orb(PSI,L_orb_row2,row,row0)


          USE CI_ALL
          USE Parallel_CI 
          USE W_INTERPARTICLE 
          USE DVR_ALL
          USE rR_hW
          USE omp_lib 
          USE LR_RAPHA

          IMPLICIT NONE         
          include 'mpif.h'

          INTEGER :: i, j, k, q, s, ll, mm, kk, orb_pos, pos, pos2, orb_pos2, myid, nprocs, ierr 
          INTEGER :: rows_proj

          COMPLEX*16, DIMENSION(ND,Morb), INTENT(IN) :: PSI
!          COMPLEX*16, DIMENSION(dimL_orb), INTENT(OUT) :: L_orb_row
          COMPLEX*16, DIMENSION(Morb,dimL_orb), INTENT(OUT) :: L_orb_row2
          COMPLEX*16, DIMENSION(dimL_orb) :: L_orb_col,Pb,L_orb_col2,Pb2
          COMPLEX*16, DIMENSION(dimL_orb) :: final
          COMPLEX*16, DIMENSION(Morb,dimL_orb) :: inter1,inter2
          COMPLEX*16, DIMENSION(dimL_orb,Morb) :: inter3, final2

          COMPLEX*16, DIMENSION(Morb,dimL_orb) :: left,inter4!,final2
          COMPLEX*16, DIMENSION(dimL_orb,Morb) :: Tinter
          COMPLEX*16, DIMENSION(dimL_orb/2,Morb) :: Loo_u_small
          COMPLEX*16, DIMENSION(dimL_orb/2,Morb) :: Loo_v_small
          COMPLEX*16, DIMENSION(dimL_orb/2,Morb) :: c_u_small
          COMPLEX*16, DIMENSION(dimL_orb/2,Morb) :: c_v_small
          COMPLEX*16, DIMENSION(Morb,Morb) :: Rho_arnoldi
          INTEGER, INTENT(IN)                          :: row          
          INTEGER, INTENT(IN)                          :: row0          
          INTEGER :: s_status(MPI_STATUS_SIZE)
          Real*8 ::  start,finish, t1, t2
          Complex*16 :: alpha,beta,zdotc
          Complex*16, allocatable :: inter(:,:), store(:)
          Complex*16, DIMENSION(dimL_orb) :: sp

          call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
          call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )

      alpha=cmplx(1.0,0.0,kind=8)
      beta=cmplx(0.0,0.0,kind=8)

!        Get pos and orb_pos, important for correct row construction!!! 
          pos=mod(row,ND)
          if(pos.eq.0) pos=ND


         do k=1, 2*Morb
              if(row.le.k*ND) then
               orb_pos=k
               exit
              end if
         end do
         if(orb_pos.gt.Morb) orb_pos=orb_pos-Morb 


!        Get inverse squareroot density matrix
         call squarerootRInv_arnoldi(Rho_arnoldi)
 

         rows_proj=(dimL_orb/2)/nprocs
         if(.not.allocated(store)) allocate(store(rows_proj)) 



         if(left_Loo_prec.eqv..FALSE.) then 



!        Need M projector columns: from the wanted column + (M-1) additional ones
!        YOU CAN ALWAYS GET M ROWS AT ONCE!!!!!!!!!!!!!!!!!

!      call CPU_TIME(start)
           do i=1, Morb
            Pb=0d0
            if(row.le.ND*Morb) then
              call get_proj_col(PSI,Pb,(i-1)*ND+pos,pos)
              left(i,:)=Conjg(Pb)
            else
              call get_proj_col(PSI,Pb,Morb*ND+(i-1)*ND+pos,pos)
              left(i,:)=Conjg(Pb)
            end if 
           end do  
!      call CPU_TIME(finish)
!      write(*,*) "step1", finish-start
           

          
          if(.not.allocated(inter)) allocate(inter(Morb,dimL_orb))
           inter=0d0 
!      call CPU_TIME(start)
!!$OMP PARALLEL
!!$OMP DO
           do j=1, dimL_orb/2
!            write(*,*) "thread",omp_get_thread_num(),"j",j
             L_orb_col=0d0 
             L_orb_col2=0d0 

!       NOT SURE IF THIS IS WORKING LIKE THIS
!          Distribute Tkin_2D for first (ND-nprocs+1) rows
!             if(row.le.ND-nprocs+1) then
              if(dim_mctdhb.eq.2.and.Tkin_prec.eqv..TRUE.) then
               if(.not.allocated(Tkin_col_p0)) allocate(Tkin_col_p0(ND))
               Tkin_col_p0=0d0
               pos2=mod(j,ND)
               if(pos2.eq.0) pos2=ND
               if(myid.eq.0) Tkin_col_p0=Tkin_2D(:,pos2)
               CALL MPI_BCAST(Tkin_col_p0(:),ND,MPI_DOUBLE_COMPLEX,0,0,MPI_COMM_WORLD,ierr) 
              end if
!             end if


             call construct_orb_col(PSI,L_orb_col,j)
      

             L_orb_col2(1:dimL_orb/2)=-Conjg(L_orb_col(dimL_orb/2+1:dimL_orb))
             L_orb_col2(dimL_orb/2+1:dimL_orb)=-Conjg(L_orb_col(1:dimL_orb/2))




!      call CPU_TIME(start)
             do i=1, Morb
              do k=1, ND 
               if(row.le.dimL_orb/2) then
                inter(i,j)=inter(i,j)+left(i,(i-1)*ND+k)*L_orb_col((i-1)*ND+k)
                inter(i,j+Morb*ND)=inter(i,j+Morb*ND)+left(i,(i-1)*ND+k)*L_orb_col2((i-1)*ND+k)
               else
                inter(i,j)=inter(i,j)+left(i,(Morb+i-1)*ND+k)*L_orb_col((Morb+i-1)*ND+k)
                inter(i,j+Morb*ND)=inter(i,j+Morb*ND)+left(i,(Morb+i-1)*ND+k)*L_orb_col2((Morb+i-1)*ND+k)
               end if  
              end do
!               inter(i,j)=dot_product(Conjg(left(i,:)),L_orb_col)
!               inter(i,j+Morb*ND)=dot_product(Conjg(left(i,:)),L_orb_col2)
             end do  
!      call CPU_TIME(finish)
!      if(j.eq.1) write(*,*) "step between 2", finish-start
           end do
!!$OMP END DO
!!$OMP END PARALLEL
!          Now inter contains M rows which look like Proj*L_orb 
!      call CPU_TIME(finish)
!      write(*,*) "step2", finish-start


!      call CPU_TIME(start)
           inter2=0d0 
!!$OMP PARALLEL
!!$OMP DO
           DO j=1, dimL_orb/2
            Pb=0d0
            Pb2=0d0
            pos2=mod(j,ND)
            if(pos2.eq.0) pos2=ND
            call get_proj_col(PSI,Pb,j,pos2)
            call get_proj_col(PSI,Pb2,j+dimL_orb/2,pos2)


            do i=1, 2*Morb
              if(j.le.i*ND) then
               orb_pos2=i
               exit
              end if
            end do
            if(orb_pos2.gt.Morb) orb_pos2=orb_pos2-Morb 

             do i=1, Morb
              do k=1,ND
                 inter2(i,j)=inter2(i,j)+inter(i,(orb_pos2-1)*ND+k)*Pb((orb_pos2-1)*ND+k)        
                 inter2(i,j+dimL_orb/2)=inter2(i,j+dimL_orb/2)+inter(i,(Morb+orb_pos2-1)*ND+k)*Pb2((Morb+orb_pos2-1)*ND+k)        
              end do
             end do  
           END DO
!!$OMP END DO
!!$OMP END PARALLEL
!          Now middle block Proj*L_orb*Proj done for Morb rows at once, stored in array inter2(:,:)  !!          
!      call CPU_TIME(finish)
!      write(*,*) "step3", finish-start




!      call CPU_TIME(start)
!        Rho from the left
          inter4=0d0
          if(row.le.Morb*ND) then
            do i=1, Morb 
             do k=1, Morb
              do j=1, ND
               do s=1, Morb
                inter4(i,j+(k-1)*ND)=inter4(i,j+(k-1)*ND)+Rho_arnoldi(i,s)*inter2(s,j+(k-1)*ND)
                inter4(i,Morb*ND+j+(k-1)*ND)=inter4(i,Morb*ND+j+(k-1)*ND)+Rho_arnoldi(i,s)*inter2(s,j+Morb*ND+(k-1)*ND)
               end do
              end do 
             end do 
            end do
          else
            do i=1, Morb 
             do k=1, Morb
              do j=1, ND
               do s=1, Morb
                inter4(i,j+(k-1)*ND)=inter4(i,j+(k-1)*ND)+Conjg(Rho_arnoldi(i,s))*inter2(s,j+(k-1)*ND)
                inter4(i,Morb*ND+j+(k-1)*ND)=inter4(i,Morb*ND+j+(k-1)*ND)+Conjg(Rho_arnoldi(i,s))*inter2(s,j+Morb*ND+(k-1)*ND)
               end do
              end do 
             end do 
            end do
          end if 
!      call CPU_TIME(finish)
!      write(*,*) "step4", finish-start




!      call CPU_TIME(start)
!        Rho from the right
          final2=0d0
           do s=1, Morb   
            do i=1, Morb 
             do j=1, ND
              do k=1, Morb
                final2(s,(i-1)*ND+j)=final2(s,(i-1)*ND+j)+Rho_arnoldi(k,i)*inter4(s,(k-1)*ND+j)
                final2(s,Morb*ND+(i-1)*ND+j)=final2(s,Morb*ND+(i-1)*ND+j)+Conjg(Rho_arnoldi(k,i))*inter4(s,Morb*ND+(k-1)*ND+j)
              end do 
             end do 
            end do
           end do
!      call CPU_TIME(finish)
!      write(*,*) "step5", finish-start


!           L_orb_row=final
           do i=1, Morb
             L_orb_row2(i,:)=final2(i,:)
           end do
 


       else

        IF(row0.le.ND-nprocs+1) THEN
!            inter3=left_Loo_hilf2
!          goto 19

         Tinter=Transpose(left_Loo_hilf2)
         do i=1, Morb       
            Loo_u(1:dimL_orb/2,(i-1)*nprocs+myid+1)=Tinter(1:dimL_orb/2,i) 
            Loo_v(1:dimL_orb/2,(i-1)*nprocs+myid+1)=Tinter(dimL_orb/2+1:dimL_orb,i) 
         end do  


         do k=1, nprocs
          do j=1, Morb
            call MPI_BCAST(Loo_u(:,1+(k-1)+(j-1)*nprocs), &
     &                dimL_orb/2, &
     &                MPI_DOUBLE_COMPLEX,k-1,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(Loo_v(:,1+(k-1)+(j-1)*nprocs), &
     &                dimL_orb/2, &
     &                MPI_DOUBLE_COMPLEX,k-1,MPI_COMM_WORLD,ierr)
          end do
         end do


       c1_mat=0d0
       c2_mat=0d0
        if(myid.eq.0) then
               do i=1, nprocs*Morb
                call mkl_zcoogemv('N', rows_proj+mod(dimL_orb/2,nprocs), projdist_vals_coo, projdist_rows_coo, &
     &                   projdist_cols_coo, proj_nonzeros(1), Loo_u(:,i), &
     &                   c1_mat(:,i))

                call mkl_zcoogemv('N', rows_proj+mod(dimL_orb/2,nprocs), Conjg(projdist_vals_coo), projdist_rows_coo, &
     &                   projdist_cols_coo, proj_nonzeros(1), Loo_v(:,i), &
     &                   c2_mat(:,i))
               end do 
        else
               do i=1, nprocs*Morb
                  call mkl_zcoogemv('N', rows_proj, projdist_vals_coo, projdist_rows_coo, &
     &                   projdist_cols_coo, proj_nonzeros(myid+1), Loo_u(:,i), &
     &                   c1_mat(:,i))

                  call mkl_zcoogemv('N', rows_proj, Conjg(projdist_vals_coo), projdist_rows_coo, &
     &                   projdist_cols_coo, proj_nonzeros(myid+1), Loo_v(:,i), &
     &                   c2_mat(:,i))
               end do
        end if


       Loo_u=0d0
       Loo_v=0d0
       if(myid.eq.0) then
          do i=1, nprocs*Morb
              Loo_u(1: &
     &                   mod(dimL_orb/2,nprocs)+rows_proj,i)=c1_mat(1: &
     &                   mod(dimL_orb/2,nprocs)+rows_proj,i)

              Loo_v(1: &
     &                   mod(dimL_orb/2,nprocs)+rows_proj,i)=c2_mat(1: &
     &                   mod(dimL_orb/2,nprocs)+rows_proj,i)
          end do
       end if
         do i=1, nprocs*Morb
          call MPI_BCAST(Loo_u(1: &
     &                   mod(dimL_orb/2,nprocs)+rows_proj,i), &
     &                rows_proj+mod(dimL_orb/2,nprocs), &
     &                MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)

          call MPI_BCAST(Loo_v(1: &
     &                   mod(dimL_orb/2,nprocs)+rows_proj,i), &
     &                rows_proj+mod(dimL_orb/2,nprocs), &
     &                MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
         end do
       do k=2, nprocs
            if(myid.eq.k-1) then
             do i=1, nprocs*Morb
              Loo_u(1+mod(dimL_orb/2,nprocs)+(k-1)*rows_proj: &
     &                   mod(dimL_orb/2,nprocs)+k*rows_proj,i)= & 
     &             c1_mat(1: &
     &                   rows_proj,i)

              Loo_v(1+mod(dimL_orb/2,nprocs)+(k-1)*rows_proj: &
     &                   mod(dimL_orb/2,nprocs)+k*rows_proj,i)= & 
     &             c2_mat(1: &
     &                   rows_proj,i)
             end do
            end if 
            do i=1, nprocs*Morb
             call MPI_BCAST(Loo_u(1+mod(dimL_orb/2,nprocs)+(k-1)*rows_proj: &
     &                   mod(dimL_orb/2,nprocs)+k*rows_proj,i), &
     &                rows_proj, &
     &                MPI_DOUBLE_COMPLEX,k-1,MPI_COMM_WORLD,ierr)

             call MPI_BCAST(Loo_v(1+mod(dimL_orb/2,nprocs)+(k-1)*rows_proj: &
     &                   mod(dimL_orb/2,nprocs)+k*rows_proj,i), &
     &                rows_proj, &
     &                MPI_DOUBLE_COMPLEX,k-1,MPI_COMM_WORLD,ierr)
            end do
       end do 




            
           do i=1, Morb
             inter3(1:dimL_orb/2,i)=Loo_u(:,(i-1)*nprocs+myid+1)
             inter3(dimL_orb/2+1:dimL_orb,i)=Loo_v(:,(i-1)*nprocs+myid+1)
           end do

!19         continue            




!        Rho from the right
          final2=0d0
           do s=1, Morb   
            do i=1, Morb 
             do j=1, ND
              do k=1, Morb
                final2((i-1)*ND+j,s)=final2((i-1)*ND+j,s)+Rho_arnoldi(k,i)*inter3((k-1)*ND+j,s)
                final2(Morb*ND+(i-1)*ND+j,s)=final2(Morb*ND+(i-1)*ND+j,s)+Conjg(Rho_arnoldi(k,i))*inter3(Morb*ND+(k-1)*ND+j,s)
              end do 
             end do 
            end do
           end do


!           L_orb_row=final
           do i=1, Morb
             L_orb_row2(i,:)=final2(:,i)
           end do
        
        ELSE
 

!            inter3=left_Loo_hilf2
!          goto 20

          if(myid.eq.0) then
           Tinter=Transpose(left_Loo_hilf2)
           do i=1, Morb       
            Loo_u_small(1:dimL_orb/2,i)=Tinter(1:dimL_orb/2,i) 
            Loo_v_small(1:dimL_orb/2,i)=Tinter(dimL_orb/2+1:dimL_orb,i) 
           end do  
          end if

            do j=1, Morb
            call MPI_BCAST(Loo_u_small(:,j), &
     &                dimL_orb/2, &
     &                MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(Loo_v_small(:,j), &
     &                dimL_orb/2, &
     &                MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
            end do



       c_u_small=0d0
       c_v_small=0d0
        if(myid.eq.0) then
               do i=1, Morb
                call mkl_zcoogemv('N', rows_proj+mod(dimL_orb/2,nprocs), projdist_vals_coo, projdist_rows_coo, &
     &                   projdist_cols_coo, proj_nonzeros(1), Loo_u_small(:,i), &
     &                   c_u_small(1:rows_proj+mod(dimL_orb/2,nprocs),i))

                call mkl_zcoogemv('N', rows_proj+mod(dimL_orb/2,nprocs), Conjg(projdist_vals_coo), projdist_rows_coo, &
     &                   projdist_cols_coo, proj_nonzeros(1), Loo_v_small(:,i), &
     &                   c_v_small(1:rows_proj+mod(dimL_orb/2,nprocs),i))
               end do

        else      
               do i=1, Morb
                call mkl_zcoogemv('N', rows_proj, projdist_vals_coo, projdist_rows_coo, &
     &                   projdist_cols_coo, proj_nonzeros(myid+1), Loo_u_small(:,i), &
     &                   c_u_small(1+myid*rows_proj+mod(dimL_orb/2,nprocs):(myid+1)*rows_proj+mod(dimL_orb/2,nprocs),i))

                call mkl_zcoogemv('N', rows_proj, Conjg(projdist_vals_coo), projdist_rows_coo, &
     &                   projdist_cols_coo, proj_nonzeros(myid+1), Loo_v_small(:,i), &
     &                   c_v_small(1+myid*rows_proj+mod(dimL_orb/2,nprocs):(myid+1)*rows_proj+mod(dimL_orb/2,nprocs),i))
               end do
        end if

      
           inter3=0d0
 
            if(myid.eq.0) then
              do j=1, Morb
                 inter3(1:rows_proj+mod(dimL_orb/2,nprocs),j)=c_u_small(1:rows_proj+mod(dimL_orb/2,nprocs),j)
                 inter3(1+dimL_orb/2:dimL_orb/2+rows_proj+mod(dimL_orb/2,nprocs),j)=c_v_small(1:rows_proj+mod(dimL_orb/2,nprocs),j)
              end do
            end if

           do k=2, nprocs
            do j=1, Morb
             if(myid.eq.k-1) then
                  store(:)= &  
     &            c_u_small(1+(k-1)*rows_proj+mod(dimL_orb/2,nprocs):k*rows_proj+mod(dimL_orb/2,nprocs),j) 

              CALL MPI_SEND(store,rows_proj, &
     &         MPI_DOUBLE_COMPLEX,0,0,MPI_COMM_WORLD,ierr)
             end if 
!              CALL MPI_BCAST(inter3(1+(k-1)*rows_proj+mod(dimL_orb/2,nprocs):k*rows_proj+mod(dimL_orb/2,nprocs),j),rows_proj, &
!     &         MPI_DOUBLE_COMPLEX,k-1,MPI_COMM_WORLD,ierr)

             if(myid.eq.0) then
              CALL MPI_RECV(store,rows_proj, &
     &         MPI_DOUBLE_COMPLEX,k-1,0,MPI_COMM_WORLD,s_status,ierr)

              inter3(1+(k-1)*rows_proj+mod(dimL_orb/2,nprocs):k*rows_proj+mod(dimL_orb/2,nprocs),j)= &
     &          store(:)       
             end if

             if(myid.eq.k-1) then
                  store(:)= &
     &            c_v_small(1+(k-1)*rows_proj+mod(dimL_orb/2,nprocs):k*rows_proj+mod(dimL_orb/2,nprocs),j) 

              CALL MPI_SEND(store,rows_proj, &
     &         MPI_DOUBLE_COMPLEX,0,1,MPI_COMM_WORLD,ierr)

             end if
!              CALL MPI_BCAST(inter3(1+dimL_orb/2+(k-1)*rows_proj+mod(dimL_orb/2,nprocs):dimL_orb/2+k*rows_proj+mod(dimL_orb/2,nprocs),j),rows_proj, &
!     &         MPI_DOUBLE_COMPLEX,k-1,MPI_COMM_WORLD,ierr)

             if(myid.eq.0) then
                CALL MPI_RECV(store,rows_proj, &
     &         MPI_DOUBLE_COMPLEX,k-1,1,MPI_COMM_WORLD,s_status,ierr)

              inter3(1+dimL_orb/2+(k-1)*rows_proj+mod(dimL_orb/2,nprocs):dimL_orb/2+k*rows_proj+mod(dimL_orb/2,nprocs),j)= &
     &           store(:)
             end if

            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            end do
           end do

!20        continue


           if(myid.eq.0) then
         !        Rho from the right
             final2=0d0
               do s=1, Morb   
                do i=1, Morb 
                 do j=1, ND
                  do k=1, Morb
                   final2((i-1)*ND+j,s)=final2((i-1)*ND+j,s)+Rho_arnoldi(k,i)*inter3((k-1)*ND+j,s)
                   final2(Morb*ND+(i-1)*ND+j,s)=final2(Morb*ND+(i-1)*ND+j,s)+Conjg(Rho_arnoldi(k,i))*inter3(Morb*ND+(k-1)*ND+j,s)
                  end do 
                 end do 
                end do
               end do

   !           L_orb_row=final
              do i=1, Morb
               L_orb_row2(i,:)=final2(:,i)
              end do

           end if
           call MPI_BARRIER(MPI_COMM_WORLD,ierr) 

        END IF 
 

       end if   ! end if of preconstruction left part upfront





!           open(unit=256, file='test.dat', status='replace', action='readwrite')
!           do i=1, dimL_orb
!             write(256,'(500(2F20.10))') Real(L_orb_row2(1,i)),Dimag(L_orb_row2(1,i))

!           end do
!          close(256)

       END SUBROUTINE construct_full_row_orb



!====== CONSTRUCT heavy left side of Loo upfront ONLY ONCE!!!!!!
       SUBROUTINE construct_left_of_Loo_upfront_sparse(PSI)


          USE CI_ALL
          USE Parallel_CI 
          USE W_INTERPARTICLE 
          USE DVR_ALL
          USE rR_hW
          USE omp_lib 
          USE LR_RAPHA
          USE LR_ARNOLDI_MOD

          IMPLICIT NONE         
          include 'mpif.h'

          INTEGER*4 :: m, i, j, k, q, s, ll, mm, kk, orb_pos, pos, pp
          INTEGER*4 :: pos2, orb_pos2, counter, infoo, start_from, jmax, myid, nprocs, ierr 
          INTEGER, DIMENSION(2) :: ind11, ind12, ind21, ind22
          INTEGER, DIMENSION(2) :: ind11a, ind12a, ind21a, ind22a
!          INTEGER*8 :: maxnonzero
          COMPLEX*16, DIMENSION(ND,Morb), INTENT(IN) :: PSI
!          COMPLEX*16, DIMENSION(dimL_orb,dimL_orb) :: left_Loo
          COMPLEX*16, allocatable :: L_orb_help2(:,:), inter_Loo(:,:)
          COMPLEX*16, DIMENSION(dimL_orb) :: L_orb_col,Pb,L_orb_col2
          Complex*16, allocatable ::  mat_val(:)
          Complex*16, allocatable ::  vals_inter(:)
          Complex*16, allocatable ::  inter_mat(:,:),c_mat(:,:),inter2_mat(:,:)
          Complex*16, allocatable ::  inter_mat_small(:,:),inter2_mat_small(:,:)
          COMPLEX*16, allocatable :: left_Loo_hilf(:,:),left_Loo_hilf_p0(:,:)
!          COMPLEX*16, DIMENSION(dimL_orb/2,2) :: left_Loo_hilf2
          integer*4, allocatable    ::  mat_csr1(:)
          integer*4, allocatable    ::  mat_csr2(:)
          integer*4, allocatable    ::  rows_inter(:)
          integer*4, allocatable    ::  cols_inter(:)
          integer, dimension(6)    ::  mjob
          COMPLEX*16, DIMENSION(Morb,Morb) :: RhoInvSQR
          Real*8 ::  start,finish,a,b,tol1_Loo, t1, t2, szgb
          logical :: exists
          INTEGER :: s_status(MPI_STATUS_SIZE), st, ii, jj, info, rows_proj
          INTEGER, allocatable :: ind(:)
          COMPLEX*16 :: alpha, beta, zdotc 
          INTEGER*4 :: Loo_maxnonzero
          CHARACTER(1), dimension(6) :: matdescra

          call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
          call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )

          allocate(ind(nprocs))
!!          Loo_maxnonzero=1500000000
          Loo_maxnonzero=maxnonzero

          alpha=Cmplx(1.0,0.0,kind=8)
          beta=Cmplx(0.0,0.0,kind=8)

          counter=0;start_from=1;jmax=1;ind=1; 

          tol1_Loo=tol1
!          tol1_Loo=1e-5
          rows_proj=(dimL_orb/2)/nprocs

          call squarerootRInv_arnoldi(RhoInvSQR)

          IF(myid.eq.0) THEN 

            if(.not.allocated(inter_mat)) &
     &       allocate(inter_mat(dimL_orb,2*nprocs))
            if(.not.allocated(inter_mat_small)) &
     &       allocate(inter_mat_small(dimL_orb,2))
            if(.not.allocated(inter2_mat_small)) &
     &       allocate(inter2_mat_small(dimL_orb/2,2))
            if(.not.allocated(inter_Loo)) &
     &       allocate(inter_Loo(dimL_orb/2,2*nprocs))


          write(*,*)"STR 5 ok"

          INQUIRE(FILE="left_Loo.dat", EXIST=exists)

          if(exists.eqv..FALSE.) then  !files do not exist       
 
            open(unit=256, file='left_Loo.dat', & 
     &        status='replace', action='readwrite')

            if(.not.allocated(Loo_vals_hilf)) then
          allocate(Loo_vals_hilf(Loo_maxnonzero), stat=ierr)
          szgb= Loo_maxnonzero*32.0d0/1024/1024/1024
          write(*,*)"Size Loo_vals_hilf Gigabytes:", szgb, Loo_maxnonzero
          if(ierr /= 0) then
          write(*,*)"allocation error for Loo_vals_hilf"
          stop "Err alloc Loo_vals_hilf"
          else
          if(ierr == 0)write(*,*)"allocation ok for Loo_vals_hilf"
          endif
           endif

            if(.not.allocated(Loo_rows_hilf)) then
            allocate(Loo_rows_hilf(Loo_maxnonzero), stat=ierr)
          szgb= Loo_maxnonzero*16*2.0/1024/1024/1024
          write(*,*)"Size Loo_rows_hilf Gigabytes:", szgb
          if(ierr /= 0) then
          write(*,*)"allocation error for Loo_rows_hilf"
          stop "Err alloc Loo_rows_hilf"
          else
          if(ierr == 0)write(*,*)"allocation ok for Loo_rows_hilf"
          endif
            endif
 
            if(.not.allocated(Loo_cols_hilf)) then
           allocate(Loo_cols_hilf(Loo_maxnonzero), stat=ierr)
          szgb= Loo_maxnonzero*16*2.0/1024/1024/1024
          write(*,*)"Size Loo_cols_hilf Gigabytes:", szgb
          if(ierr /= 0) then
          write(*,*)"allocation error for Loo_cols_hilf"
          stop "Err alloc Loo_cols_hilf"
          else
          if(ierr == 0)write(*,*)"allocation ok for Loo_cols_hilf"
          endif
            endif


            start_from=1
            ind(1)=1

          else

          write(*,*)"STR 6"
            open(unit=256, file='left_Loo.dat', &
     &        status='old', action='readwrite')

!!            if(.not.allocated(Loo_vals_hilf)) &
!!     &       allocate(Loo_vals_hilf(Loo_maxnonzero))
            if(.not.allocated(Loo_vals_hilf)) then
          allocate(Loo_vals_hilf(Loo_maxnonzero), stat=ierr)
          szgb= Loo_maxnonzero*32.0/1024/1024/1024
          write(*,*)"Size Loo_vals_hilf Gigabytes:", szgb, Loo_maxnonzero, maxnonzero
          if(ierr /= 0) then
          write(*,*)"allocation error for Loo_vals_hilf"
          stop "Err alloc Loo_vals_hilf"
          else
          if(ierr == 0)write(*,*)"allocation ok for Loo_vals_hilf"
          endif
           endif
!!            if(.not.allocated(Loo_rows_hilf)) &
!!     &       allocate(Loo_rows_hilf(Loo_maxnonzero))
            if(.not.allocated(Loo_rows_hilf)) then
            allocate(Loo_rows_hilf(Loo_maxnonzero), stat=ierr)
          szgb= Loo_maxnonzero*32.0/1024/1024/1024
          write(*,*)"Size Loo_rows_hilf Gigabytes:", szgb
          if(ierr /= 0) then
          write(*,*)"allocation error for Loo_rows_hilf"
          stop "Err alloc Loo_rows_hilf"
          else
          if(ierr == 0)write(*,*)"allocation ok for Loo_rows_hilf"
          endif
            endif
!!            if(.not.allocated(Loo_cols_hilf)) &
!!     &       allocate(Loo_cols_hilf(Loo_maxnonzero))
            if(.not.allocated(Loo_cols_hilf)) then
           allocate(Loo_cols_hilf(Loo_maxnonzero), stat=ierr)
          szgb= Loo_maxnonzero*32.0/1024/1024/1024
          write(*,*)"Size Loo_cols_hilf Gigabytes:", szgb
          if(ierr /= 0) then
          write(*,*)"allocation error for Loo_cols_hilf"
          stop "Err alloc Loo_cols_hilf"
          else
          if(ierr == 0)write(*,*)"allocation ok for Loo_cols_hilf"
          endif
            endif

  
            jmax=1; ind=1;
            do 
              read(256,*,iostat=st) i,j,a,b
              if(st.ne.0) exit
               
              Loo_vals_hilf(ind)=Cmplx(a,b,kind=8)
              Loo_rows_hilf(ind)=i
              Loo_cols_hilf(ind)=j

              if(j.le.dimL_orb/2) then 
                if(j.gt.jmax) jmax=j  
              end if
              ind(1)=ind(1)+1  
            end do



            start_from=jmax 

          end if
          END IF


          CALL MPI_BCAST(exists,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(start_from,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

          if(start_from.eq.dimL_orb/2) goto 658
          if(.not.allocated(KSL_hilf2_p0)) allocate(KSL_hilf2_p0(ND)) 
          if(.not.allocated(KSL_conjg_hilf2_p0)) allocate(KSL_conjg_hilf2_p0(ND)) 
          if(.not.allocated(left_Loo_hilf2)) &
     &                     allocate(left_Loo_hilf2(dimL_orb/2,2))
          if(.not.allocated(left_Loo_hilf)) &
     &                     allocate(left_Loo_hilf(dimL_orb/2,2))
          if(.not.allocated(left_Loo_hilf_p0)) &
     &                     allocate(left_Loo_hilf_p0(dimL_orb/2,2))

            if(.not.allocated(inter2_mat)) &
     &       allocate(inter2_mat(dimL_orb/2,2*nprocs))

            if(myid.eq.0) then
             if(.not.allocated(c_mat)) &
     &        allocate(c_mat(rows_proj+mod(dimL_orb/2,nprocs),2*nprocs))
            else 
             if(.not.allocated(c_mat)) &
     &        allocate(c_mat(rows_proj,2*nprocs))
            end if            

            matdescra(1)='H'
            matdescra(2)='L'
            matdescra(3)='N'
            matdescra(4)='F'

          if(.not.allocated(vals_inter)) allocate(vals_inter(dimL_orb)) 
          if(.not.allocated(rows_inter)) allocate(rows_inter(dimL_orb)) 
          if(.not.allocated(cols_inter)) allocate(cols_inter(dimL_orb)) 

!         Calculate Proj*L_orb
           j=start_from 
           do while(j.le.dimL_orb/2-nprocs+1) 
             L_orb_col=0d0 
             L_orb_col2=0d0
             left_Loo_hilf=0d0 
             left_Loo_hilf2=0d0 
             vals_inter=0d0
             cols_inter=0
             rows_inter=0
             inter2_mat=0d0
             c_mat=0d0


       call CPU_TIME(t1)

!            BROADCAST KSL entries to each process
             if(KSL_prec.eqv..TRUE.) then
              KSL_dist=0d0
              KSL_conjg_dist=0d0

             IF(myid.eq.0) then
              pos=mod(j,ND)
              if(pos.eq.0) pos=ND
              do mm=1, Morb
               do ll=1, Morb
                  k=KSL_mat_ind_csc(pos,mm,ll)
                  if(k.ne.KSL_mat_ind_csc(pos+1,mm,ll)) then
                  pp=pos+1
                  do k=KSL_mat_ind_csc(pos,mm,ll), KSL_mat_ind_csc(pp,mm,ll)-1 
                    KSL_dist(KSL_mat_rows_csc(k,mm,ll),mm,ll)= &
     &                           KSL_mat_vals_csc(k,mm,ll)
                  end do 
                  end if
                  k=KSL_conjgmat_ind_csc(pos,mm,ll)
                  if(k.ne.KSL_conjgmat_ind_csc(pos+1,mm,ll)) then
                  pp=pos+1
                  do k=KSL_conjgmat_ind_csc(pos,mm,ll), &
     &                             KSL_conjgmat_ind_csc(pp,mm,ll)-1 
                    KSL_conjg_dist(KSL_conjgmat_rows_csc(k,mm,ll),mm,ll)= &
     &                           KSL_conjgmat_vals_csc(k,mm,ll)
                  end do 
                  end if
               end do
              end do 




             END IF




             do m=2, nprocs
              call MPI_BARRIER(MPI_COMM_WORLD,ierr)
              pos2=mod(j+m-1,ND)
              if(pos2.eq.0) pos2=ND
              do mm=1, Morb
                do ll=1, Morb
                 KSL_hilf2_p0=0d0
                 KSL_conjg_hilf2_p0=0d0
                 if(myid.eq.0) then
                  k=KSL_mat_ind_csc(pos2,mm,ll)
                  if(k.ne.KSL_mat_ind_csc(pos2+1,mm,ll)) then
                  pp=pos2+1
                  do k=KSL_mat_ind_csc(pos2,mm,ll), KSL_mat_ind_csc(pp,mm,ll)-1 
                     KSL_hilf2_p0(KSL_mat_rows_csc(k,mm,ll))= &
     &                           KSL_mat_vals_csc(k,mm,ll)
                  end do 
                  end if
                  k=KSL_conjgmat_ind_csc(pos2,mm,ll)
                  if(k.ne.KSL_conjgmat_ind_csc(pos2+1,mm,ll)) then
                  pp=pos2+1
                  do k=KSL_conjgmat_ind_csc(pos2,mm,ll), &
     &                             KSL_conjgmat_ind_csc(pp,mm,ll)-1 
                     KSL_conjg_hilf2_p0(KSL_conjgmat_rows_csc(k,mm,ll))= &
     &                           KSL_conjgmat_vals_csc(k,mm,ll)
                  end do 
                  end if
                 end if

                  if(myid.eq.m-1) CALL MPI_RECV(KSL_hilf2_p0(:),ND,MPI_DOUBLE_COMPLEX,0,0,MPI_COMM_WORLD,s_status,ierr) 
                  if(myid.eq.0) CALL MPI_SEND(KSL_hilf2_p0(:),ND,MPI_DOUBLE_COMPLEX,m-1,0,MPI_COMM_WORLD,ierr)
                  if(myid.eq.m-1) CALL MPI_RECV(KSL_conjg_hilf2_p0(:),ND,MPI_DOUBLE_COMPLEX,0,1,MPI_COMM_WORLD,s_status,ierr) 
                  if(myid.eq.0) CALL MPI_SEND(KSL_conjg_hilf2_p0(:),ND,MPI_DOUBLE_COMPLEX,m-1,1,MPI_COMM_WORLD,ierr)
                  if(myid.eq.m-1) KSL_dist(:,mm,ll)=KSL_hilf2_p0(:)  
                  if(myid.eq.m-1) KSL_conjg_dist(:,mm,ll)=KSL_conjg_hilf2_p0(:) 
                end do
              end do
             end do

             end if  


!            BROADCAST Tkin_2D entries to each process
             if(dim_mctdhb.eq.2.and.Tkin_prec.eqv..TRUE.) then
               if(.not.allocated(Tkin_col_p0)) allocate(Tkin_col_p0(ND))
               Tkin_col_p0=0d0
               do m=2, nprocs
                 pos2=mod(j+m-1,ND)
                 if(pos2.eq.0) pos2=ND
                 if(myid.eq.0) Tkin_col_p0=Tkin_2D(:,pos2)
                 if(myid.eq.m-1) CALL MPI_RECV(Tkin_col_p0(:),ND,MPI_DOUBLE_COMPLEX,0,0,MPI_COMM_WORLD,s_status,ierr) 
                 if(myid.eq.0) CALL MPI_SEND(Tkin_col_p0(:),ND,MPI_DOUBLE_COMPLEX,m-1,0,MPI_COMM_WORLD,ierr)
                 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
               end do 
               if(myid.eq.0) then 
                 pos2=mod(j,ND)
                 if(pos2.eq.0) pos2=ND
                 Tkin_col_p0=Tkin_2D(:,pos2)
               end if 
             end if
       call CPU_TIME(t2)


             call MPI_BARRIER(MPI_COMM_WORLD,ierr)


       call CPU_TIME(t1)

!            CALL ORBITAL COLUMN CONSTRUCTION
             call construct_orb_col(PSI,L_orb_col,j+myid)
             L_orb_col2(1:dimL_orb/2)=-Conjg(L_orb_col(dimL_orb/2+1:dimL_orb))

       call CPU_TIME(t2)




       call CPU_TIME(t1)
       inter2_mat=0d0
       inter2_mat(1:dimL_orb/2,1+myid)=L_orb_col(1:dimL_orb/2)  
       inter2_mat(1:dimL_orb/2,1+myid+nprocs)=L_orb_col2(1:dimL_orb/2)  

       do k=1, nprocs
            call MPI_BCAST(inter2_mat(:,1+(k-1)), &
     &                dimL_orb/2, &
     &                MPI_DOUBLE_COMPLEX,k-1,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(inter2_mat(:,1+nprocs+(k-1)), &
     &                dimL_orb/2, &
     &                MPI_DOUBLE_COMPLEX,k-1,MPI_COMM_WORLD,ierr)
       end do
       call CPU_TIME(t2)


       c_mat=0d0
       call CPU_TIME(t1)
       if(myid.eq.0) then
               do i=1, 2*nprocs
                call mkl_zcoogemv('N', rows_proj+mod(dimL_orb/2,nprocs), projdist_vals_coo, projdist_rows_coo, &
     &                   projdist_cols_coo, proj_nonzeros(1), inter2_mat(:,i), &
     &                   c_mat(:,i))
               end do 


       else


         do i=1, 2*nprocs
                call mkl_zcoogemv('N', rows_proj, projdist_vals_coo, projdist_rows_coo, &
     &                   projdist_cols_coo, proj_nonzeros(myid+1), inter2_mat(:,i), &
     &                   c_mat(:,i))
         end do

       end if
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call CPU_TIME(t2)

       call CPU_TIME(t1)
       inter2_mat=0d0
       if(myid.eq.0) then
          do i=1, 2*nprocs
              inter2_mat(1: &
     &                   mod(dimL_orb/2,nprocs)+rows_proj,i)=c_mat(1: &
     &                   mod(dimL_orb/2,nprocs)+rows_proj,i)
          end do
       end if
         do i=1, 2*nprocs
          call MPI_BCAST(inter2_mat(1: &
     &                   mod(dimL_orb/2,nprocs)+rows_proj,i), &
     &                rows_proj+mod(dimL_orb/2,nprocs), &
     &                MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
         end do
       do k=2, nprocs
            if(myid.eq.k-1) then
             do i=1, 2*nprocs
              inter2_mat(1+mod(dimL_orb/2,nprocs)+(k-1)*rows_proj: &
     &                   mod(dimL_orb/2,nprocs)+k*rows_proj,i)= & 
     &             c_mat(1: &
     &                   rows_proj,i)
             end do
            end if 
            do i=1, 2*nprocs
             call MPI_BCAST(inter2_mat(1+mod(dimL_orb/2,nprocs)+(k-1)*rows_proj: &
     &                   mod(dimL_orb/2,nprocs)+k*rows_proj,i), &
     &                rows_proj, &
     &                MPI_DOUBLE_COMPLEX,k-1,MPI_COMM_WORLD,ierr)
            end do
       end do 
       call CPU_TIME(t2)
       call MPI_Barrier(MPI_COMM_WORLD,ierr)


       left_Loo_hilf(:,1)=inter2_mat(:,1+myid)
       left_Loo_hilf(:,2)=inter2_mat(:,1+myid+nprocs)
       left_Loo_hilf2=0d0 


         do k=2, nprocs
              ind(k)=1
         end do
       
       call CPU_TIME(t1)
       do i=1, dimL_orb/2 
             pos=mod(i,ND)
             if(pos.eq.0) pos=ND
             do kk=1, Morb
              if(i.le.kk*ND) then
               orb_pos=kk
               exit
              end if
             end do
               do jj=1, Morb
                 left_Loo_hilf2(i,1)=left_Loo_hilf2(i,1)+ &
     &                         RhoInvSQR(orb_pos,jj)*left_Loo_hilf(pos+(jj-1)*ND,1)
                 left_Loo_hilf2(i,2)=left_Loo_hilf2(i,2)+ &
     &                         RhoInvSQR(orb_pos,jj)*left_Loo_hilf(pos+(jj-1)*ND,2)
               end do 

             if(myid.eq.0) then
               if(abs(left_Loo_hilf2(i,1)).gt.tol1_Loo) then
                 Loo_vals_hilf(ind(1))=left_Loo_hilf2(i,1)
                 Loo_cols_hilf(ind(1))=j
                 Loo_rows_hilf(ind(1))=i
                 ind(1)=ind(1)+1
               end if
               if(abs(left_Loo_hilf2(i,2)).gt.tol1_Loo) then
                 Loo_vals_hilf(ind(1))=left_Loo_hilf2(i,2)
                 Loo_cols_hilf(ind(1))=j+dimL_orb/2
                 Loo_rows_hilf(ind(1))=i
                 ind(1)=ind(1)+1
               end if
             else
               if(abs(left_Loo_hilf2(i,1)).gt.tol1_Loo) then
                 vals_inter(ind(myid+1))=left_Loo_hilf2(i,1)
                 cols_inter(ind(myid+1))=j+myid
                 rows_inter(ind(myid+1))=i
                 ind(myid+1)=ind(myid+1)+1
               end if
               if(abs(left_Loo_hilf2(i,2)).gt.tol1_Loo) then
                 vals_inter(ind(myid+1))=left_Loo_hilf2(i,2)
                 cols_inter(ind(myid+1))=j+myid+dimL_orb/2
                 rows_inter(ind(myid+1))=i
                 ind(myid+1)=ind(myid+1)+1
               end if
             end if

       end do


            do k=2, nprocs
              if(myid.eq.0) CALL MPI_RECV(ind(k),1,MPI_INTEGER,k-1,2,MPI_COMM_WORLD,s_status,ierr) 
              if(myid.eq.k-1) CALL MPI_SEND(ind(k),1,MPI_INTEGER,0,2,MPI_COMM_WORLD,ierr)

             if(ind(k).ne.1) then
              if(myid.eq.0) CALL MPI_RECV(vals_inter(1:ind(k)-1),ind(k)-1,MPI_DOUBLE_COMPLEX,k-1,0,MPI_COMM_WORLD,s_status,ierr) 
              if(myid.eq.k-1) CALL MPI_SEND(vals_inter(1:ind(k)-1),ind(k)-1,MPI_DOUBLE_COMPLEX,0,0,MPI_COMM_WORLD,ierr)
              if(myid.eq.0) CALL MPI_RECV(cols_inter(1:ind(k)-1),ind(k)-1,MPI_INTEGER,k-1,1,MPI_COMM_WORLD,s_status,ierr) 
              if(myid.eq.k-1) CALL MPI_SEND(cols_inter(1:ind(k)-1),ind(k)-1,MPI_INTEGER,0,1,MPI_COMM_WORLD,ierr)
              if(myid.eq.0) CALL MPI_RECV(rows_inter(1:ind(k)-1),ind(k)-1,MPI_INTEGER,k-1,3,MPI_COMM_WORLD,s_status,ierr) 
              if(myid.eq.k-1) CALL MPI_SEND(rows_inter(1:ind(k)-1),ind(k)-1,MPI_INTEGER,0,3,MPI_COMM_WORLD,ierr)


              if(myid.eq.0) then
               do i=1, ind(k)-1
                   Loo_vals_hilf(ind(1))=vals_inter(i)
                   Loo_cols_hilf(ind(1))=cols_inter(i)
                   Loo_rows_hilf(ind(1))=rows_inter(i)
                   ind(1)=ind(1)+1
               end do
               vals_inter=0d0
               cols_inter=0
               rows_inter=0 
              end if
             end if
            end do  



       call CPU_TIME(t2)
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)




           
       call CPU_TIME(t2)
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            j=j+nprocs
           end do


           if(j.gt.dimL_orb/2-nprocs+1.and.j.le.dimL_orb/2.and.myid.eq.0) then

              do while(j.le.dimL_orb/2) 
                L_orb_col=0d0 
                L_orb_col2=0d0 
                left_Loo_hilf=0d0 
                left_Loo_hilf2=0d0 

!              PUT KSL-mat entries in KSL_dist if preconstructed
                if(KSL_prec.eqv..TRUE.) then
                 pos2=mod(j,ND)
                 if(pos2.eq.0) pos2=ND
                 KSL_dist=0d0
                 KSL_conjg_dist=0d0 
                 do mm=1, Morb
                   do ll=1, Morb
                    k=KSL_mat_ind_csc(pos2,mm,ll)
                    if(k.ne.KSL_mat_ind_csc(pos2+1,mm,ll)) then
                    pp=pos2+1
                    do k=KSL_mat_ind_csc(pos2,mm,ll), KSL_mat_ind_csc(pp,mm,ll)-1 
                       KSL_dist(KSL_mat_rows_csc(k,mm,ll),mm,ll)= &
     &                           KSL_mat_vals_csc(k,mm,ll)
                    end do 
                    end if 
                    k=KSL_conjgmat_ind_csc(pos2,mm,ll)
                    if(k.ne.KSL_conjgmat_ind_csc(pos2+1,mm,ll)) then
                    pp=pos2+1
                    do k=KSL_conjgmat_ind_csc(pos2,mm,ll), &
     &                             KSL_conjgmat_ind_csc(pp,mm,ll)-1 
                       KSL_conjg_dist(KSL_conjgmat_rows_csc(k,mm,ll),mm,ll)= &
     &                           KSL_conjgmat_vals_csc(k,mm,ll)
                    end do 
                    end if
                   end do
                 end do 
                end if


!              PUT Tkin_2D entries in Tkin_col_p0 if preconstructed
               if(dim_mctdhb.eq.2.and.Tkin_prec.eqv..TRUE.) then
                 pos2=mod(j,ND)
                 if(pos2.eq.0) pos2=ND
                 Tkin_col_p0(:)=Tkin_2D(:,pos2)
               end if

                call construct_orb_col(PSI,L_orb_col,j)
                L_orb_col2(1:dimL_orb/2)=-Conjg(L_orb_col(dimL_orb/2+1:dimL_orb))
                L_orb_col2(dimL_orb/2+1:dimL_orb)=-Conjg(L_orb_col(1:dimL_orb/2))


                inter2_mat_small=0d0
                call mkl_zcoogemv('N', dimL_orb/2, proj_vals_coo, proj_rows_coo, &
     &                   proj_cols_coo, proj_nonzero, L_orb_col(1:dimL_orb/2), &
     &                   inter2_mat_small(:,1))

                call mkl_zcoogemv('N', dimL_orb/2, proj_vals_coo, proj_rows_coo, &
     &                   proj_cols_coo, proj_nonzero, L_orb_col2(1:dimL_orb/2), &
     &                   inter2_mat_small(:,2))



            do i=1, dimL_orb/2 
             pos=mod(i,ND)
             if(pos.eq.0) pos=ND
             do kk=1, Morb
              if(i.le.kk*ND) then
               orb_pos=kk
               exit
              end if
             end do
               do jj=1, Morb
                 left_Loo_hilf2(i,1)=left_Loo_hilf2(i,1)+ &
     &                         RhoInvSQR(orb_pos,jj)*inter2_mat_small(pos+(jj-1)*ND,1)
                 left_Loo_hilf2(i,2)=left_Loo_hilf2(i,2)+ &
     &                         RhoInvSQR(orb_pos,jj)*inter2_mat_small(pos+(jj-1)*ND,2)
               end do 
 

                if(myid.eq.0) then
                  if(abs(left_Loo_hilf2(i,1)).gt.tol1_Loo) then
                    Loo_vals_hilf(ind(1))=left_Loo_hilf2(i,1)
                    Loo_cols_hilf(ind(1))=j
                    Loo_rows_hilf(ind(1))=i
                    ind(1)=ind(1)+1
                  end if
                  if(abs(left_Loo_hilf2(i,2)).gt.tol1_Loo) then
                    Loo_vals_hilf(ind(1))=left_Loo_hilf2(i,2)
                    Loo_cols_hilf(ind(1))=j+dimL_orb/2
                    Loo_rows_hilf(ind(1))=i
                    ind(1)=ind(1)+1
                  end if
                end if

            end do


               j=j+1
              end do


           end if
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)


658       continue
        IF(MYID.eq.0) then

!       Conversion of Loo from sparse COO to CSR format
        mjob(1)=2
        mjob(2)=1
        mjob(3)=1
        mjob(4)=2
        mjob(5)=ind(1)
        mjob(6)=0


        counter=ind(1)-1 
        write(*,*) ind(1)-1

           if(.not.allocated(Loo_vals_csr)) & 
     &                allocate(Loo_vals_csr(counter))
           if(.not.allocated(Loo_col_csr)) &
     &                allocate(Loo_col_csr(counter))
           if(.not.allocated(Loo_ind_csr)) &
     &                allocate(Loo_ind_csr(dimL_orb/2+1))

         

       
        Loo_col_csr=0; Loo_vals_csr=0d0; Loo_ind_csr=0;
        call mkl_zcsrcoo(mjob, dimL_orb/2, Loo_vals_csr, Loo_col_csr, Loo_ind_csr, &
     &             counter, Loo_vals_hilf,  &
     &             Loo_rows_hilf, Loo_cols_hilf, info) 

        write(*,*) info, counter


            if(allocated(Loo_vals_hilf)) &
     &       deallocate(Loo_vals_hilf)
            if(allocated(Loo_rows_hilf)) &
     &       deallocate(Loo_rows_hilf)
            if(allocated(Loo_cols_hilf)) &
     &       deallocate(Loo_cols_hilf)
            if(allocated(inter_mat)) &
     &       deallocate(inter_mat)
            if(allocated(inter2_mat)) &
     &       deallocate(inter2_mat)
            if(allocated(inter_mat_small)) &
     &       deallocate(inter_mat_small)
            if(allocated(inter2_mat_small)) &
     &       deallocate(inter2_mat_small)
            if(allocated(inter_Loo)) &
     &       deallocate(inter_Loo)
            if(allocated(c_mat)) &
     &       deallocate(c_mat)


          close(256)

       END IF
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)





         if(allocated(L_orb_help2)) &
     &    deallocate(L_orb_help2)       
         if(allocated(KSL_hilf2_p0)) &
     &    deallocate(KSL_hilf2_p0)
         if(allocated(KSL_conjg_hilf2_p0)) &
     &    deallocate(KSL_conjg_hilf2_p0)
         if(allocated(left_Loo_hilf)) &
     &    deallocate(left_Loo_hilf)
         if(allocated(left_Loo_hilf_p0)) &
     &    deallocate(left_Loo_hilf_p0)
         if(allocated(left_Loo_hilf2)) &
     &    deallocate(left_Loo_hilf2)
          if(allocated(vals_inter)) deallocate(vals_inter) 
          if(allocated(rows_inter)) deallocate(rows_inter) 
          if(allocated(cols_inter)) deallocate(cols_inter) 

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

       END SUBROUTINE construct_left_of_Loo_upfront_sparse



!===CONSTRUCT ORB MATRIX COLUMN==============================================================
       SUBROUTINE construct_orb_col(PSI,L_orb_col,col)

          USE CI_ALL
          USE Parallel_CI 
          USE W_INTERPARTICLE 
          USE DVR_ALL
          USE rR_hW
          USE LR_RAPHA
          USE LR_ARNOLDI_MOD

          IMPLICIT NONE         
          include 'mpif.h'

          INTEGER :: i, j, k, q, s, ll, mm, kk, orb_pos, pos, counter, x_pos, y_pos, z_pos 
          INTEGER  :: i1,i2,j1,j2,q1,q2,myid,nprocs,ierr
          COMPLEX*16, DIMENSION(ND,Morb), INTENT(IN) :: PSI
          COMPLEX*16, DIMENSION(dimL_orb), INTENT(OUT) :: L_orb_col
          INTEGER, INTENT(IN)                          :: col          
          Real*8  :: start,finish 
          COMPLEX*16, DIMENSION(dimL_orb/2) :: h2_col
          COMPLEX*16, DIMENSION(dimL_orb/2) :: A, B, mu_col
          COMPLEX*16, DIMENSION(dimL_orb/2) :: Lu_col,Lv_col
          COMPLEX*16, DIMENSION(ND) :: Tkin_col
          COMPLEX*16, DIMENSION(ND) :: KSL_col,KSL1_col ! ADDED BY STR to KSL - non-local potential
          COMPLEX*16, DIMENSION(Morb,Morb) :: mu_PSI, RhoInvSQR, rmu, Rho_ij
          COMPLEX*16, DIMENSION(Morb,Morb,Morb,Morb) :: Rho_ijkl
          COMPLEX*16, DIMENSION(NDX*NDY*NDZ)  :: WSL,WSL1 !ADDED BY STR
          COMPLEX*16, DIMENSION(ND,ND)  :: ham


          call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
          call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )
 
           xlambda0=xlambda_0

!     Build matrix


         A=0.0d0; B=0.0d0; L_orb_col=0.0d0;
         mu_col=0.0d0; h2_col=0.0d0; Lu_col=0d0; Lv_col=0d0;


!      pos and orb_pos  necessary for assembling matrix 
         pos=mod(col,ND)
         if(pos.eq.0) pos=ND
!         write(*,*) "pos",pos

         do i=1, 2*Morb
              if(col.le.i*ND) then
               orb_pos=i
               exit
              end if
         end do
         if(orb_pos.gt.Morb) orb_pos=orb_pos-Morb 
!         write(*,*) "orb_pos",orb_pos


         counter=1
         Tkin_col=0d0


!     1D kinetic energy
      call CPU_TIME(start)
         if(dim_mctdhb.eq.1) then 
          IF(Time_DVRMETHODX==4) THEN ! FFT - explicit construction Tkin matrix
            IF(preconstr.eqv..FALSE.) then    
                call T_FFT_GENERAL(Tkin_col,pos)
            ELSE  
                Tkin_col(:)=Tkin_1D(:,pos) 
            END IF
          ELSE  ! HO_DVR CASE
             Tkin_col(:)=Op_x(:,pos)
          END IF
         end if

        
      call CPU_TIME(start)
!        2D kinetic energy
         if(dim_mctdhb.eq.2) then 

          IF(Tkin_prec.eqv..FALSE.) then      !FOR REALLY LARGE PROBLEMS
                                      !CALCULATE column every time again

           IF(Time_DVRMETHODX.eq.4.and.Time_DVRMETHODY.eq.4) THEN ! FFT - explicit construction Tkin matrix
                call T_FFT_GENERAL(Tkin_col,pos)
           ELSE 
            ham=0d0
            if(mod(col,NDX).eq.0) then
              y_pos=col/NDX
            else
              y_pos=(col/NDX)+1 
            end if

            x_pos=mod(col,NDX)
            if(x_pos.eq.0) x_pos=NDX 

            DO j1=1,NDY
             j2=y_pos 
             if(j1.eq.j2) then
              q1=1
             else
              q1=0
             endif
             Do i1=1,NDX
              i2=x_pos
              if(i1.eq.i2) then
               q2=1
              else
               q2=0
              endif
              I=NDX*(j1-1)+i1
              Tkin_col(I)=q1*Op_X(i1,i2)+q2*Op_Y(j1,j2)
             endDo
            ENDDO
           END IF 

          ELSE         !SMALLER PROBLEMS WHERE Tkin_2D already available
!            if(myid.eq.0) Tkin_col(:)=Tkin_2D(:,pos) 
             Tkin_col(:)=Tkin_col_p0(:) 
          ENDIF 
         end if
      call CPU_TIME(finish)
!      write(*,*) "Tkin",finish-start 


        if(dim_mctdhb.eq.3) then 
         do i=1, NDX
          do j=1, NDY
           do k=1, NDZ
            Tkin_col(counter)=Tkin_col(counter)+Op_x(pos,i)+Op_y(pos,j)+Op_z(pos,k)
!            if (counter.eq.pos) Tkin_col(counter)=Tkin_col(counter)+VTRAP_EXT(counter) 
            counter=counter+1 
           end do
          end do
         end do 
        end if


!      Add external potential on diagonal
       Tkin_col(pos)=Tkin_col(pos)+VTRAP_EXT(pos)
      call CPU_TIME(finish)
!      write(*,*) "Tkin",finish-start 
 
!       END IF    ! end if FFT or not


!       Get Rho_ij
         call Get_Full_Rij(Rho_ij)

!       Get Rho_ijikl
         call Get_Full_Rijkl(Rho_ijkl)

!       Assemble h2_col
         do q=1, Morb
          do k=1, ND
           h2_col((q-1)*ND+k)=Rho_ij(q,orb_pos)*Tkin_col(k)
          end do
         end do


!     get Lagrange multipliers
              ZmuR= MAtmul(Zmu,Transpose(-Rho_ij)) !AIS 18JUL14 Zmu is already available after execution  Func.F
 
               DO q=1,Morb
                     mu_col((q-1)*ND+pos)= &
     &                   ZmuR(q,orb_pos)
               END DO

!     ADD LOCAL AND EXCHANGE POTENTIALS
 
      if(xlambda0.ne.0d0) then
      call CPU_TIME(start)
!      write(*,*) "preconstr",preconstr 
            DO k=1,Morb
!               DO q=1,Morb
           
                     counter=0 
                     DO s=1,Morb
                        DO ll=1,Morb

                          counter=counter+1
                          if(preconstr.eqv..FALSE.) then
                           CALL  Get_WSL_omp(WSL,psi(:,s),psi(:,ll))  
                           CALL  Get_KSL_col(KSL_col,psi(:,s),psi(:,ll),col) 
                           if(orb_real) then
                              KSL1_col(:)=KSL_col(:)
                           else
                             CALL  Get_KSL_col(KSL1_col,CONJG(psi(:,s)),psi(:,ll),col) 
                           end if
                          else
                           WSL=WSL_mat(:,s,ll)
                           if(KSL_prec.eqv..TRUE.) then 
                             KSL_col(:)=KSL_dist(:,s,ll)
                             KSL1_col(:)=KSL_conjg_dist(:,s,ll)
!ORIG                             KSL_col(:)=KSL_mat(:,pos,s,ll)
!ORIG                             KSL1_col(:)=KSL_conjg_mat(:,pos,s,ll)
                           else
                             CALL  Get_KSL_col(KSL_col,psi(:,s),psi(:,ll),col) 
                             if(orb_real) then
                               KSL1_col(:)=KSL_col(:)
                             else 
                               CALL  Get_KSL_col(KSL1_col,CONJG(psi(:,s)),psi(:,ll),col) 
                             end if
                           end if
!                           if(orb_real) then
!                             KSL1_col(:)=KSL_col(:)
!                           else 
!                             CALL  Get_KSL_col(KSL1_col,CONJG(psi(:,s)),psi(:,ll),col) 
!                           end if
                          end if 

                    A((k-1)*ND+pos)= A((k-1)*ND+pos) +Rho_ijkl(k,s,ll,orb_pos)*WSL(pos) 
                  DO i=1,ND
!                  DO j=1,ND
 
                           A((k-1)*ND+i)= &
     &                        A((k-1)*ND+i) &
     &                        +Rho_ijkl(k,s,ll,orb_pos)*KSL_col(i) 


                           B((k-1)*ND+i)= &
     &                        B((k-1)*ND+i) &
     &                        +Rho_ijkl(k,orb_pos,s,ll)*KSL1_col(i)  
!                  END DO
                  END DO
                        
                        END DO
                     END DO

!               END DO
            END DO  
      end if
      call CPU_TIME(finish)
!      write(*,*) "INt",finish-start 
    

!     Assemble Lu and Lv coulumns
            Lu_col=h2_col-mu_col+A
            Lv_col=B

!     Assemble full coulumn
       if(col.le.Morb*ND) then
           L_orb_col(1:dimL_orb/2)=Lu_col(:)
           L_orb_col(dimL_orb/2+1:dimL_orb)=-Conjg(Lv_col(:))
       else
           L_orb_col(1:dimL_orb/2)=Lv_col(:)
           L_orb_col(dimL_orb/2+1:dimL_orb)=-Conjg(Lu_col(:))
       end if
         

      END SUBROUTINE construct_orb_col

!===========================================================================================

!      Column of projector matrix
       SUBROUTINE get_proj_col(PSI,Pb,col,pos)
 
           USE CI_ALL
           USE Parallel_CI
           USE DVR_ALL
 
           IMPLICIT NONE

           INTEGER :: i, j, k, q, s, ll, mm, col, pos, orb_pos

           COMPLEX*16, DIMENSION(ND,Morb), INTENT(IN) :: PSI
           COMPLEX*16, DIMENSION(dimL_orb), INTENT(OUT) :: Pb
           COMPLEX*16, DIMENSION(ND) :: Pbb 
           COMPLEX*16 :: zdotc 

           Pb=0.0d0; Pbb=0d0;

!       Get orb_pos
         do k=1, 2*Morb
              if(col.le.k*ND) then
               orb_pos=k
               exit
              end if
         end do
         if(orb_pos.gt.Morb) orb_pos=orb_pos-Morb 

!         do i=1, 2*Morb
!           if(col.ge.(i-1)*ND.and.col.le.i*ND) then
!             orb_pos=i
!             goto 456
!           end if
!         end do
!456      continue 

!         write(*,*) "orb_pos line 545",orb_pos,pos
!         write(*,*) PSI

 
!        Get column elements
          DO i=1,ND
                   Pbb(i)=-SUM(PSI(i,:) &
     &                 *CONJG(PSI(pos,:)))
          END DO


       
            if(col.gt.Morb*ND) then
              Pb(dimL_orb/2+(orb_pos-1)*ND+1:dimL_orb/2+orb_pos*ND)=Conjg(Pbb(:)) 
            else
              Pb(1+(orb_pos-1)*ND:orb_pos*ND)=Pbb(:)
            end if 
            Pb(col)=Pb(col)+1.0d0    ! ADD UNITY ON DIAGONAL ELEMENT         
          
 
       END SUBROUTINE get_proj_col

!=====================================================================================


      SUBROUTINE get_mu_test(mu_PSI)

          USE CI_ALL
          USE Parallel_CI 
          USE DVR_ALL

          IMPLICIT NONE         

          INTEGER :: i, j, k, q, s, ll, mm 

          COMPLEX*16, DIMENSION(Morb,Morb), INTENT(OUT) :: mu_PSI
          COMPLEX*16, DIMENSION(Morb,Morb) :: rmu
          COMPLEX*16 :: en_PSI
          COMPLEX*16, DIMENSION(Morb,Morb) :: Rho_ij
          COMPLEX*16, DIMENSION(Morb,Morb) :: H_ij
          COMPLEX*16, DIMENSION(Morb,Morb) :: RhoINVSQR
          COMPLEX*16, DIMENSION(Morb,Morb,Morb,Morb) :: Rho_ijkl
          COMPLEX*16, DIMENSION(Morb,Morb,Morb,Morb) :: W_ijkl

          mu_PSI=0.0d0; rmu=0.0d0; en_PSI=0.0d0

          call Get_Full_Rij(Rho_ij)
          call Get_Full_Rijkl(Rho_ijkl)
          call Get_Full_Hij(H_ij)
          call Get_Full_Wijkl(W_ijkl)

          DO k=1,Morb
             DO j=1,Morb

                DO q=1,Morb

                   mu_PSI(k,j)=mu_PSI(k,j) + Rho_ij(k,q)*H_IJ(j,q) 
                END DO


                DO q=1,Morb
                   DO s=1,Morb
                      DO ll=1,Morb

                         mu_PSI(k,j)=mu_PSI(k,j) + Rho_ijkl(k,s,q,ll)*W_ijkl(j,s,q,ll)       
                      END DO
                   END DO
                END DO


             END DO
          END DO

          DO k=1,Morb
             DO j=1,Morb

                   en_PSI=en_PSI + Rho_ij(k,j)*H_IJ(k,j) 

                DO q=1,Morb
                   DO s=1,Morb

                         en_PSI=en_PSI + 0.5d0*Rho_ijkl(k,j,q,s)*W_ijkl(k,j,q,s)       
                   END DO
                END DO

             END DO
          END DO

       !call  squarerootRInv_arnoldi(RhoInvSQR)
       !mu_PSI= MAtmul(mu_PSI,Transpose(RhoInvSQR))
          

         ! Test rhoinv * mu
          DO k=1,Morb
             DO j=1,Morb
                   rmu(k,j)=mu_PSI(k,j)
             END DO
          END DO


      END SUBROUTINE get_mu_test
!=====================================================================================

!      Inverse squareroot one-body density matrix, same as before!!!!
       SUBROUTINE squarerootRInv_arnoldi(RhoInvSQR)
 
           USE CI_ALL
           USE Parallel_CI
           USE DVR_ALL
           USE rR_hW
  

           IMPLICIT NONE
 
           COMPLEX*16, DIMENSION(Morb,Morb), INTENT(OUT) :: RhoInvSQR
 
           COMPLEX*16, DIMENSION(Morb,Morb) :: V, VI, Vhelp
           Real*8, DIMENSION(Morb) :: w
 
           INTEGER :: i, j, k, info, infoLU, infoInv
           INTEGER, PARAMETER :: lwork=1000, lrwork=1000, liwork=1000, lworkInv=1000
           COMPLEX*16, DIMENSION(lwork) :: work
           REAL*8, DIMENSION(lrwork) :: rwork
           INTEGER, DIMENSION(liwork) :: iwork
           COMPLEX*16, DIMENSION(lworkInv) :: workInv
           INTEGER, DIMENSION(Morb) :: ipiv
 
 !     Eigensystem of InvZRIJ: w, V 
 
 
           V=InvZRIJ
!           write(*,*) myid,"V Arn",V
           CALL ZHEEVD('V','U', Morb, V, Morb, w, work, lwork, rwork, lrwork, iwork, &
      &       liwork, info)
 
           if(info.ne.0) WRITE(6,*) 'info ZHEEVD in squarerootRInv:', info, 'work(1): ', &
      &              work(1), 'lwork: ', lwork, 'rwork(1): ', rwork(1), &
      &             'lrwork: ', lrwork, 'iwork(1): ', iwork(1), &
      &             'liwork: ', liwork
 
 
!     Inverse of V
!     LU
           VI=V
 
           CALL ZGETRF( Morb, Morb, VI, Morb, ipiv, infoLU )
 
           if(info.ne.0) WRITE(6,*) 'info ZGETRF in squarerootRInv:', infoLU
 
!     INV
           CALL ZGETRI( Morb, VI, Morb, ipiv, workInv, lworkInv, infoInv )
 
            if(info.ne.0) WRITE(6,*) 'info ZGETRI in squarerootRInv:', infoInv, 'work(1):', workInv(1)
 

!     Calculate squareroot of inverse matrix
 
 
           RhoInvSQR=0.0d0
           DO i=1,Morb
              DO j=1,Morb
                 Vhelp(i,j)=V(i,j)*SQRT(w(j))
              END DO
           END DO
 
          CALL ZGEMM('N', 'N', Morb, Morb, Morb, DCMPLX(1.0d0),Vhelp, &
      &                      Morb, VI, Morb, DCMPLX(0.0d0), RhoInvSQR, Morb)
 
 
        END SUBROUTINE squarerootRInv_arnoldi

!===========================================================================
!      Squareroot one-body density matrix, same as before!!!!
       SUBROUTINE squarerootR_arnoldi(RhoSQR)
 
           USE CI_ALL
           USE Parallel_CI
           USE DVR_ALL
           USE rR_hW
  

           IMPLICIT NONE
 
           COMPLEX*16, DIMENSION(Morb,Morb), INTENT(OUT) :: RhoSQR
 
           COMPLEX*16, DIMENSION(Morb,Morb) :: V, VI, Vhelp
           Real*8, DIMENSION(Morb) :: w
 
           INTEGER :: i, j, k, info, infoLU, infoInv
           INTEGER, PARAMETER :: lwork=1000, lrwork=1000, liwork=1000, lworkInv=1000
           COMPLEX*16, DIMENSION(lwork) :: work
           REAL*8, DIMENSION(lrwork) :: rwork
           INTEGER, DIMENSION(liwork) :: iwork
           COMPLEX*16, DIMENSION(lworkInv) :: workInv
           INTEGER, DIMENSION(Morb) :: ipiv
 
 !     Eigensystem of InvZRIJ: w, V 

        call Get_Full_Rij(RhoSQR)  
 
           V=RhoSQR
!           write(*,*) myid,"V Arn",V
           CALL ZHEEVD('V','U', Morb, V, Morb, w, work, lwork, rwork, lrwork, iwork, &
      &       liwork, info)
 
           if(info.ne.0) WRITE(6,*) 'info ZHEEVD in squarerootRInv:', info, 'work(1): ', &
      &              work(1), 'lwork: ', lwork, 'rwork(1): ', rwork(1), &
      &             'lrwork: ', lrwork, 'iwork(1): ', iwork(1), &
      &             'liwork: ', liwork
 
 
!     Inverse of V
!     LU
           VI=V
 
           CALL ZGETRF( Morb, Morb, VI, Morb, ipiv, infoLU )
 
           if(info.ne.0) WRITE(6,*) 'info ZGETRF in squarerootRInv:', infoLU
 
!     INV
           CALL ZGETRI( Morb, VI, Morb, ipiv, workInv, lworkInv, infoInv )
 
            if(info.ne.0) WRITE(6,*) 'info ZGETRI in squarerootRInv:', infoInv, 'work(1):', workInv(1)
 

!     Calculate squareroot of inverse matrix
 
 
           RhoSQR=0.0d0
           DO i=1,Morb
              DO j=1,Morb
                 Vhelp(i,j)=V(i,j)*SQRT(w(j))
              END DO
           END DO
 
          CALL ZGEMM('N', 'N', Morb, Morb, Morb, DCMPLX(1.0d0),Vhelp, &
      &                      Morb, VI, Morb, DCMPLX(0.0d0), RhoSQR, Morb)
 
 
        END SUBROUTINE squarerootR_arnoldi


!=====================================================================================================

       SUBROUTINE construct_CO_matrix_col_bare(PSI,VIN,COmat_bare_col,Nc,col)

          USE CI_ALL
          USE Parallel_CI 
          USE DVR_ALL
       USE   W_INTERPARTICLE !ADD BY STR
          IMPLICIT NONE
          COMPLEX*16, INTENT(IN) :: VIN(Nc)
          COMPLEX*16  :: COmat_col(2*NC),COmat_col2(2*NC)
          COMPLEX*16, INTENT(OUT) :: COmat_bare_col(2*NC)
          COMPLEX*16  :: PL(2*ND*Morb)
          COMPLEX*16, DIMENSION(ND,Morb), INTENT(IN) :: PSI
          COMPLEX*16 :: Vhelp1(Nc), Vhelp2(Nc)
          COMPLEX*16 :: Vhelp1a(Nc), Vhelp2a(Nc)
          COMPLEX*16 :: Vhelp1b(Nc), Vhelp2b(Nc)
!          COMPLEX*16, DIMENSION(ND,Morb) :: h2_PSI
!          COMPLEX*16 :: COmat_help(2*NC)
          COMPLEX*16, DIMENSION(Morb,Morb) :: Rho_arnoldi
          INTEGER, DIMENSION(2) :: ind11, ind12, ind21, ind22
          INTEGER, DIMENSION(2) :: ind11a, ind12a, ind21a, ind22a
          COMPLEX*16 :: VIN_help(Nc)
          COMPLEX*16 :: rhohelp
         COMPLEX*16, DIMENSION(NDX*NDY*NDZ)  :: WSL !ADDED BY STR
!==========================================================
          INTEGER :: i,j,k,s,l,q,Nc,col,pos,CI_pos,orb_pos
!===========================================================
           xlambda0=xlambda_0

!         Initialize

          COmat_col2=0d0 
          COmat_col=0.0d0; Vhelp1a=0.0d0; Vhelp2a=0.0d0; Vhelp1b=0.0d0; Vhelp2b=0.0d0
          CALL get_h2_PSI(PSI,h2_PSI)

          pos=mod(col,ND)
          if(pos.eq.0) pos=ND  

         do i=1, 2*Morb
              if(col.le.i*ND) then
               orb_pos=i
               exit
              end if
         end do
         if(orb_pos.gt.Morb) orb_pos=orb_pos-Morb 




!         UPPER LEFT         
!RB          DO q=1,Morb
             DO k=1,Morb

!               1-body
                VIN_help=VIN; Vhelp1a=0.0d0
                CALL Produce_Cij(VIN_help,Vhelp1a,k,orb_pos,0,0,Nc) !Correct - proved by STR
!                  write(*,*) "Vhelp1a(row)",row,Vhelp1a(CI_pos)
!RB FEB29                rhohelp=SUM(Conjg(VIN_help)*Vhelp1a)


                DO i=1,Nc
!RB                   DO j=1,ND
                      COmat_col(i)=COmat_col(i)+Vhelp1a(i)*CONJG(h2_PSI(pos,k))
!RB                   END DO
                END DO



!               2-body

                DO s=1,Morb
                   DO l=1,Morb
!                      IF(k==l.AND.s==k) THEN
!                      IF(q==s.AND.k==l.AND.q==k) GOTO 111
                      VIN_help=VIN; Vhelp2a=0.0d0
                      CALL Produce_Cij(VIN_help,Vhelp2a,k,s,orb_pos,l,Nc) !Correct - proved by STR
!RB FEB29                rhohelp=SUM(Conjg(VIN)*Vhelp2a)
!                  write(*,*) "Vhelp2a(row)",row,Vhelp2a(CI_pos)
                 

            if(ND.gt.128*128) then
               CALL  Get_WSL_omp(WSL,psi(:,s),psi(:,l))  
            else
               WSL=WSL_mat(:,s,l)
            end if 
!            CALL  Get_WSL_omp(WSL,psi(:,s),psi(:,l))

                      DO i=1,Nc
!RB                         DO j=1,ND
             
!                  COmat(i,(q-1)*ND+j)=COmat(i,(q-1)*ND+j)+xlambda_0*Vhelp2a(i)*CONJG(PSI(j,k))*CONJG(PSI(j,s))*PSI(j,l)
                            COmat_col(i)=COmat_col(i)+Vhelp2a(i)*CONJG(PSI(pos,k))*WSL(pos)
                             
!RB                         END DO
                      END DO

                   END DO
                END DO

             END DO
!RB          END DO



!         upper right
          CALL get_h2_PSI(CONJG(PSI),h2_PSI)

!RB          DO q=1,Morb
             DO k=1,Morb


!               1-body
                VIN_help=VIN; Vhelp1b=0.0d0
                CALL Produce_Cij(VIN_help,Vhelp1b,orb_pos,k,0,0,Nc) !Correct - proved by STR
!                  write(*,*) "Vhelp1b(row)",row,Vhelp1b(CI_pos)
!RB FEB29                rhohelp=SUM(Conjg(VIN)*Vhelp1b)


                DO i=1,Nc
!RB                   DO j=1,ND
             
                      COmat_col2(i)=COmat_col2(i)+Vhelp1b(i)*CONJG(h2_PSI(pos,k))
!RB                   END DO
                END DO

!               2-body

                DO s=1,Morb
                   DO l=1,Morb
                      VIN_help=VIN; Vhelp2b=0.0d0
!                  write(*,*) "Vhelp2b(row)",row,Vhelp2b(CI_pos)
                      CALL Produce_Cij(VIN_help,Vhelp2b,orb_pos,s,k,l,Nc)
!RB FEB29                      rhohelp=SUM(Conjg(VIN)*Vhelp2b)
            
            if(ND.gt.128*128) then
               CALL  Get_WSL_omp(WSL,psi(:,s),psi(:,l))  
            else
               WSL=WSL_mat(:,s,l)
            end if 
!            CALL  Get_WSL_omp(WSL,psi(:,s),psi(:,l))

                      DO i=1,Nc
!RB                         DO j=1,ND
             
                        COmat_col2(i)=COmat_col2(i)+Vhelp2b(i)*PSI(pos,k)*WSL(pos)

!RB                         END DO
                      END DO

                   END DO
                END DO

             END DO


!RB          END DO

          if(col.le.dimL_orb/2) then
            COmat_bare_col(1:NC)=COmat_col
            COmat_bare_col(NC+1:2*NC)=-Conjg(COmat_col2)
          else
            COmat_bare_col(1:NC)=COmat_col2
            COmat_bare_col(NC+1:2*NC)=-Conjg(COmat_col)
          end if
          COmat_col=COmat_bare_col 


!           open(unit=256, file='test.dat', status='replace', action='readwrite')
!           do q=1, Morb
!           do i=1, 2*NC
!            write(256,'(500(2F20.10))') (Real(Tkin(i,j)),Dimag(Tkin(i,j)), &
!     &        j=1,nd,1)
!             write(256,'(500(2F20.10))') Real(COmat_col(i)),Dimag(COmat_col(i))
!           end do
!           end do 
!          close(256)


       END SUBROUTINE construct_CO_matrix_col_bare

!==================================================================================================



       SUBROUTINE construct_OC_matrix_row(PSI,VIN,OCmat_rows,Nc,row)

          USE CI_ALL
          USE Parallel_CI 
          IMPLICIT NONE
          COMPLEX*16, INTENT(IN) :: VIN(Nc)
          COMPLEX*16, INTENT(OUT) :: OCmat_rows(Morb,2*Nc)
!          COMPLEX*16, DIMENSION(ND,Morb), INTENT(IN) :: h2_PSI,h2_PSI_conjg
!          COMPLEX*16, INTENT(OUT) :: OCmat_row(2*Nc)
          COMPLEX*16, DIMENSION(ND,Morb), INTENT(IN) :: PSI
          COMPLEX*16, DIMENSION(Morb,2*Morb*ND) :: left
          COMPLEX*16, DIMENSION(Morb,2*NC) :: inter
          COMPLEX*16, DIMENSION(2*NC) :: final
          COMPLEX*16, DIMENSION(Morb,2*NC) :: final2
          COMPLEX*16 :: COmat_bare_col(2*NC)
          COMPLEX*16, DIMENSION(dimL_orb) :: Pb,COmat_row,COmat_bare_row,OCmat_bare_col1, OCmat_bare_col2
          COMPLEX*16, DIMENSION(Morb,Morb) :: Rho_arnoldi
          INTEGER :: i,j,k,s,l,q,Nc,pos,orb_pos,row
          LOGICAL ::  only_bare
          Real*8 :: start,finish 
!===========================================================


          OCmat_rows=0.0d0; OCmat_bare_col1=0d0; OCmat_bare_col2=0d0;


          pos=mod(row,ND)
          if(pos.eq.0) pos=ND
!          write(*,*) "pos",pos

         if(left_Loc_prec.eqv..FALSE.) then
         

!        Get inverse squareroot matrix
         call squarerootRInv_arnoldi(Rho_arnoldi)




         do i=1, 2*Morb
              if(row.le.i*ND) then
               orb_pos=i
               exit
              end if
         end do
         if(orb_pos.gt.Morb) orb_pos=orb_pos-Morb 



!       Get needed projector rows (Morb rows, stored in array left(:,:))
           do i=1, Morb
            Pb=0d0
            if(row.le.ND*Morb) then
              call get_proj_col(PSI,Pb,(i-1)*ND+pos,pos)
              left(i,:)=Conjg(Pb)
            else
              call get_proj_col(PSI,Pb,Morb*ND+(i-1)*ND+pos,pos)
              left(i,:)=Conjg(Pb)
            end if 
           end do  

           
!       Multiply the Morb projectors rows with columns of bare CO_matrix
!                call CPU_TIME(start) 
           inter=0d0 
           only_bare=.TRUE.
           do j=1, NC
             COmat_bare_row=0d0 
!             call construct_CO_matrix_col_bare(PSI,VIN,COmat_bare_col,Nc,j)
             call construct_CO_matrix_row(PSI,VIN,only_bare,COmat_row,COmat_bare_row,Nc,j)

             OCmat_bare_col1(1:dimL_orb/2)=Conjg(COmat_bare_row(1:dimL_orb/2))
             OCmat_bare_col1(dimL_orb/2+1:dimL_orb)=-Conjg(COmat_bare_row(dimL_orb/2+1:dimL_orb))


             OCmat_bare_col2(1:dimL_orb/2)=COmat_bare_row(dimL_orb/2+1:dimL_orb)
             OCmat_bare_col2(dimL_orb/2+1:dimL_orb)=-COmat_bare_row(1:dimL_orb/2)


             do i=1, Morb
!              do k=1, dimL_orb/2
!                inter(i,j)=inter(i,j)+Conjg(left(i,k))*OCmat_bare_col1(k)
!                inter(i,j+NC)=inter(i,j+NC)+Conjg(left(i,k))*OCmat_bare_col2(k)
                inter(i,j)=dot_product(Conjg(left(i,:)),OCmat_bare_col1)
                inter(i,j+NC)=dot_product(Conjg(left(i,:)),OCmat_bare_col2)
!              end do
             end do  

           end do
!                call CPU_TIME(finish)
!                write(*,*) finish-start




!        Rho from the left
          final2=0d0
          if(row.le.Morb*ND) then
            do s=1, Morb 
             do k=1, Morb
              do j=1, NC
                final2(s,j)=final2(s,j)+Rho_arnoldi(s,k)*inter(k,j)
                final2(s,NC+j)=final2(s,NC+j)+Rho_arnoldi(s,k)*inter(k,j+NC)
              end do 
             end do 
            end do
          else
            do s=1, Morb 
             do k=1, Morb
              do j=1, NC
                final2(s,j)=final2(s,j)+Conjg(Rho_arnoldi(s,k))*inter(k,j)
                final2(s,NC+j)=final2(s,NC+j)+Conjg(Rho_arnoldi(s,k))*inter(k,j+NC)
              end do 
             end do 
            end do
          end if 



          do i=1, Morb
            OCmat_rows(i,:)=final2(i,:)
          end do 

!         OCmat_row=final



         else

          do i=1, Morb
            if(row.le.dimL_orb/2) then
!ORIG              OCmat_rows(i,:)=left_Loc((i-1)*ND+pos,:)
              OCmat_rows(i,:)=left_Loc_hilf2(i,:) 
            else
!ORIG              OCmat_rows(i,:)=left_Loc((Morb+i-1)*ND+pos,:)
              OCmat_rows(i,:)=left_Loc_hilf2(i,:) 
            end if
          end do 


         end if    ! end if left_Loc_prec  



!           open(unit=256, file='test.dat', status='replace', action='readwrite')
!           do i=1, 2*NC
!            write(256,'(500(2F20.10))') (Real(Tkin(i,j)),Dimag(Tkin(i,j)), &
!     &        j=1,nd,1)
!             write(256,'(500(2F20.10))') Real(OCmat_row(i)),Dimag(OCmat_row(i))
!           end do
!          close(256)


       END SUBROUTINE construct_OC_matrix_row



!==================================================================================================

!=====  ONLY USE this ROUTINE IF MATRIX SMALL ENOUGH...
!=====  AND IF M >= 4 or 5 ...
       SUBROUTINE construct_Loc_matrix_full(PSI,VIN,NC)

          USE CI_ALL
          USE Parallel_CI 
          IMPLICIT NONE
          include 'mpif.h'

          COMPLEX*16, Dimension(Nc), INTENT(IN) :: VIN
!          COMPLEX*16, DIMENSION(ND,Morb), INTENT(IN) :: h2_PSI,h2_PSI_conjg
!          COMPLEX*16, INTENT(OUT) :: OCmat_row(2*Nc)
          COMPLEX*16, DIMENSION(ND,Morb), INTENT(IN) :: PSI
          COMPLEX*16 :: COmat_bare_col(2*NC)
          COMPLEX*16, DIMENSION(dimL_orb) :: Pb,COmat_row,COmat_bare_row,OCmat_bare_col1, OCmat_bare_col2
          COMPLEX*16, DIMENSION(Morb,Morb) :: RhoInvSQR
          INTEGER :: i,j,k,s,l,q,Nc,pos,orb_pos,row,start_from, myid, nprocs, ierr, jmax
          LOGICAL ::  only_bare, exists
          Real*8 :: start,finish,a,b 
          Complex*16, allocatable :: OCmat_help(:,:)
          INTEGER, DIMENSION(2) :: ind11, ind12, ind21, ind22
          INTEGER, DIMENSION(2) :: ind11a, ind12a, ind21a, ind22a
          COMPLEX*16, DIMENSION(dimL_orb/2,2) :: left_Loc_hilf,left_Loc_hilf_p0
          INTEGER :: s_status(MPI_STATUS_SIZE), st


          call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
          call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )

          if(myid.eq.0) then 

           INQUIRE(FILE="left_Loc_inter.dat", EXIST=exists)
!          CALL MPI_BCAST(exists,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)

!          left_Loc=0d0


          if(exists.eqv..FALSE.) then  !files do not exist       
 
            open(unit=256, file='left_Loc_inter.dat',status='replace', action='readwrite')

            start_from=1

          else

            open(unit=256, file='left_Loc_inter.dat',status='old', action='readwrite')

            jmax=1   
            do 
              read(256,*,iostat=st) i,j,a,b
              if(st.ne.0) exit
              left_Loc(i,j)=Cmplx(a,b,kind=8)
              if(j.le.NC) then 
                if(j.gt.jmax) jmax=j  
              end if
            end do
            start_from=jmax  

          end if
          end if

          CALL MPI_BCAST(exists,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(start_from,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

          if(start_from.eq.NC) goto 698

          j=start_from

!         Calculate Proj*L_oc
           only_bare=.TRUE.
           do while(j.le.NC-nprocs+1) 
!           do j=start_from, NC
             COmat_bare_row=0d0 
             COmat_row=0d0 
!             call construct_CO_matrix_col_bare(PSI,VIN,COmat_bare_col,Nc,j)
             call construct_CO_matrix_row(PSI,VIN,only_bare,COmat_row,COmat_bare_row,Nc,j+myid)

             OCmat_bare_col1(1:dimL_orb/2)=Conjg(COmat_bare_row(1:dimL_orb/2))
             OCmat_bare_col1(dimL_orb/2+1:dimL_orb)=-Conjg(COmat_bare_row(dimL_orb/2+1:dimL_orb))


             OCmat_bare_col2(1:dimL_orb/2)=COmat_bare_row(dimL_orb/2+1:dimL_orb)
             OCmat_bare_col2(dimL_orb/2+1:dimL_orb)=-COmat_bare_row(1:dimL_orb/2)




            do i=1, dimL_orb/2 
             Pb=0d0  
             pos=mod(i,ND)
             if(pos.eq.0) pos=ND
             call get_proj_col(PSI,Pb,i,pos)
             Pb=Conjg(Pb)

!             left_Loc(i,j+myid)=dot_product(Conjg(Pb),OCmat_bare_col1) 
!             left_Loc(i,j+myid+NC)=dot_product(Conjg(Pb),OCmat_bare_col2) 

             left_Loc_hilf(i,1)=dot_product(Conjg(Pb),OCmat_bare_col1) 
             left_Loc_hilf(i,2)=dot_product(Conjg(Pb),OCmat_bare_col2) 
              
            end do 

             if(myid.eq.0) then
               left_Loc(:,j)=left_Loc_hilf(:,1)
               left_Loc(:,j+NC)=left_Loc_hilf(:,2)

               do i=1, dimL_orb/2
                 if(abs(left_Loc(i,j)).gt.tol1) then
                   write(256,'(2I10,10x,2F20.10)') i,j,Real(left_Loc(i,j)),Dimag(left_Loc(i,j))
                 end if
                 if(abs(left_Loc(i,j+NC)).gt.tol1) then
                   write(256,'(2I10,10x,2F20.10)') i,j+NC,Real(left_Loc(i,j+NC)),Dimag(left_Loc(i,j+NC))
                 end if
               end do
             end if 

            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            do k=2, nprocs

              if(myid.eq.k-1) then
                 left_Loc_hilf_p0(:,1)=left_Loc_hilf(:,1) 
                 left_Loc_hilf_p0(:,2)=left_Loc_hilf(:,2) 
              end if


              if(myid.eq.0) CALL MPI_RECV(left_Loc_hilf_p0(:,1),dimL_orb/2,MPI_DOUBLE_COMPLEX,k-1,0,MPI_COMM_WORLD,s_status,ierr) 
              if(myid.eq.k-1) CALL MPI_SEND(left_Loc_hilf_p0(:,1),dimL_orb/2,MPI_DOUBLE_COMPLEX,0,0,MPI_COMM_WORLD,ierr)
              if(myid.eq.0) CALL MPI_RECV(left_Loc_hilf_p0(:,2),dimL_orb/2,MPI_DOUBLE_COMPLEX,k-1,1,MPI_COMM_WORLD,s_status,ierr) 
              if(myid.eq.k-1) CALL MPI_SEND(left_Loc_hilf_p0(:,2),dimL_orb/2,MPI_DOUBLE_COMPLEX,0,1,MPI_COMM_WORLD,ierr)

!              CALL MPI_BCAST(left_Loc_hilf_p0(:,1),dimL_orb/2,MPI_DOUBLE_COMPLEX,k-1,MPI_COMM_WORLD,ierr)
!              CALL MPI_BCAST(left_Loc_hilf_p0(:,2),dimL_orb/2,MPI_DOUBLE_COMPLEX,k-1,MPI_COMM_WORLD,ierr)

!ORIG              CALL MPI_BCAST(left_Loc(:,j+k-1),dimL_orb/2,MPI_DOUBLE_COMPLEX,k-1,MPI_COMM_WORLD,ierr)
!ORIG              CALL MPI_BCAST(left_loc(:,j+k-1+NC),dimL_orb/2,MPI_DOUBLE_COMPLEX,k-1,MPI_COMM_WORLD,ierr)
 
             if(myid.eq.0) then
               left_Loc(:,j+k-1)=left_Loc_hilf_p0(:,1)
               left_Loc(:,j+k-1+NC)=left_Loc_hilf_p0(:,2)

               do i=1, dimL_orb/2
                 if(abs(left_Loc(i,j+k-1)).gt.tol1) then
                   write(256,'(2I10,10x,2F20.10)') i,j+k-1,Real(left_Loc(i,j+k-1)),Dimag(left_Loc(i,j+k-1))
                 end if
                 if(abs(left_Loc(i,j+k-1+NC)).gt.tol1) then
                   write(256,'(2I10,10x,2F20.10)') i,j+k-1+NC,Real(left_Loc(i,j+k-1+NC)),Dimag(left_Loc(i,j+k-1+NC))
                 end if
               end do
             end if 
             call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            end do  


            j=j+nprocs
           end do


           if(j.gt.NC-nprocs+1.and.j.le.NC.and.myid.eq.0) then

           do while(j.le.NC) 
!           do j=start_from, NC
             COmat_bare_row=0d0 
             COmat_row=0d0 
!             call construct_CO_matrix_col_bare(PSI,VIN,COmat_bare_col,Nc,j)
             call construct_CO_matrix_row(PSI,VIN,only_bare,COmat_row,COmat_bare_row,Nc,j)

             OCmat_bare_col1(1:dimL_orb/2)=Conjg(COmat_bare_row(1:dimL_orb/2))
             OCmat_bare_col1(dimL_orb/2+1:dimL_orb)=-Conjg(COmat_bare_row(dimL_orb/2+1:dimL_orb))


             OCmat_bare_col2(1:dimL_orb/2)=COmat_bare_row(dimL_orb/2+1:dimL_orb)
             OCmat_bare_col2(dimL_orb/2+1:dimL_orb)=-COmat_bare_row(1:dimL_orb/2)




            do i=1, dimL_orb/2 
             Pb=0d0  
             pos=mod(i,ND)
             if(pos.eq.0) pos=ND
             call get_proj_col(PSI,Pb,i,pos)
             Pb=Conjg(Pb)

             left_Loc(i,j)=dot_product(Conjg(Pb),OCmat_bare_col1) 
             left_Loc(i,j+NC)=dot_product(Conjg(Pb),OCmat_bare_col2) 

              
            end do 


 
!             if(myid.eq.0) then
               do i=1, dimL_orb/2
                 if(abs(left_Loc(i,j)).gt.tol1) then
                   write(256,'(2I10,10x,2F20.10)') i,j,Real(left_Loc(i,j)),Dimag(left_Loc(i,j))
                 end if
                 if(abs(left_Loc(i,j+NC)).gt.tol1) then
                   write(256,'(2I10,10x,2F20.10)') i,j+NC,Real(left_Loc(i,j+NC)),Dimag(left_Loc(i,j+NC))
                 end if
               end do
!             end if 
!             call MPI_BARRIER(MPI_COMM_WORLD,ierr)


            j=j+1
           end do

           end if

698      continue

         if(myid.eq.0) close(256)


!        Get inverse squareroot density matrix
         call squarerootRInv_arnoldi(RhoInvSQR)
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        IF(MYID.eq.0) then

          allocate(OCmat_help(dimL_orb/2,2*NC))
          OCmat_help=0.0d0;
          DO k=1,Morb
             DO s=1,Nc
                DO i=1,ND
                   DO j=1,Morb

                  ! for the blocks
                  ind11=(/(j-1)*ND+i,s/)
                  ind12=(/(j-1)*ND+i,Nc+s/)
              !    ind21=(/Morb*ND+(j-1)*ND+i,s/)
              !    ind22=(/Morb*ND+(j-1)*ND+i,Nc+s/)
                  ind11a=(/(k-1)*ND+i,s/)
                  ind12a=(/(k-1)*ND+i,Nc+s/)
              !    ind21a=(/Morb*ND+(k-1)*ND+i,s/)
              !    ind22a=(/Morb*ND+(k-1)*ND+i,Nc+s/)

                  OCmat_help(ind11a(1),ind11a(2))=OCmat_help(ind11a(1),ind11a(2))+RhoInvSQR(k,j)*left_Loc(ind11(1),ind11(2))

                  OCmat_help(ind12a(1),ind12a(2))=OCmat_help(ind12a(1),ind12a(2))+RhoInvSQR(k,j)*left_Loc(ind12(1),ind12(2))

             !     OCmat_help(ind21a(1),ind21a(2))=OCmat_help(ind21a(1),ind21a(2))+CONJG(RhoInvSQR(k,j))*left_Loc(ind21(1),ind21(2))

             !     OCmat_help(ind22a(1),ind22a(2))=OCmat_help(ind22a(1),ind22a(2))+CONJG(RhoInvSQR(k,j))*left_Loc(ind22(1),ind22(2))

                  END DO
               END DO
            END DO
         END DO

         left_Loc=OCmat_help
         deallocate(OCmat_help)

       END IF
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)



       END SUBROUTINE construct_Loc_matrix_full


!=====  ONLY USE this ROUTINE IF MATRIX SMALL ENOUGH...
!=====  AND IF M >= 4 or 5 ...
       SUBROUTINE construct_Loc_matrix_full_sparse(PSI,VIN,NC)

          USE CI_ALL
          USE Parallel_CI 
          USE LR_ARNOLDI_MOD
          USE LR_RAPHA
          IMPLICIT NONE
          include 'mpif.h'

          COMPLEX*16, Dimension(Nc), INTENT(IN) :: VIN
!          COMPLEX*16, DIMENSION(ND,Morb), INTENT(IN) :: h2_PSI,h2_PSI_conjg
!          COMPLEX*16, INTENT(OUT) :: OCmat_row(2*Nc)
          COMPLEX*16, DIMENSION(ND,Morb), INTENT(IN) :: PSI
          COMPLEX*16 :: COmat_bare_col(2*NC)
          COMPLEX*16, DIMENSION(dimL_orb) :: Pb,COmat_row,COmat_bare_row,OCmat_bare_col1, OCmat_bare_col2
          COMPLEX*16, DIMENSION(Morb,Morb) :: RhoInvSQR
          INTEGER*4 :: i,j,k,s,l,q,Nc,pos,orb_pos,row,start_from, myid, nprocs, ierr, jmax, ind2, ll, jj
          LOGICAL ::  only_bare, exists
          Real*8 :: start,finish,a,b, tol1_Loc, t1, t2, szgb
          Complex*16, allocatable :: OCmat_help(:,:)
          INTEGER, DIMENSION(2) :: ind11, ind12, ind21, ind22
          INTEGER, DIMENSION(2) :: ind11a, ind12a, ind21a, ind22a
          COMPLEX*16, DIMENSION(dimL_orb/2,2) :: left_Loc_hilf,left_Loc_hilf_p0
          COMPLEX*16, DIMENSION(2*NC) :: test
          COMPLEX*16, allocatable :: Loc_vals_hilf2(:)
          Complex*16, allocatable ::  inter_mat_small(:,:),inter2_mat_small(:,:)
          INTEGER*4, allocatable :: Loc_rows_hilf2(:),Loc_cols_hilf2(:)
          Complex*16, allocatable ::  vals_inter(:)
          Complex*16, allocatable ::  inter2_mat(:,:)
          Complex*16, allocatable ::  c_mat(:,:)
          integer*4, allocatable    ::  rows_inter(:)
          integer*4, allocatable    ::  cols_inter(:)
!          COMPLEX*16, allocatable :: Loc_vals_csr(:)
          COMPLEX*16, allocatable :: Loc_inter(:,:), Loc_inter2(:,:)
!          INTEGER*8, allocatable :: Loc_ind_csr(:),Loc_col_csr(:)
          INTEGER :: s_status(MPI_STATUS_SIZE), st, info, rows_proj
          INTEGER, DIMENSION(6) :: mjob
          COMPLEX*16 :: hilf, hilf2, zdotc
          INTEGER*4 :: counter, Loc_maxnonzero
          INTEGER, allocatable :: ind(:)


          call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
          call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )

          tol1_Loc=tol1
          allocate(ind(nprocs))
          rows_proj=(dimL_orb/2)/nprocs
!!          Loc_maxnonzero=2000000000
          Loc_maxnonzero=maxnonzero

          if(myid.eq.0) then 

           INQUIRE(FILE="left_Loc.dat", EXIST=exists)

            if(.not.allocated(inter2_mat_small)) &
     &       allocate(inter2_mat_small(dimL_orb/2,2))


          if(exists.eqv..FALSE.) then  !files do not exist       

!!            if(.not.allocated(Loc_vals_hilf)) &
!!     &       allocate(Loc_vals_hilf(Loc_maxnonzero))
            if(.not.allocated(Loc_vals_hilf)) then
          allocate(Loc_vals_hilf(Loc_maxnonzero), stat=ierr)
          szgb= Loc_maxnonzero*32.0d0/1024/1024/1024
          write(*,*)"Size Loc_vals_hilf Gigabytes:", szgb, Loc_maxnonzero
          if(ierr /= 0) then
          write(*,*)"allocation error for Loc_vals_hilf"
          stop "Err alloc Loc_vals_hilf"
          else
          if(ierr == 0)write(*,*)"allocation ok for Loc_vals_hilf"
          endif
           endif

!!            if(.not.allocated(Loc_rows_hilf)) &
!!     &       allocate(Loc_rows_hilf(Loc_maxnonzero))
            if(.not.allocated(Loc_rows_hilf)) then
            allocate(Loc_rows_hilf(Loc_maxnonzero), stat=ierr)
          szgb= Loc_maxnonzero*16*2.0/1024/1024/1024
          write(*,*)"Size Loc_rows_hilf Gigabytes:", szgb
          if(ierr /= 0) then
          write(*,*)"allocation error for Loc_rows_hilf"
          stop "Err alloc Loo_rows_hilf"
          else
          if(ierr == 0)write(*,*)"allocation ok for Loo_rows_hilf"
          endif
            endif
 
!!            if(.not.allocated(Loc_cols_hilf)) &
!!     &       allocate(Loc_cols_hilf(Loc_maxnonzero))
            if(.not.allocated(Loc_cols_hilf)) then
           allocate(Loc_cols_hilf(Loc_maxnonzero), stat=ierr)
          szgb= Loc_maxnonzero*16*2.0/1024/1024/1024
          write(*,*)"Size Loc_cols_hilf Gigabytes:", szgb
          if(ierr /= 0) then
          write(*,*)"allocation error for Loc_cols_hilf"
          stop "Err alloc Loc_cols_hilf"
          else
          if(ierr == 0)write(*,*)"allocation ok for Loc_cols_hilf"
          endif
            endif

!!            if(.not.allocated(Loc_vals_hilf)) &
!!     &       allocate(Loc_vals_hilf(Loc_maxnonzero))
!!            if(.not.allocated(Loc_rows_hilf)) &
!!     &       allocate(Loc_rows_hilf(Loc_maxnonzero))
!!            if(.not.allocated(Loc_cols_hilf)) &
!!     &       allocate(Loc_cols_hilf(Loc_maxnonzero))
 
            open(unit=256, file='left_Loc.dat',status='replace', action='readwrite')

            start_from=1
            ind=1

          else

            open(unit=256, file='left_Loc.dat',status='old', action='readwrite')

!!            if(.not.allocated(Loc_vals_hilf)) &
!!     &       allocate(Loc_vals_hilf(Loc_maxnonzero))
            if(.not.allocated(Loc_vals_hilf)) then
          allocate(Loc_vals_hilf(Loc_maxnonzero), stat=ierr)
          szgb= Loc_maxnonzero*32.0d0/1024/1024/1024
          write(*,*)"Size Loc_vals_hilf Gigabytes:", szgb, Loc_maxnonzero
          if(ierr /= 0) then
          write(*,*)"allocation error for Loc_vals_hilf"
          stop "Err alloc Loc_vals_hilf"
          else
          if(ierr == 0)write(*,*)"allocation ok for Loc_vals_hilf"
          endif
           endif

!!            if(.not.allocated(Loc_rows_hilf)) &
!!     &       allocate(Loc_rows_hilf(Loc_maxnonzero))
            if(.not.allocated(Loc_rows_hilf)) then
            allocate(Loc_rows_hilf(Loc_maxnonzero), stat=ierr)
          szgb= Loc_maxnonzero*16*2.0/1024/1024/1024
          write(*,*)"Size Loc_rows_hilf Gigabytes:", szgb
          if(ierr /= 0) then
          write(*,*)"allocation error for Loc_rows_hilf"
          stop "Err alloc Loo_rows_hilf"
          else
          if(ierr == 0)write(*,*)"allocation ok for Loo_rows_hilf"
          endif
            endif
 
!!            if(.not.allocated(Loc_cols_hilf)) &
!!     &       allocate(Loc_cols_hilf(Loc_maxnonzero))
            if(.not.allocated(Loc_cols_hilf)) then
           allocate(Loc_cols_hilf(Loc_maxnonzero), stat=ierr)
          szgb= Loc_maxnonzero*16*2.0/1024/1024/1024
          write(*,*)"Size Loc_cols_hilf Gigabytes:", szgb
          if(ierr /= 0) then
          write(*,*)"allocation error for Loc_cols_hilf"
          stop "Err alloc Loc_cols_hilf"
          else
          if(ierr == 0)write(*,*)"allocation ok for Loc_cols_hilf"
          endif
            endif

!!            if(.not.allocated(Loc_vals_hilf)) &
!!     &       allocate(Loc_vals_hilf(Loc_maxnonzero))
!!            if(.not.allocated(Loc_rows_hilf)) &
!!     &       allocate(Loc_rows_hilf(Loc_maxnonzero))
!!            if(.not.allocated(Loc_cols_hilf)) &
!!     &       allocate(Loc_cols_hilf(Loc_maxnonzero))

            write(*,*) "Loc_maxnonzero", Loc_maxnonzero, maxnonzero

            jmax=1   
            ind=1 
            do 
              read(256,*,iostat=st) i,j,a,b
              if(st.ne.0) exit
               Loc_rows_hilf(ind)=i
               Loc_cols_hilf(ind)=j
               Loc_vals_hilf(ind)=Cmplx(a,b,kind=8) 
               ind(1)=ind(1)+1
               if(j.gt.jmax) jmax=j
!              left_Loc(i,j)=Cmplx(a,b,kind=8)
!               if(mod(j,NC).gt.jmax) jmax=mod(j,NC)  
!               if(mod(j,NC).eq.NC) jmax=NC
            end do
            start_from=jmax-NC  
            write(*,*) jmax 

          end if
          end if

          CALL MPI_BCAST(exists,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(start_from,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

          if(myid.ne.0) ind=1

          if(start_from.eq.NC) goto 698

          j=start_from
         call squarerootRInv_arnoldi(RhoInvSQR)


          if(.not.allocated(vals_inter)) allocate(vals_inter(dimL_orb)) 
          if(.not.allocated(rows_inter)) allocate(rows_inter(dimL_orb)) 
          if(.not.allocated(cols_inter)) allocate(cols_inter(dimL_orb)) 


            if(.not.allocated(inter2_mat)) &
     &       allocate(inter2_mat(dimL_orb/2,2*nprocs))

            if(myid.eq.0) then
             if(.not.allocated(c_mat)) &
     &        allocate(c_mat(rows_proj+mod(dimL_orb/2,nprocs),2*nprocs))
            else 
             if(.not.allocated(c_mat)) &
     &        allocate(c_mat(rows_proj,2*nprocs))
            end if            


!         Calculate Proj*L_oc
           only_bare=.TRUE.
       call CPU_TIME(t1)
           do while(j.le.NC-nprocs+1) 
             COmat_bare_row=0d0 
             COmat_row=0d0 
             vals_inter=0d0
             cols_inter=0
             rows_inter=0

             call construct_CO_matrix_row(PSI,VIN,only_bare,COmat_row,COmat_bare_row,Nc,j+myid)

             OCmat_bare_col1(1:dimL_orb/2)=Conjg(COmat_bare_row(1:dimL_orb/2))
!             OCmat_bare_col1(dimL_orb/2+1:dimL_orb)=-Conjg(COmat_bare_row(dimL_orb/2+1:dimL_orb))


             OCmat_bare_col2(1:dimL_orb/2)=COmat_bare_row(dimL_orb/2+1:dimL_orb)
!             OCmat_bare_col2(dimL_orb/2+1:dimL_orb)=-COmat_bare_row(1:dimL_orb/2)


!            do i=1, dimL_orb/2 
!             Pb=0d0  
!             pos=mod(i,ND)
!             if(pos.eq.0) pos=ND
!             call get_proj_col(PSI,Pb,i,pos)
!             Pb=Conjg(Pb)
!                left_Loc_hilf(i,1)=zdotc(dimL_orb,Pb,1,OCmat_bare_col1,1) 
!                left_Loc_hilf(i,2)=zdotc(dimL_orb,Pb,1,OCmat_bare_col2,1) 
!            end do 


       inter2_mat=0d0
       inter2_mat(1:dimL_orb/2,1+myid)=OCmat_bare_col1(1:dimL_orb/2)  
       inter2_mat(1:dimL_orb/2,1+myid+nprocs)=OCmat_bare_col2(1:dimL_orb/2)  

       do k=1, nprocs
            call MPI_BCAST(inter2_mat(:,1+(k-1)), &
     &                dimL_orb/2, &
     &                MPI_DOUBLE_COMPLEX,k-1,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(inter2_mat(:,1+nprocs+(k-1)), &
     &                dimL_orb/2, &
     &                MPI_DOUBLE_COMPLEX,k-1,MPI_COMM_WORLD,ierr)
       end do


       c_mat=0d0
        if(myid.eq.0) then
               do i=1, 2*nprocs
                call mkl_zcoogemv('N', rows_proj+mod(dimL_orb/2,nprocs), projdist_vals_coo, projdist_rows_coo, &
     &                   projdist_cols_coo, proj_nonzeros(1), inter2_mat(:,i), &
     &                   c_mat(:,i))
               end do 

       else

         do i=1, 2*nprocs
                call mkl_zcoogemv('N', rows_proj, projdist_vals_coo, projdist_rows_coo, &
     &                   projdist_cols_coo, proj_nonzeros(myid+1), inter2_mat(:,i), &
     &                   c_mat(:,i))
         end do

        end if
!       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!       stop       



       inter2_mat=0d0
       if(myid.eq.0) then
          do i=1, 2*nprocs
              inter2_mat(1: &
     &                   mod(dimL_orb/2,nprocs)+rows_proj,i)=c_mat(1: &
     &                   mod(dimL_orb/2,nprocs)+rows_proj,i)
          end do
       end if
         do i=1, 2*nprocs
          call MPI_BCAST(inter2_mat(1: &
     &                   mod(dimL_orb/2,nprocs)+rows_proj,i), &
     &                rows_proj+mod(dimL_orb/2,nprocs), &
     &                MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
         end do
       do k=2, nprocs
            if(myid.eq.k-1) then
             do i=1, 2*nprocs
              inter2_mat(1+mod(dimL_orb/2,nprocs)+(k-1)*rows_proj: &
     &                   mod(dimL_orb/2,nprocs)+k*rows_proj,i)= & 
     &             c_mat(1: &
     &                   rows_proj,i)
             end do
            end if 
            do i=1, 2*nprocs
             call MPI_BCAST(inter2_mat(1+mod(dimL_orb/2,nprocs)+(k-1)*rows_proj: &
     &                   mod(dimL_orb/2,nprocs)+k*rows_proj,i), &
     &                rows_proj, &
     &                MPI_DOUBLE_COMPLEX,k-1,MPI_COMM_WORLD,ierr)
            end do
       end do 


       left_Loc_hilf(:,1)=inter2_mat(:,1+myid)
       left_Loc_hilf(:,2)=inter2_mat(:,1+myid+nprocs)

              do k=2, nprocs
                ind(k)=1
              end do 
            
              do i=1, ND
                do k=1, Morb
                 hilf=0d0; hilf2=0d0; 
                 do ll=1, Morb
                  hilf=hilf+RhoInvSQR(k,ll)*left_Loc_hilf(i+(ll-1)*ND,1)  
                  hilf2=hilf2+RhoInvSQR(k,ll)*left_Loc_hilf(i+(ll-1)*ND,2)  
                 end do

             IF(myid.eq.0) then
                 if(abs(hilf).gt.tol1_Loc) then
                   Loc_rows_hilf(ind(1))=i+(k-1)*ND
                   Loc_cols_hilf(ind(1))=j
                   Loc_vals_hilf(ind(1))=hilf
!                   write(256,'(2I10,10x,2F20.10)') i+(k-1)*ND,j,Real(hilf),Dimag(hilf)
                   ind(1)=ind(1)+1 
                 end if
                 if(abs(hilf2).gt.tol1_Loc) then
                   Loc_rows_hilf(ind(1))=i+(k-1)*ND
                   Loc_cols_hilf(ind(1))=j+NC
                   Loc_vals_hilf(ind(1))=hilf2
!                   write(256,'(2I10,10x,2F20.10)') i+(k-1)*ND,j+NC,Real(hilf2),Dimag(hilf2)
                   ind(1)=ind(1)+1 
                 end if
             ELSE 
                 if(abs(hilf).gt.tol1_Loc) then
                   rows_inter(ind(myid+1))=i+(k-1)*ND
                   cols_inter(ind(myid+1))=j+myid
                   vals_inter(ind(myid+1))=hilf
                   ind(myid+1)=ind(myid+1)+1 
                 end if
                 if(abs(hilf2).gt.tol1_Loc) then
                   rows_inter(ind(myid+1))=i+(k-1)*ND
                   cols_inter(ind(myid+1))=j+myid+NC
                   vals_inter(ind(myid+1))=hilf2
                   ind(myid+1)=ind(myid+1)+1 
                 end if
             END IF 
                end do
              end do


            call MPI_BARRIER(MPI_COMM_WORLD,ierr)

            do k=2, nprocs

!              if(myid.eq.k-1) then
!                 left_Loc_hilf_p0(:,1)=left_Loc_hilf(:,1) 
!                 left_Loc_hilf_p0(:,2)=left_Loc_hilf(:,2) 
!              end if

              if(myid.eq.0) CALL MPI_RECV(ind(k),1,MPI_INTEGER,k-1,2,MPI_COMM_WORLD,s_status,ierr) 
              if(myid.eq.k-1) CALL MPI_SEND(ind(k),1,MPI_INTEGER,0,2,MPI_COMM_WORLD,ierr)

             if(ind(k).ne.1) then
              if(myid.eq.0) CALL MPI_RECV(vals_inter(1:ind(k)-1),ind(k)-1,MPI_DOUBLE_COMPLEX,k-1,0,MPI_COMM_WORLD,s_status,ierr) 
              if(myid.eq.k-1) CALL MPI_SEND(vals_inter(1:ind(k)-1),ind(k)-1,MPI_DOUBLE_COMPLEX,0,0,MPI_COMM_WORLD,ierr)
              if(myid.eq.0) CALL MPI_RECV(cols_inter(1:ind(k)-1),ind(k)-1,MPI_INTEGER,k-1,1,MPI_COMM_WORLD,s_status,ierr) 
              if(myid.eq.k-1) CALL MPI_SEND(cols_inter(1:ind(k)-1),ind(k)-1,MPI_INTEGER,0,1,MPI_COMM_WORLD,ierr)
              if(myid.eq.0) CALL MPI_RECV(rows_inter(1:ind(k)-1),ind(k)-1,MPI_INTEGER,k-1,3,MPI_COMM_WORLD,s_status,ierr) 
              if(myid.eq.k-1) CALL MPI_SEND(rows_inter(1:ind(k)-1),ind(k)-1,MPI_INTEGER,0,3,MPI_COMM_WORLD,ierr)


 
             if(myid.eq.0) then


              do i=1, ind(k)-1
                   Loc_rows_hilf(ind(1))=rows_inter(i)
                   Loc_cols_hilf(ind(1))=cols_inter(i)
                   Loc_vals_hilf(ind(1))=vals_inter(i)
!                   write(256,'(2I10,10x,2F20.10)') rows_inter(i),cols_inter(i), &
!     &                     Real(vals_inter(i)),Dimag(vals_inter(i))
                   ind(1)=ind(1)+1 
              end do
               vals_inter=0d0
               cols_inter=0
               rows_inter=0 
             end if


             end if 
             call MPI_BARRIER(MPI_COMM_WORLD,ierr)
            end do  

!            if(myid.eq.0) write(*,*) "Loc nonzeroes after",j, "cols:",ind(1) 
            j=j+nprocs
           end do


           if(j.gt.NC-nprocs+1.and.j.le.NC.and.myid.eq.0) then

           do while(j.le.NC) 
             COmat_bare_row=0d0 
             COmat_row=0d0 
             call construct_CO_matrix_row(PSI,VIN,only_bare,COmat_row,COmat_bare_row,Nc,j)

             OCmat_bare_col1(1:dimL_orb/2)=Conjg(COmat_bare_row(1:dimL_orb/2))
!             OCmat_bare_col1(dimL_orb/2+1:dimL_orb)=-Conjg(COmat_bare_row(dimL_orb/2+1:dimL_orb))


             OCmat_bare_col2(1:dimL_orb/2)=COmat_bare_row(dimL_orb/2+1:dimL_orb)
!             OCmat_bare_col2(dimL_orb/2+1:dimL_orb)=-COmat_bare_row(1:dimL_orb/2)


                inter2_mat_small=0d0
                call mkl_zcoogemv('N', dimL_orb/2, proj_vals_coo, proj_rows_coo, &
     &                   proj_cols_coo, proj_nonzero, OCmat_bare_col1(1:dimL_orb/2), &
     &                   inter2_mat_small(:,1))

                call mkl_zcoogemv('N', dimL_orb/2, proj_vals_coo, proj_rows_coo, &
     &                   proj_cols_coo, proj_nonzero, OCmat_bare_col2(1:dimL_orb/2), &
     &                   inter2_mat_small(:,2))

!            do i=1, dimL_orb/2 
!             Pb=0d0  
!             pos=mod(i,ND)
!             if(pos.eq.0) pos=ND
!             call get_proj_col(PSI,Pb,i,pos)
!             Pb=Conjg(Pb)
!                left_Loc_hilf(i,1)=zdotc(dimL_orb,Pb,1,OCmat_bare_col1,1) 
!                left_Loc_hilf(i,2)=zdotc(dimL_orb,Pb,1,OCmat_bare_col2,1) 
!            end do 


              left_Loc_hilf=0d0
              left_Loc_hilf(:,1)=inter2_mat_small(:,1)
              left_Loc_hilf(:,2)=inter2_mat_small(:,2)
 

              do i=1, ND
                do k=1, Morb
                 hilf=0d0; hilf2=0d0; 
                 do ll=1, Morb
                  hilf=hilf+RhoInvSQR(k,ll)*left_Loc_hilf(i+(ll-1)*ND,1)  
                  hilf2=hilf2+RhoInvSQR(k,ll)*left_Loc_hilf(i+(ll-1)*ND,2)  
                 end do
                 if(abs(hilf).gt.tol1_Loc) then
                   Loc_rows_hilf(ind(1))=i+(k-1)*ND
                   Loc_cols_hilf(ind(1))=j
                   Loc_vals_hilf(ind(1))=hilf
!                   write(256,'(2I10,10x,2F20.10)') i+(k-1)*ND,j,Real(hilf),Dimag(hilf)
                   ind(1)=ind(1)+1 
                 end if
                 if(abs(hilf2).gt.tol1_Loc) then
                   Loc_rows_hilf(ind(1))=i+(k-1)*ND
                   Loc_cols_hilf(ind(1))=j+NC
                   Loc_vals_hilf(ind(1))=hilf2
!                   write(256,'(2I10,10x,2F20.10)') i+(k-1)*ND,j+NC,Real(hilf2),Dimag(hilf2)
                   ind(1)=ind(1)+1 
                 end if
                end do
              end do
            


!            if(myid.eq.0) write(*,*) "Loc nonzeroes after",j, "cols:",ind(1) 
            j=j+1
           end do

           end if

       call CPU_TIME(t2)
       if(myid.eq.0) write(*,*) "5th:",t2-t1
       call MPI_Barrier(MPI_COMM_WORLD,ierr)


698      continue

         if(myid.eq.0) close(256)
         Loc_nonzeroes=ind(1)-1
         goto 12

!        Get inverse squareroot density matrix
!         call squarerootRInv_arnoldi(RhoInvSQR)
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        IF(MYID.eq.0) then

!       Conversion of left_Loc_inter from sparse COO to CSR format
        mjob(1)=2
        mjob(2)=1
        mjob(3)=1
        mjob(4)=2
        mjob(5)=ind(1)
        mjob(6)=0


        counter=ind(1)-1 
        write(*,*) ind(1)-1

           if(.not.allocated(Loc_vals_csr)) & 
     &                allocate(Loc_vals_csr(counter))
           if(.not.allocated(Loc_col_csr)) &
     &                allocate(Loc_col_csr(counter))
           if(.not.allocated(Loc_ind_csr)) &
     &                allocate(Loc_ind_csr(dimL_orb/2+1))

         

       
        Loc_col_csr=0; Loc_vals_csr=0d0; Loc_ind_csr=0;
        call mkl_zcsrcoo(mjob, dimL_orb/2, Loc_vals_csr, Loc_col_csr, Loc_ind_csr, &
     &             counter, Loc_vals_hilf,  &
     &             Loc_rows_hilf, Loc_cols_hilf, info) 

        write(*,*) info, counter


         if(allocated(Loc_rows_hilf)) deallocate(Loc_rows_hilf) 
         if(allocated(Loc_cols_hilf)) deallocate(Loc_cols_hilf) 
         if(allocated(Loc_vals_hilf)) deallocate(Loc_vals_hilf) 

       END IF

12      continue
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)

            if(allocated(inter2_mat_small)) &
     &       deallocate(inter2_mat_small)
          if(allocated(vals_inter)) deallocate(vals_inter) 
          if(allocated(rows_inter)) deallocate(rows_inter) 
          if(allocated(cols_inter)) deallocate(cols_inter) 
            if(allocated(inter2_mat)) &
     &       deallocate(inter2_mat)
            if(allocated(c_mat)) &
     &       deallocate(c_mat)



       END SUBROUTINE construct_Loc_matrix_full_sparse



!==================================================================================================

       SUBROUTINE construct_CO_matrix_row(PSI,VIN,only_bare,COmat_row,COmat_bare_row,Nc,row)

          USE CI_ALL
          USE Parallel_CI 
          USE DVR_ALL
          USE W_INTERPARTICLE 
          USE LR_RAPHA

          IMPLICIT NONE
          INCLUDE 'mpif.h'
          COMPLEX*16, INTENT(IN) :: VIN(Nc)
          COMPLEX*16, INTENT(OUT) :: COmat_row(2*ND*Morb)
          COMPLEX*16, INTENT(OUT) :: COmat_bare_row(2*ND*Morb)
          COMPLEX*16  :: PL(2*ND*Morb)
          COMPLEX*16, DIMENSION(ND,Morb), INTENT(IN) :: PSI
          COMPLEX*16 :: Vhelp1(Nc), Vhelp2(Nc)
          COMPLEX*16 :: Vhelp1b(Nc), Vhelp2b(Nc)
          COMPLEX*16 :: Vhelp1a(Nc), Vhelp2a(Nc)
!          COMPLEX*16, DIMENSION(ND,Morb), INTENT(IN) :: h2_PSI,h2_PSI_conjg
          COMPLEX*16 :: COmat_help(2*ND*Morb)
          COMPLEX*16, DIMENSION(Morb,Morb) :: Rho_arnoldi
          INTEGER, DIMENSION(2) :: ind11, ind12, ind21, ind22
          INTEGER, DIMENSION(2) :: ind11a, ind12a, ind21a, ind22a
          COMPLEX*16, DIMENSION(dimL_orb) :: Pb, Pbb, PLL, inter, final, Pb2
          COMPLEX*16, DIMENSION(dimL_orb/2) :: P
          COMPLEX*16, allocatable :: Lco_inter(:,:),A(:,:),B(:,:),A_small(:,:),B_small(:,:)
          COMPLEX*16, allocatable :: store(:)
          COMPLEX*16 :: VIN_help(Nc)
          COMPLEX*16 :: rhohelp
         COMPLEX*16, DIMENSION(NDX*NDY*NDZ)  :: WSL !ADDED BY STR
          LOGICAL, INTENT(IN)                :: only_bare
          Real*8 :: start,finish
          INTEGER :: s_status(MPI_STATUS_SIZE)
!==========================================================
          INTEGER :: i,j,k,s,l,q,Nc,row,pos,CI_pos,orb_pos2,myid,ierr,nprocs,row0,rows_proj
!===========================================================
           xlambda0=xlambda_0

!         Initialize

          call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
          call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )

          rows_proj=(diml_orb/2)/nprocs 
          if(.not.allocated(store)) allocate(store(rows_proj))


          COmat_row=0.0d0; Vhelp1a=0.0d0; Vhelp2a=0.0d0; Vhelp1b=0.0d0; Vhelp2b=0.0d0

          CI_pos=mod(row,NC)
!           write(*,*) "CI_pos",row,NC
          if(CI_pos.eq.0) CI_pos=NC 
!           write(*,*) "CI_pos",CI_pos,row,NC


!         UPPER LEFT + UPPER RIGHT        
          if(preconstr.eqv..FALSE.) then 
!             allocate(h2_PSI(NDX*NDY*NDZ,Morb)) 
!             allocate(h2_PSI_conjg(NDX*NDY*NDZ,Morb)) 
             CALL get_h2_PSI(PSI,h2_PSI)
             CALL get_h2_PSI(Conjg(PSI),h2_PSI_conjg)
          end if

          DO q=1,Morb
             DO k=1,Morb

!               1-body
               if(V_CIJ_prec.eqv..FALSE.) then 
                VIN_help=VIN; Vhelp1a=0.0d0
                CALL Produce_Cij(VIN_help,Vhelp1a,k,q,0,0,Nc) !Correct - proved by STR
                VIN_help=VIN; Vhelp1b=0.0d0
                CALL Produce_Cij(VIN_help,Vhelp1b,q,k,0,0,Nc) !Correct - proved by STR
               else
                 Vhelp1a=0.0d0
                 Vhelp1b=0.0d0
!                 if(myid.eq.0) then
                    Vhelp1a=V_PrdCIJ2(:,k,q)  
                    Vhelp1b=V_PrdCIJ2(:,q,k)  
!                 end if
!                 call MPI_BCAST(Vhelp1a,ND,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
!                 call MPI_BCAST(Vhelp1b,ND,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
               end if

!                  write(*,*) "Vhelp1a(row)",row,Vhelp1a(CI_pos)
!RB Feb29                rhohelp=SUM(Conjg(VIN_help)*Vhelp1a)

!RB                DO i=1,Nc

                 if(abs(Vhelp1a(CI_pos)).ne.0d0.and.abs(Vhelp1b(CI_pos)).ne.0d0) then 
                   DO j=1,ND
                      COmat_row((q-1)*ND+j)=COmat_row((q-1)*ND+j)+Vhelp1a(CI_pos)*CONJG(h2_PSI(j,k))
                      COmat_row(Morb*ND+(q-1)*ND+j)=COmat_row(Morb*ND+(q-1)*ND+j)+Vhelp1b(CI_pos)*CONJG(h2_PSI_conjg(j,k))
                   END DO
                 else
                    if(abs(Vhelp1a(CI_pos)).ne.0d0) then 
                      DO j=1,ND
                        COmat_row((q-1)*ND+j)=COmat_row((q-1)*ND+j)+Vhelp1a(CI_pos)*CONJG(h2_PSI(j,k))
                      END DO
                    end if
                    if(abs(Vhelp1b(CI_pos)).ne.0d0) then 
                      DO j=1,ND
                        COmat_row(Morb*ND+(q-1)*ND+j)=COmat_row(Morb*ND+(q-1)*ND+j)+Vhelp1b(CI_pos)*CONJG(h2_PSI_conjg(j,k))
                      END DO
                    end if
                 end if 
!RB                END DO


!               2-body
              if(xlambda0.ne.0d0) then
!                     call CPU_TIME(start)
                DO s=1,Morb
                   DO l=1,Morb
!                      IF(k==l.AND.s==k) THEN
!                      IF(q==s.AND.k==l.AND.q==k) GOTO 111
                    if(V_CIJ_prec.eqv..FALSE.) then 
                      VIN_help=VIN; Vhelp2a=0.0d0
                      CALL Produce_Cij(VIN_help,Vhelp2a,k,s,q,l,Nc) !Correct - proved by STR
                      VIN_help=VIN; Vhelp2b=0.0d0
                      CALL Produce_Cij(VIN_help,Vhelp2b,q,s,k,l,Nc)
                    else
                      Vhelp2a=0.0d0
                      Vhelp2b=0.0d0
!                      if(myid.eq.0) then
                        Vhelp2a=V_PrdCIJ(:,k,s,q,l)  
                        Vhelp2b=V_PrdCIJ(:,q,s,k,l)  
!                      end if
!                      call MPI_BCAST(Vhelp2a,ND,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
!                      call MPI_BCAST(Vhelp2b,ND,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
                    end if

!RB Feb29                rhohelp=SUM(Conjg(VIN)*Vhelp2a)
!                  write(*,*) "Vhelp2a(row)",row,Vhelp2a(CI_pos)
                 
            if(WSL_prec.eqv..FALSE.) then
               CALL  Get_WSL_omp(WSL,psi(:,s),psi(:,l))  
            else
               WSL=WSL_mat(:,s,l)
            end if 
!            CALL  Get_WSL_omp(WSL,psi(:,s),psi(:,l))

!RB                      DO i=1,Nc


                 if(abs(Vhelp2a(CI_pos)).ne.0d0.and.abs(Vhelp2b(CI_pos)).ne.0d0) then 
                         DO j=1,ND
                            COmat_row((q-1)*ND+j)=COmat_row((q-1)*ND+j)+Vhelp2a(CI_pos)*CONJG(PSI(j,k))*WSL(j)
                            COmat_row(Morb*ND+(q-1)*ND+j)=COmat_row(Morb*ND+(q-1)*ND+j)+Vhelp2b(CI_pos)*PSI(j,k)*WSL(j)
                         END DO
                 else
                    if(abs(Vhelp2a(CI_pos)).ne.0d0) then 
                         DO j=1,ND
                            COmat_row((q-1)*ND+j)=COmat_row((q-1)*ND+j)+Vhelp2a(CI_pos)*CONJG(PSI(j,k))*WSL(j)
                         END DO
                    end if 
                    if(abs(Vhelp2b(CI_pos)).ne.0d0) then 
                         DO j=1,ND
                            COmat_row(Morb*ND+(q-1)*ND+j)=COmat_row(Morb*ND+(q-1)*ND+j)+Vhelp2b(CI_pos)*PSI(j,k)*WSL(j)
                         END DO
                    end if 
                 end if

!RB                      END DO

                   END DO
                END DO
!                     call CPU_TIME(finish)
!                     write(*,*) finish-start 
               end if

             END DO
          END DO

!          if(preconstr.eqv..FALSE.) then 
!             deallocate(h2_PSI) 
!             deallocate(h2_PSI_conjg) 
!          end if

 
          if(row.le.NC) then
            COmat_bare_row=COmat_row
          else
            COmat_bare_row(1:dimL_orb/2)=-Conjg(COmat_row(dimL_orb/2+1:dimL_orb))
            COmat_bare_row(dimL_orb/2+1:dimL_orb)=-Conjg(COmat_row(1:dimL_orb/2))
          end if
          COmat_row=COmat_bare_row 
!        NOW bare row calculated!!!

       if(only_bare) then
         return
       else  

        if(myid.eq.0) row0=row
        call MPI_BCAST(row0,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)        
 
        IF(row0.le.NC-nprocs+1) THEN 

!        inter=COmat_row
!        goto 21 

         if(.not.allocated(Lco_inter)) allocate(Lco_inter(dimL_orb/2,2*nprocs)) 
         if(myid.eq.0) then
           if(.not.allocated(A)) allocate(A(rows_proj+mod(dimL_orb/2,nprocs),nprocs)) 
           if(.not.allocated(B)) allocate(B(rows_proj+mod(dimL_orb/2,nprocs),nprocs)) 
         else
           if(.not.allocated(A)) allocate(A(rows_proj,nprocs)) 
           if(.not.allocated(B)) allocate(B(rows_proj,nprocs)) 
         end if

         Lco_inter=0d0
         Lco_inter(1:dimL_orb/2,myid+1)=COmat_row(1:dimL_orb/2) 
         Lco_inter(1:dimL_orb/2,nprocs+myid+1)=COmat_row(1+dimL_orb/2:dimL_orb/2) 

         do k=1, nprocs
            call MPI_BCAST(Lco_inter(:,k), &
     &                dimL_orb/2, &
     &                MPI_DOUBLE_COMPLEX,k-1,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(Lco_inter(:,k+nprocs), &
     &                dimL_orb/2, &
     &                MPI_DOUBLE_COMPLEX,k-1,MPI_COMM_WORLD,ierr)
         end do


       A=0d0
       B=0d0
        if(myid.eq.0) then
               do i=1, nprocs
                call mkl_zcoogemv('N', rows_proj+mod(dimL_orb/2,nprocs), projdist_vals_coo, projdist_rows_coo, &
     &                   projdist_cols_coo, proj_nonzeros(1), Lco_inter(:,i), &
     &                   A(:,i))

                call mkl_zcoogemv('N', rows_proj+mod(dimL_orb/2,nprocs), Conjg(projdist_vals_coo), projdist_rows_coo, &
     &                   projdist_cols_coo, proj_nonzeros(1), Lco_inter(:,i+nprocs), &
     &                   B(:,i))
               end do 
        else
               do i=1, nprocs
                  call mkl_zcoogemv('N', rows_proj, projdist_vals_coo, projdist_rows_coo, &
     &                   projdist_cols_coo, proj_nonzeros(myid+1), Lco_inter(:,i), &
     &                   A(:,i))

                  call mkl_zcoogemv('N', rows_proj, Conjg(projdist_vals_coo), projdist_rows_coo, &
     &                   projdist_cols_coo, proj_nonzeros(myid+1), Lco_inter(:,i+nprocs), &
     &                   B(:,i))
               end do
        end if


       Lco_inter=0d0
       if(myid.eq.0) then
          do i=1, nprocs
              Lco_inter(1: &
     &                   mod(dimL_orb/2,nprocs)+rows_proj,i)=A(1: &
     &                   mod(dimL_orb/2,nprocs)+rows_proj,i)

              Lco_inter(1: &
     &                   mod(dimL_orb/2,nprocs)+rows_proj,i+nprocs)=B(1: &
     &                   mod(dimL_orb/2,nprocs)+rows_proj,i)
          end do
       end if
         do i=1, nprocs
          call MPI_BCAST(Lco_inter(1: &
     &                   mod(dimL_orb/2,nprocs)+rows_proj,i), &
     &                rows_proj+mod(dimL_orb/2,nprocs), &
     &                MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)

          call MPI_BCAST(Lco_inter(1: &
     &                   mod(dimL_orb/2,nprocs)+rows_proj,i+nprocs), &
     &                rows_proj+mod(dimL_orb/2,nprocs), &
     &                MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
         end do
       do k=2, nprocs
            if(myid.eq.k-1) then
             do i=1, nprocs
              Lco_inter(1+mod(dimL_orb/2,nprocs)+(k-1)*rows_proj: &
     &                   mod(dimL_orb/2,nprocs)+k*rows_proj,i)= & 
     &             A(1: &
     &                   rows_proj,i)

              Lco_inter(1+mod(dimL_orb/2,nprocs)+(k-1)*rows_proj: &
     &                   mod(dimL_orb/2,nprocs)+k*rows_proj,i+nprocs)= & 
     &             B(1: &
     &                   rows_proj,i)
             end do
            end if 
            do i=1, nprocs
             call MPI_BCAST(Lco_inter(1+mod(dimL_orb/2,nprocs)+(k-1)*rows_proj: &
     &                   mod(dimL_orb/2,nprocs)+k*rows_proj,i), &
     &                rows_proj, &
     &                MPI_DOUBLE_COMPLEX,k-1,MPI_COMM_WORLD,ierr)

             call MPI_BCAST(Lco_inter(1+mod(dimL_orb/2,nprocs)+(k-1)*rows_proj: &
     &                   mod(dimL_orb/2,nprocs)+k*rows_proj,i+nprocs), &
     &                rows_proj, &
     &                MPI_DOUBLE_COMPLEX,k-1,MPI_COMM_WORLD,ierr)
            end do
       end do 



          inter(1:dimL_orb/2)=Lco_inter(:,1+myid)
          inter(1+dimL_orb/2:dimL_orb)=Lco_inter(:,1+nprocs+myid)

!21       continue 


!      get squareroot of inverse of rho
       CALL squarerootRInv_arnoldi(Rho_arnoldi)

!      Multiply with inverse sqaureroot of Rho
          final=0.0d0;

            do i=1, Morb 
             do j=1, ND
              do k=1, Morb
                final((i-1)*ND+j)=final((i-1)*ND+j)+Rho_arnoldi(k,i)*inter((k-1)*ND+j)
                final(Morb*ND+(i-1)*ND+j)=final(Morb*ND+(i-1)*ND+j)+Conjg(Rho_arnoldi(k,i))*inter(Morb*ND+(k-1)*ND+j)
              end do 
             end do 
            end do

         COmat_row=final

!           open(unit=256, file='test.dat', status='replace', action='readwrite')
!           do q=1, Morb
!           do i=1, dimL_orb
!            write(256,'(500(2F20.10))') (Real(Tkin(i,j)),Dimag(Tkin(i,j)), &
!     &        j=1,nd,1)
!             write(256,'(500(2F20.10))') Real(COmat_row(i)),Dimag(COmat_row(i))
!           end do
!           end do 
!          close(256)

        ELSE

!        inter=COmat_row
!        goto 22 


         if(myid.eq.0) then
           if(.not.allocated(A_small)) allocate(A_small(rows_proj+mod(dimL_orb/2,nprocs),2)) 
           if(.not.allocated(B_small)) allocate(B_small(rows_proj+mod(dimL_orb/2,nprocs),2)) 
         else
           if(.not.allocated(A_small)) allocate(A_small(rows_proj,2)) 
           if(.not.allocated(B_small)) allocate(B_small(rows_proj,2)) 
         end if

       A_small=0d0
       B_small=0d0
        if(myid.eq.0) then
                call mkl_zcoogemv('N', rows_proj+mod(dimL_orb/2,nprocs), projdist_vals_coo, projdist_rows_coo, &
     &                   projdist_cols_coo, proj_nonzeros(1), COmat_row(1:dimL_orb/2), &
     &                   A_small(:,1))

                call mkl_zcoogemv('N', rows_proj+mod(dimL_orb/2,nprocs), Conjg(projdist_vals_coo), projdist_rows_coo, &
     &                   projdist_cols_coo, proj_nonzeros(1), COmat_row(1+dimL_orb/2:dimL_orb), &
     &                   A_small(:,2))
        else
                  call mkl_zcoogemv('N', rows_proj, projdist_vals_coo, projdist_rows_coo, &
     &                   projdist_cols_coo, proj_nonzeros(myid+1), COmat_row(1:dimL_orb/2), &
     &                   A_small(:,1))

                  call mkl_zcoogemv('N', rows_proj, Conjg(projdist_vals_coo), projdist_rows_coo, &
     &                   projdist_cols_coo, proj_nonzeros(myid+1), COmat_row(1+dimL_orb/2:dimL_orb), &
     &                   A_small(:,2))
        end if



            if(myid.eq.0) then
                 inter(1:rows_proj+mod(dimL_orb/2,nprocs))=A_small(:,1)
                 inter(1+dimL_orb/2:dimL_orb/2+rows_proj+mod(dimL_orb/2,nprocs))=A_small(:,2)
            end if

           do k=2, nprocs
             if(myid.eq.k-1) then
                  store(:)= &  
     &            A_small(:,1) 

              CALL MPI_SEND(store,rows_proj, &
     &         MPI_DOUBLE_COMPLEX,0,0,MPI_COMM_WORLD,ierr)
             end if 
!              CALL MPI_BCAST(inter3(1+(k-1)*rows_proj+mod(dimL_orb/2,nprocs):k*rows_proj+mod(dimL_orb/2,nprocs),j),rows_proj, &
!     &         MPI_DOUBLE_COMPLEX,k-1,MPI_COMM_WORLD,ierr)

             if(myid.eq.0) then
              CALL MPI_RECV(store,rows_proj, &
     &         MPI_DOUBLE_COMPLEX,k-1,0,MPI_COMM_WORLD,s_status,ierr)

              inter(1+(k-1)*rows_proj+mod(dimL_orb/2,nprocs):k*rows_proj+mod(dimL_orb/2,nprocs))= &
     &          store(:)       
             end if

             if(myid.eq.k-1) then
                  store(:)= &
     &            A_small(:,2) 

              CALL MPI_SEND(store,rows_proj, &
     &         MPI_DOUBLE_COMPLEX,0,1,MPI_COMM_WORLD,ierr)

             end if
!              CALL MPI_BCAST(inter3(1+dimL_orb/2+(k-1)*rows_proj+mod(dimL_orb/2,nprocs):dimL_orb/2+k*rows_proj+mod(dimL_orb/2,nprocs),j),rows_proj, &
!     &         MPI_DOUBLE_COMPLEX,k-1,MPI_COMM_WORLD,ierr)

             if(myid.eq.0) then
                CALL MPI_RECV(store,rows_proj, &
     &         MPI_DOUBLE_COMPLEX,k-1,1,MPI_COMM_WORLD,s_status,ierr)

              inter(1+dimL_orb/2+(k-1)*rows_proj+mod(dimL_orb/2,nprocs):dimL_orb/2+k*rows_proj+mod(dimL_orb/2,nprocs))= &
     &           store(:)
             end if

            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
           end do

!22     continue

       if(myid.eq.0) then
!      get squareroot of inverse of rho
       CALL squarerootRInv_arnoldi(Rho_arnoldi)

!      Multiply with inverse sqaureroot of Rho
          final=0.0d0;

            do i=1, Morb 
             do j=1, ND
              do k=1, Morb
                final((i-1)*ND+j)=final((i-1)*ND+j)+Rho_arnoldi(k,i)*inter((k-1)*ND+j)
                final(Morb*ND+(i-1)*ND+j)=final(Morb*ND+(i-1)*ND+j)+Conjg(Rho_arnoldi(k,i))*inter(Morb*ND+(k-1)*ND+j)
              end do 
             end do 
            end do

           COmat_row=final
        end if

            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        END IF

        end if 

       END SUBROUTINE construct_CO_matrix_row




!==========action of h2 on PSI      
       SUBROUTINE get_h2_PSI(PSI,h2_PSI_local)

          USE CI_ALL
          USE Parallel_CI 
          USE DVR_ALL
          USE LR_RAPHA

          IMPLICIT NONE
          include 'mpif.h'

          COMPLEX*16, DIMENSION(ND,Morb), INTENT(IN) :: PSI
          COMPLEX*16, DIMENSION(ND,Morb), INTENT(OUT) :: h2_PSI_local
          COMPLEX*16, DIMENSION(ND)          :: Tkin_col,Tkin_row
!          COMPLEX*16, DIMENSION(dimL_orb,dimL_orb) :: L_orb
          COMPLEX*16, DIMENSION(dimL_orb/2) :: h2
          INTEGER :: i,j,q,myid,nprocs,ierr
          INTEGER :: s_status(MPI_STATUS_SIZE)

          call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
          call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )

          if(.not.allocated(Tkin_col_p0)) &
     &       allocate(Tkin_col_p0(ND)) 

          h2_PSI_local=0.0d0

!          write(*,*) myid, allocated(Tkin_col_p0)
!          if(myid.eq.0) write(*,*) Tkin_2D(1,1)
!          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!          stop  

!       FOR all PEs
          DO i=1, ND
            h2=0d0 
        
            if(dim_mctdhb.eq.2.and.Tkin_prec.eqv..TRUE.) then
              Tkin_col_p0=0d0
                if(myid.eq.0) then 
                   Tkin_col_p0(:)=Tkin_2D(:,i)
                end if
                CALL MPI_BCAST(Tkin_col_p0,  &
     &               ND,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD, &
     &               ierr)
            end if
             
            
            CALL h2_help(h2,i)




             DO q=1,Morb
!              DO i=1,ND
                DO j=1,ND
      
!                   h2_PSI(i,q)=h2_PSI(i,q)+h2(i,j)*PSI(j,q) 
                   h2_PSI_local(i,q)=h2_PSI_local(i,q)+h2(j)*PSI(j,q) 
           
                END DO
             END DO
          END DO



       END SUBROUTINE get_h2_PSI

!==========================================================================================

       SUBROUTINE h2_help(h2_out,row)

          USE CI_ALL
          USE Parallel_CI 
          USE W_INTERPARTICLE !ADD BY STR
          USE DVR_ALL
          USE rR_hW
          USE LR_RAPHA 
           
          IMPLICIT NONE         

          INTEGER :: i, j, k, q, s, ll, mm, kk, row, pos, counter, x_pos, y_pos, z_pos 
          INTEGER  :: i1,i2,j1,j2,q1,q2

!          COMPLEX*16, DIMENSION(dimL_orb/2), INTENT(OUT) :: h2_out
          COMPLEX*16, DIMENSION(dimL_orb/2), INTENT(OUT) :: h2_out
          COMPLEX*16, DIMENSION(ND,ND) :: ham
           
          COMPLEX*16, DIMENSION(ND) :: Tkin_row,Tkin_col
 
!     Non-interacting
         h2_out=0.0d0;


         pos=mod(row,ND)
         if(pos.eq.0) pos=ND
!         write(*,*) "pos",pos


!          write(*,*) "DVR",Time_DVRMETHODX
!          IF(Time_DVRMETHODX.eq.4) then 
!            call T_FFT_GENERAL(Tkin_col,pos)
!            Tkin_col=Conjg(Tkin_col)
!            Tkin_col(pos)=Tkin_col(pos)+VTRAP_EXT(pos)
!          end if
!          ELSE 
!          Tkin_row=Op_X
!          ENDIF

         

         counter=1
         Tkin_col=0d0

         if(dim_mctdhb.eq.1) then 
            IF(Time_DVRMETHODX==4) THEN ! FFT - explicit construction Tkin matrix
              if(Tkin_prec.eqv..TRUE.) then
                Tkin_col(:)=Tkin_1D(:,pos) 
              else
                call T_FFT_GENERAL(Tkin_col,pos)
              end if
            ELSE  ! HO-DVR CASE
!           do i=1, NDX
!            Tkin_col(counter)=Tkin_col(counter)+Op_x(i,pos)
!            counter=counter+1 
!           end do
                Tkin_col(:)=Op_x(:,pos)
            END IF 
         end if


         if(dim_mctdhb.eq.2) then 
          if(Tkin_prec.eqv..FALSE.) then

           if(Time_DVRMETHODX==1) then
            ham=0d0
            if(mod(row,NDX).eq.0) then
              y_pos=row/NDX
            else
              y_pos=(row/NDX)+1 
            end if

            x_pos=mod(row,NDX)
            if(x_pos.eq.0) x_pos=NDX 

          Do j1=1,NDY
!          Do j2=1,NDY
          j2=y_pos 
          IF(j1.eq.j2) then
          q1=1
          else
          q1=0
          endIF
             Do i1=1,NDX
!             Do i2=1,NDX
             i2=x_pos
             IF(i1.eq.i2) then
             q2=1
             else
             q2=0
             endIF
             I=NDX*(j1-1)+i1
!             J=NDX*(j2-1)+i2
             Tkin_col(I)=q1*Op_X(i1,i2)+q2*Op_Y(j1,j2)
!             ham(I,J)=q1*Op_X(i1,i2)+q2*Op_Y(j1,j2)
!         IF(I.eq.J) ham(I,I)=ham(I,I)+Vext(I)
!      write(6,'(i2,i2,a10,4i2)') I,J,"i1i2_j1j2",q1*i1,q1*i2,q2*j1,q2*j2
!            endDo
             endDo
!          endDo
          endDo
           else
             if(Time_DVRMETHODX==4) &
     &                    call T_FFT_GENERAL(Tkin_col,row)
           end if
          else 
!ORIG            Tkin_col(:)=Tkin_2D(:,pos) 
            Tkin_col=Tkin_col_p0 
          end if 
         end if


!         if(dim_mctdhb.eq.2) then 
!          do j=1, NDY
!           do i=1, NDX

!            if(mod(row,NDX).eq.0) then
!              y_pos=row/NDX
!            else
!              y_pos=(row/NDX)+1 
!            end if

!            x_pos=mod(row,NDX)
!            if(x_pos.eq.0) x_pos=NDX 


!            Tkin_col(counter)=Tkin_col(counter)+Op_x(i,x_pos)+Op_y(j,y_pos)
!            if (counter.eq.pos) Tkin_col(counter)=Tkin_col(counter)+VTRAP_EXT(counter) 
!            counter=counter+1 
!           end do
!          end do
!         end if

        if(dim_mctdhb.eq.3) then 
         do i=1, NDX
          do j=1, NDY
           do k=1, NDZ
            Tkin_col(counter)=Tkin_col(counter)+Op_x(pos,i)+Op_y(pos,j)+Op_z(pos,k)
!            if (counter.eq.pos) Tkin_col(counter)=Tkin_col(counter)+VTRAP_EXT(counter) 
            counter=counter+1 
           end do
          end do
         end do 
        end if 

       Tkin_col(pos)=Tkin_col(pos)+VTRAP_EXT(pos)

         Tkin_row=Conjg(Tkin_col) 

!              Tkin=Op_X !STR October 2012 modifications
!              DO j=1,ND
!              Tkin(j,j)=Tkin(j,j)+REAL(VTRAP_EXT(j))
!              END DO

         DO k=1,Morb
             h2_out((k-1)*ND+1:k*ND)=Tkin_row
!             h2_out=Tkin_row
         END DO

       END SUBROUTINE h2_help


!=============================================================================================================



       subroutine Produce_Cij(VIN,VOUT,i,j,k,l,Nc)
       USE CI_All
       USE Parallel_CI
       USE omp_lib
       USE CI_prod
       implicit NONE
!=================== MPI 
       INCLUDE 'mpif.h'
       INTEGER ::  ierr,MYID,IPRC,numprocs
!==========================================================
       COMPLEX*16, INTENT(IN) :: VIN(Nc)
       COMPLEX*16, INTENT(OUT) :: VOUT(Nc)
!==========================================================
       integer :: nn,Iorb,Jorb,ihelp,ihelphelp,jhelp,jhelphelp,khelp,lhelp,conj,suml,sumr
       integer :: maxsil,NL,Iter_Total,Nterms
       INTEGER :: II, KK, PP, cI, cJ, cK, cL, I_current_term
       INTEGER :: myII, myPP, mycI, mycJ, mycK, mycL
       real*8 :: xbprefac
       INTEGER, INTENT(IN) :: i, j, k, l, Nc

!=====1-body===============================================
       IF(k==0.OR.l==0) THEN

!         discriminate ij and ji

          conj=0
          ihelp=i; jhelp=j;
          IF(i>j) THEN
             ihelp=j; jhelp=i;
             conj=1
          END IF

!         create index for mapping

!          II=1
!          DO Iorb=1,Morb
!             DO Jorb=Iorb,Morb
!                IF(ihelp==Iorb.and.jhelp==Jorb) GOTO 119
!                II=II+1
!             EndDO      
!          EndDO    

       !   write(*,*) MaxTrm1b,TERM_REQ_1B


!         Find index II from ihelp and jhelp
                VOUT=0d0 !Added by str
!              write(6,*)"Im in Produce _CI"
          DO II=1,MaxTrm1b
             I_current_term= TERM_REQ_1B(II)
             PP=TERM_INDEX_1B(I_current_term)
             cJ= INT(PP/100); cI=PP-cJ*100; 
             IF(cI==ihelp.AND.cJ==jhelp) THEN
              myII=II; MyPP=PP; mycI=cI; mycJ=cJ;
             END IF
          END DO
        II=myII; PP=MyPP; cI=mycI; cJ=mycJ;
!              write(6,*)"Im in Produce II, PP ",II,PP
!              write(6,*)"Im in Produce I,J",cI,CJ
!              write(6,*)"Prefactors CI 1 1",Prefactors_1b(:,II)
!              write(6,*)"Im in Produce conj",conj
!STRSTR       WRITE(6,*) 'II, ij', II, 'cI:', cI, 'cJ:', cJ, 'I_curr', I_current_term



!         loop over configurations

          IF(conj==0) THEN
             Do KK=1,Nc 
                PP=Ind_CI_1b_Arnoldi(KK,II)
                xbprefac=Prefactors_1b_Arnoldi(KK,II)
                VOUT(PP)=VOUT(PP)+VIN(KK)* DSQRT(xbprefac) ! STR corrected
!                write(6,*)" qq", KK, "v",PP," times ",xbprefac
!                write(6,*)"Vin", VIN(KK), " Vout", VOUT(PP)
!                pause
             EndDO       
          ELSE 
             Do KK=1,Nc
                PP=Ind_CI_1b_Arnoldi(KK,II)
                xbprefac=Prefactors_1b_Arnoldi(KK,II)
                VOUT(KK)=VOUT(KK)+VIN(PP)* DSQRT(xbprefac)
             EndDO     
          END IF
!              write(6,*)"INd CI 1 1",Ind_CI_1b(:,1)
!              write(6,*)"VIN",VIN(:)
!              write(6,*)"VOUT",VOUT(:)
!              write(6,*)"Prefactors CI 1 1",Prefactors_1b(:,3)
!=====2-body===============================================
       ELSE
     
!         discriminate ijkl: ij,ji and kl,lk

          ihelp=i; jhelp=j;
          IF(i>j) THEN
             ihelp=j; jhelp=i;
          END IF
          khelp=k; lhelp=l;
          IF(k>l) THEN
             khelp=l; lhelp=k;
          END IF

!         discriminate ijkl: ij and kl

          suml=i+j; sumr=k+l;
          conj=0
          IF(ihelp>khelp ) THEN
             ihelphelp=ihelp; jhelphelp=jhelp
             ihelp=khelp; jhelp=lhelp; khelp=ihelphelp; lhelp=jhelphelp;
             conj=1
          END IF
          IF(ihelp==khelp .AND. jhelp>lhelp ) THEN
             ihelphelp=ihelp; jhelphelp=jhelp
             ihelp=khelp; jhelp=lhelp; khelp=ihelphelp; lhelp=jhelphelp;
             conj=1
          END IF

!         create index for mapping
          DO II=1,MaxTrm2b
             I_current_term= TERM_REQ_2B(II)
             PP=TERM_INDEX_2B(I_current_term)
             cL= INT(PP/1000000)
             cK= INT((PP-cL*1000000)/10000)
             cJ= INT((PP-cL*1000000-cK*10000)/100)
             cI= PP-cL*1000000-cK*10000-cJ*100
             IF(cI==ihelp.AND.cJ==jhelp.AND.cK==khelp.AND.cL==lhelp) THEN
              myII=II; MyPP=PP; mycI=cI; mycJ=cJ; mycK=cK; mycL=cL
             END IF
          END DO
        II=myII; PP=MyPP; cI=mycI; cJ=mycJ; cK=mycK; cL=mycL
!STR STR   WRITE(6,*) 'II',II,'PP',PP, 'cI', cI, 'cJ', cJ, 'cK', cK, 'cL', cL, 'conj', conj

!         loop over configurations

          IF(conj==0) THEN
             Do KK=1,Nc 
                PP=Ind_CI_2b_Arnoldi(KK,II)
                xbprefac=Prefactors_2b_Arnoldi(KK,II)
                VOUT(PP)=VOUT(PP)+VIN(KK)* DSQRT(xbprefac) !STR
             EndDO       
          ELSE 
             Do KK=1,Nc
                PP=Ind_CI_2b_Arnoldi(KK,II)
                xbprefac=Prefactors_2b_Arnoldi(KK,II)
                VOUT(KK)=VOUT(KK)+VIN(PP)* DSQRT(xbprefac) ! STR
             EndDO     
          END IF

       END IF



       end subroutine Produce_Cij



!========================================================================================

       subroutine Get_Full_Rijkl(Rho_ijkl)
       USE   PASS_ARG
       USE   SHARED_DIMS
       USE   rR_hW
       USE   CI_All
!       USE   PROP_MB
       USE   W_INTERPARTICLE
       USE   DVR_ALL
       USE   Parallel_Orb
       integer :: i,P,cK,cJ,cL,cI
       complex*16 :: rho_jk(Morb,Morb),check
       complex*16 :: rho_ijkl(Morb,Morb,Morb,Morb)
!=====================================================================
!      write(6,*) "IN=",MaxTrm2b
         Do I=1,MaxTrm2b
          P=TERM_INDEX_2B(I)
!================ Unpack cI cJ cK cL from P
          cL= INT(P/1000000)
          cK= INT((P-cL*1000000)/10000)
          cJ= INT((P-cL*1000000-cK*10000)/100)
          cI= P-cL*1000000-cK*10000-cJ*100
        Rho_ijkl(cI,cJ,cK,cL)=ZRIJKL(I) 
        Rho_ijkl(cI,cJ,cL,cK)=ZRIJKL(I) 
        Rho_ijkl(cJ,cI,cK,cL)=ZRIJKL(I) 
        Rho_ijkl(cJ,cI,cL,cK)=ZRIJKL(I) 
      IF(cI.ne.cK) then
       Rho_ijkl(cK,cL,cI,cJ)=Conjg(ZRIJKL(I))
       Rho_ijkl(cK,cL,cJ,cI)=Conjg(ZRIJKL(I))
       Rho_ijkl(cL,cK,cI,cJ)=Conjg(ZRIJKL(I))
       Rho_ijkl(cL,cK,cJ,cI)=Conjg(ZRIJKL(I))
      else
      IF(cJ.ne.cL) then
        Rho_ijkl(cK,cL,cI,cJ)=Conjg(ZRIJKL(I))
        Rho_ijkl(cK,cL,cJ,cI)=Conjg(ZRIJKL(I))
        Rho_ijkl(cL,cK,cI,cJ)=Conjg(ZRIJKL(I))
        Rho_ijkl(cL,cK,cJ,cI)=Conjg(ZRIJKL(I))
      endif
      endif
      EndDo

      end subroutine Get_Full_Rijkl



       subroutine Get_Full_Rij(Rho_JK)
       USE   PASS_ARG
       USE   SHARED_DIMS
       USE   rR_hW
       USE   CI_All
!       USE   PROP_MB
       USE   W_INTERPARTICLE
       USE   DVR_ALL
       USE   Parallel_Orb
       integer :: i,P,cK,cJ,icntr,Iorb,Jorb
       complex*16 :: rho_jk(Morb,Morb)
!=====================================================================
!      write(6,*) "IN=",MaxTrm1b
!         ZRIJ=ZRIG+(3.0,5.0)
        Do I=1,MaxTrm1b
        P=TERM_INDEX_1B(I)
!================ Unpack cI cJ cK cL from P
        cK= INT(P/100)
        cJ= P-cK*100
        Rho_JK(cJ,cK)=ZRIJ(I) 
        IF(cK.ne.cJ) Rho_JK(cK,cJ)=Conjg(ZRIJ(I))
!        write(6,*) cJ,cK,"NEW Rho",Rho_JK(cJ,cK)
        EnDdo
!        write(6,*) "NEW", Rho_JK
        
!========== BOTH ARE CORRECT !!!!!!!!!!!!!!
        icntr=1
        DO Iorb=1,Morb
        DO Jorb=Iorb,Morb

        Rho_JK(Iorb,Jorb)=ZRIJ(icntr)
        IF(Jorb.ne.Iorb) Rho_JK(Jorb,Iorb)=Conjg(ZRIJ(icntr)) ! Correct
        icntr=icntr+1
        EndDO
        EndDO
!        write(6,*) "OLD", Rho_JK
!        write(6,*) "OLD Rho",Rho_JK
!          pause
      end subroutine Get_Full_Rij


      subroutine Get_Full_Wijkl(W_ijkl)
       USE   PASS_ARG
       USE   SHARED_DIMS
       USE   rR_hW
       USE   CI_All
!       USE   PROP_MB
       USE   W_INTERPARTICLE
       USE   DVR_ALL
       USE   Parallel_Orb
       integer :: i,q,s,r,P,cK,cJ,cL,cI
       complex*16 :: W_ijkl(Morb,Morb,Morb,Morb)
!=====================================================================
         Do I=1,MaxTrm2b
          P=TERM_INDEX_2B(I)
!================ Unpack cI cJ cK cL from P
          cL= INT(P/1000000)
          cK= INT((P-cL*1000000)/10000)
          cJ= INT((P-cL*1000000-cK*10000)/100)
          cI= P-cL*1000000-cK*10000-cJ*100
        IF(mod(cL+cK+cJ+cI,2)==0) THEN
          W_ijkl(cI,cJ,cK,cL)=WIJKL(I) 
          W_ijkl(cI,cJ,cL,cK)=WIJKL(I) 
          W_ijkl(cJ,cI,cK,cL)=WIJKL(I) 
          W_ijkl(cJ,cI,cL,cK)=WIJKL(I) 
          IF(cI.ne.cK) then
             W_ijkl(cK,cL,cI,cJ)=Conjg(WIJKL(I))
             W_ijkl(cK,cL,cJ,cI)=Conjg(WIJKL(I))
             W_ijkl(cL,cK,cI,cJ)=Conjg(WIJKL(I))
             W_ijkl(cL,cK,cJ,cI)=Conjg(WIJKL(I))
          else
          IF(cJ.ne.cL) then
             W_ijkl(cK,cL,cI,cJ)=Conjg(WIJKL(I))
             W_ijkl(cK,cL,cJ,cI)=Conjg(WIJKL(I))
             W_ijkl(cL,cK,cI,cJ)=Conjg(WIJKL(I))
             W_ijkl(cL,cK,cJ,cI)=Conjg(WIJKL(I))
          endif
          endif
        ELSE
          W_ijkl(cI,cJ,cK,cL)=WIJKL(I)!*0.5d0 
          W_ijkl(cI,cJ,cL,cK)=WIJKL(I)!*0.5d0 
          W_ijkl(cJ,cI,cK,cL)=WIJKL(I)!*0.5d0  
          W_ijkl(cJ,cI,cL,cK)=WIJKL(I)!*0.5d0  
          IF(cI.ne.cK) then
             W_ijkl(cK,cL,cI,cJ)=Conjg(WIJKL(I))!*0.5d0  
             W_ijkl(cK,cL,cJ,cI)=Conjg(WIJKL(I))!*0.5d0  
             W_ijkl(cL,cK,cI,cJ)=Conjg(WIJKL(I))!*0.5d0  
             W_ijkl(cL,cK,cJ,cI)=Conjg(WIJKL(I))!*0.5d0  
!       WRITE(6,*) 'h11check', cI,cJ,cK,cL
          else
          IF(cJ.ne.cL) then
             W_ijkl(cK,cL,cI,cJ)=Conjg(WIJKL(I))!*0.5d0  
             W_ijkl(cK,cL,cJ,cI)=Conjg(WIJKL(I))!*0.5d0  
             W_ijkl(cL,cK,cI,cJ)=Conjg(WIJKL(I))!*0.5d0  
             W_ijkl(cL,cK,cJ,cI)=Conjg(WIJKL(I))!*0.5d0  
          endif
          endif
        END IF
      EndDo

      end subroutine Get_Full_Wijkl



      subroutine Get_Full_Hij(H_ij)
       USE   PASS_ARG
       USE   SHARED_DIMS
       USE   rR_hW
       USE   CI_All
!       USE   PROP_MB
       USE   W_INTERPARTICLE
       USE   DVR_ALL
       USE   Parallel_Orb
       integer :: icntr, Iorb, Jorb
       complex*16, DIMENSION(Morb,Morb) :: H_ij
           icntr=1  
        DO Iorb=1,Morb
           DO Jorb=Iorb,Morb
              H_ij(Iorb,Jorb)=HIJ(icntr)       
              H_ij(Jorb,Iorb)=Conjg(H_ij(Iorb,Jorb))        
              icntr=icntr+1                                     
           EndDO   
        EndDO    
      end subroutine Get_Full_Hij



END MODULE ARNOLDI_MOD

