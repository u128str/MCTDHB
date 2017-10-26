       subroutine computevtilde_FFTWIMEST(vtilde,wtilde)

       use SHARED_DIMS 
       use DVR_ALL
       use W_INTERPARTICLE
!       use MKL_DFTI
       use omp_lib

      implicit none

       double complex                  :: vtilde(ndx*ndy*ndz)
       double complex                  :: wtilde(ndx*ndy*ndz)
       double precision                :: x,y,z,sigma,dx,dy,dz,xall,zall
       double precision                :: yall  
       integer                         :: m,i,j,k

       include 'fftw3.f'
       integer*8 plan,plany
       integer iret,nthr
!       COMPLEX*16 :: dummy(NDY),dummyX(NDX)
       INTEGER ::  ierr
       INTEGER :: FFTW_USE_WISDOM,isuccess
!       COMPLEX*16,   DIMENSION(NDX*NDY*NDZ) :: psi
!---{{{
c=========================================================
!#if FFTW
!      write(*,*)'THE IMEST ALGORITHM IN V3.0 IS IMPLEMENTED WITH FFTW!'
      nthr = OMP_GET_MAX_THREADS()
!      write(6,*)"IN THE COMPUTATION OMP_GET_MAX_THREADS is=", nthr
!       stop
!#endif
cccccccccccccccccc
cccccc define grid spacing
       dx=0d0
       dy=0d0
       dz=0d0
      dx=ort_x(2)-ort_x(1)     
      IF(DIM_MCTDHB.ge.2) dy=ort_y(2)-ort_y(1)
      IF(DIM_MCTDHB.ge.3) dz=ort_z(2)-ort_z(1)
CCCCCCCCCCCCCCCCCC
CCCCCCC definition of spatial grid

      xall=ort_x(NDX)-ort_x(1)+dx     
      yall=ort_y(NDY)-ort_y(1)+dy
      zall=ort_z(NDZ)-ort_z(1)+dz     

ccccccc  do the mkl fft voodoo
      if (DIM_MCTDHB.eq.1) then 
         if(MOD(NDX,2)/=0) stop "For IMEST works only even NDX"

      call dfftw_plan_dft_1d(plan,NDX,wtilde,wtilde,
     .                       FFTW_FORWARD,
     .                       FFTW_ESTIMATE)
      call dfftw_execute_dft(plan,wtilde,wtilde)
      call dfftw_destroy_plan(plan)
      else if (DIM_MCTDHB.eq.2) then 
         if(MOD(NDX,2)/=0) stop "For IMEST works only even NDX"
         if(MOD(NDY,2)/=0) stop "For IMEST works only even NDY"
      call dfftw_plan_dft_2d(plan,NDX,NDY,wtilde,wtilde,
     .                FFTW_FORWARD,FFTW_ESTIMATE,fftw_unaligned)
      call dfftw_execute_dft(plan,wtilde,wtilde)
      call dfftw_destroy_plan(plan)
      else if (DIM_MCTDHB.eq.3) then 
         if(MOD(NDX,2)/=0) stop "For IMEST works only even NDX"
         if(MOD(NDY,2)/=0) stop "For IMEST works only even NDY"
         if(MOD(NDZ,2)/=0) stop "For IMEST works only even NDZ"
      call dfftw_plan_dft_3d(plan,NDX,NDY,NDZ,wtilde,wtilde,
     .                FFTW_FORWARD,FFTW_ESTIMATE,fftw_unaligned)
      call dfftw_execute_dft(plan,wtilde,wtilde)
      call dfftw_destroy_plan(plan)
      endif 
       if (DIM_MCTDHB.eq.1) then 
!--- {{{ 1d 
         if(ndx*ndy*ndz.ne.ndx) then
           write(6,*)"this is not good" 
           stop
         endif
      do i=1,ndx
         vtilde(i) = wtilde(i)*(-1.d0)**(i-1)
      enddo

!---}}}
       else if (DIM_MCTDHB.eq.2) then
!---{{{ 2d 
       if(ndx*ndy*ndz.ne.ndx*ndy) then
           write(6,*)"this is not good" 
         stop
       endif

         do m=1,NDX*NDY*NDZ
             call get_ijk_from_m(m,ndx,ndy,i,j,k)
           vtilde(m) =  (-1.d0)**(i-1+j-1)*wtilde(m) 
         enddo


!---}}}
       else if (DIM_MCTDHB.eq.3) then
!---{{{
c       stop

         do m=1,NDX*NDY*NDZ
           call get_ijk_from_m(m,ndx,ndy,i,j,k)
           vtilde(m) =(-1.d0)**(i-1 + j-1 + k-1) * wtilde(m) 
         enddo

!---}}}
       else
         write(6,*)"more than 3d?"
         stop
       endif

!---}}}
       end subroutine computevtilde_FFTWIMEST

       subroutine get_ijk_from_m(m,gdimx,gdimy,i,j,k)
!the routine returns i,j,k given m---{{{
!of a vector that stores the ijk-th element in
!the position m=i+gdimx(j-1)+gdimy*gdimx*(k-1)
![(ijk) start with (111)]
       implicit none
       integer     :: m,gdimx,gdimy
       integer     :: i,j,k
       integer                 :: n

         n = m
         k = (n-1)/(gdimx*gdimy) + 1 
         n =  n-(k-1)*gdimx*gdimy
         j = (n-1)/gdimx + 1
         n =  n-(j-1)*gdimx
         i =  n
       end subroutine get_ijk_from_m 

!---}}}

