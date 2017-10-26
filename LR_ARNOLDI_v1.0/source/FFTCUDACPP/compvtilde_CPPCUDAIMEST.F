       subroutine computevtilde_CPPCUDAIMEST(vtilde,wtilde)

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

!       include 'fftw3.f'
       integer*8 plan,plany
       integer iret,nthr
       COMPLEX*16 :: ortkx(NDX)
       INTEGER ::  ierr
       INTEGER :: FFTW_USE_WISDOM,isuccess
!       COMPLEX*16,   DIMENSION(NDX*NDY*NDZ) :: psi
!---{{{
!=========================================================
      write(6,*)"THE IMEST ALGORITHM IN V3.0 WITH CUDACPP!",DIM_MCTDHB
!      nthr = OMP_GET_MAX_THREADS()
!      write(6,*)"IN THE COMPUTATION OMP_GET_MAX_THREADS is=", nthr
      if (DIM_MCTDHB.eq.1) then 
         if(MOD(NDX,2)/=0) stop "For IMEST works only even NDX"
!          call kernel_imestfft(wtilde, NDX,NDY,NDZ, DIM_MCTDHB)
           
      write(6,*)"BEFORE kernel_imestfft_T!"
           call kernel_imestfft_T(wtilde, NDX,NDY,NDZ, 1 )


      write(6,*)"DONE THE IMEST ALGORITHM IN V3.0 WITH CUDACPP!"
      else if (DIM_MCTDHB.eq.2) then 
         if(MOD(NDX,2)/=0) stop "For IMEST works only even NDX"
         if(MOD(NDY,2)/=0) stop "For IMEST works only even NDY"

          call kernel_imestfft_T(wtilde, NDX, NDY ,NDZ, 2)

      else if (DIM_MCTDHB.eq.3) then 
         if(MOD(NDX,2)/=0) stop "For IMEST works only even NDX"
         if(MOD(NDY,2)/=0) stop "For IMEST works only even NDY"
         if(MOD(NDZ,2)/=0) stop "For IMEST works only even NDZ"
        
         call kernel_imestfft_T(wtilde, NDX, NDY ,NDZ, 3)

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
! c       stop

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
       end subroutine computevtilde_CPPCUDAIMEST

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

