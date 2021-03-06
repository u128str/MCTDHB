       subroutine computevtilde_MKLIMEST(vtilde,wtilde)
       use SHARED_DIMS 
       use DVR_ALL
       use W_INTERPARTICLE
      Use MKL_DFTI
!---{{{

       implicit none

       double complex                  :: vtilde(ndx*ndy*ndz)
       double complex                  :: wtilde(ndx*ndy*ndz)
       double precision                :: x,y,z,sigma,dx,dy,dz,xall,zall
       double precision                :: yall  
       integer                         :: m,i,j,k

      integer*8 plan
cccccc MKL STUFF
      INTEGER :: STRIDE(2),INFO
      INTEGER :: L(2),G(3)
      type(DFTI_DESCRIPTOR), POINTER :: Desc_Handle_Dim1
      Integer :: Status
c=========================================================
#if FFTW
      write(*,*) 'THE IMEST ALGORITHM IS NOT IMPLEMENTED WITH FFTW!!!'
#endif
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
         Status = DftiCreateDescriptor(Desc_Handle_Dim1, DFTI_DOUBLE,
     &           DFTI_COMPLEX, 1, NDX )
         Status = DftiCommitDescriptor(Desc_Handle_Dim1)
         Status = DftiComputeForward(Desc_Handle_Dim1,wtilde)
      else if (DIM_MCTDHB.eq.2) then 
         if(MOD(NDX,2)/=0) stop "For IMEST works only even NDX"
         if(MOD(NDY,2)/=0) stop "For IMEST works only even NDY"
         L(1)=NDX
         L(2)=NDY
         Status = DftiCreateDescriptor(Desc_Handle_Dim1, DFTI_DOUBLE,
     &        DFTI_COMPLEX, 2, L )
         Status = DftiCommitDescriptor( Desc_Handle_Dim1) 
         Status = DftiComputeForward(Desc_Handle_Dim1,wtilde)
      else if (DIM_MCTDHB.eq.3) then 
         if(MOD(NDX,2)/=0) stop "For IMEST works only even NDX"
         if(MOD(NDY,2)/=0) stop "For IMEST works only even NDY"
         if(MOD(NDZ,2)/=0) stop "For IMEST works only even NDZ"
         G(1)=NDX
         G(2)=NDY
         G(3)=NDZ
         Status = DftiCreateDescriptor(Desc_Handle_Dim1, DFTI_DOUBLE,
     &        DFTI_COMPLEX, 3, G )
         Status = DftiCommitDescriptor( Desc_Handle_Dim1) 
         Status = DftiComputeForward(Desc_Handle_Dim1,wtilde)
      endif 
      Status = DftiFreeDescriptor(Desc_Handle_Dim1)

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
       end subroutine computevtilde_MKLIMEST

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

