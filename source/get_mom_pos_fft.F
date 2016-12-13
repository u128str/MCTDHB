! Subroutined needed to work PROPERLY with FFT in Kinetics and in Wsl evaluations ! THanks to Kaspar and Axel   
      subroutine getpos(xi,xf,gdim,pos)
      USE SHARED_DIMS, ONLY : PI
       !---{{{ returns pos grid
       implicit none

       integer,intent(IN)              :: gdim
       double precision,intent(IN)     :: xi,xf
       double precision,intent(OUT)    :: pos(gdim)

       double precision                :: dx
       integer                         :: i

       dx             =     (xf-xi)/(1.d0*gdim)
       do i=1,gdim
         pos(i) = xi+(i-1)*dx
       enddo
       !---}}}
       end subroutine getpos

      subroutine getmom(xi,xf,gdim,mom)
      USE SHARED_DIMS, ONLY : PI
       !---{{{ returns momentum grid in fft order
       implicit none

       integer,intent(IN)              :: gdim
       double precision,intent(IN)     :: xi,xf
       double precision,intent(OUT)    :: mom(gdim)

       double precision                :: dp
       integer                         :: i

       dp             =     2.d0*pi/(xf-xi)
       do i=1,gdim
         mom(i) = (i-1)*dp
       enddo

       do i=gdim/2+2,gdim !put momenta in FFT
         mom(i)=1.d0*(i-gdim-1)*dp
       end do

       !---}}}
       end subroutine getmom

