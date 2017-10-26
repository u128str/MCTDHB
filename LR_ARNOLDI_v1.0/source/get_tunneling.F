c#  This routine calculates the number of particles
c#  inside the center and the outer ring for the 2D crater potential        


      SUBROUTINE  get_tunneling(Phi_Vec,t,t_next,t_print)
      USE   SHARED_DIMS
      USE   W_INTERPARTICLE
      USE   DVR_ALL
      USE   rR_hW
      USE USR_PAR
      IMPLICIT NONE
      integer                    ::      i,j,k,ind,status_open
      integer                    ::      ierr,Iorb,Jorb
      real*8                     ::      t,r,x,y,z,t_print,t_next
      real*8                     ::      occ_in,occ_out,occ_nirv
      COMPLEX*16,   DIMENSION(NDX*NDY*NDZ,Morb) :: phi_vec
      COMPLEX*16 :: rho
      logical :: exists


c!=========================================================================

      INQUIRE(FILE="tunneling.dat", EXIST=exists)


      if(exists) then
          open(unit=123, file='tunneling.dat', status='old',
     .position='append', action='readwrite', iostat=status_open)
          if(status_open /= 0) write(*,*)'tunneling.dat error'
      else
          open(unit=123, file='tunneling.dat',
     .status='new', action='readwrite', iostat=status_open)
           if(status_open /= 0) write(*,*)'tunneling.dat error!'
      end if




      occ_in=0.0
      occ_out=0.0
      occ_nirv=0.0

      ind=1
      do K=1, NDZ
      do J=1, NDY
      do I=1, NDX       

      rho=cmplx(0.0,0.0,kind=8)
        do Iorb=1,Morb
         do Jorb=1, MOrb

            rho=rho+conjg(Phi_vec(ind,Iorb))*AllZRIJ(Iorb,Jorb)
     .*Phi_vec(ind,Jorb)
       
         end do
        end do


        r=DSQRT(ort_X(i)**2+ort_Y(j)**2)
        x=ort_X(i)
        y=ort_y(j)
        z=ort_z(k)


        if(sqrt(x**2+y**2) <= Rb) then
            occ_in=occ_in+real(rho)
        else
           if(sqrt(x**2+y**2)<=Rc) then
              occ_out=occ_out+real(rho)
           else 
              occ_nirv=occ_nirv+real(rho)
           end if
        end if
 
 
         ind=ind+1


      end do
      end do
      end do
     


      write(*,*) "Occ IN: ",occ_in
      write(*,*) "Occ OUT: ",occ_out



!      if(abs(t_next-t_print).le.1.0d-8.or.mod(t,0.1).eq.0.0) then
!      if(ABS(t-t_print).le.1.0d-8) then

      INQUIRE(FILE="tunneling.dat", EXIST=exists)
 
 
      if(exists) then
           open(unit=123, file='tunneling.dat', status='old',
     .position='append', action='readwrite', iostat=status_open)
          if(status_open /= 0) write(*,*)'tunneling.dat error'
      else
          open(unit=123, file='tunneling.dat',
     .status='new', action='readwrite', iostat=status_open)
            if(status_open /= 0) write(*,*)'tunneling.dat error!'
      end if

      write(123,526) t,occ_in,occ_out,occ_nirv     
 
      close(123)    

!      end if
526    format(4F26.16)

       END SUBROUTINE get_tunneling

