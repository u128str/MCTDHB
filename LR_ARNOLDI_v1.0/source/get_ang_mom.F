c#  This routine gives you the elements L_xyz(Iorb,Jorb)
c#  to compute total angular momentum L=Sum_i,j rho_ij*L_xyz(i,j)        


      SUBROUTINE  get_ang_mom(Phi_Vec,LPhi,dir,t,t_next,t_print)
      USE   SHARED_DIMS
      USE   W_INTERPARTICLE
      USE   DVR_ALL
      USE   rR_hW
      USE USR_PAR
      IMPLICIT NONE
      integer                    ::      i,j,k,ind,status_open
      integer                    ::      ierr,dir,Iorb,Jorb,jj,k1,l1
      real*8                     ::      x,y,z,t,t_next,t_print
      complex*16 compt,ang_mom,Z1
      COMPLEX*16,   DIMENSION(NDX*NDY*NDZ,Morb) :: phi_vec
      COMPLEX*16,   DIMENSION(Morb,Morb) :: LPHI
      COMPLEX*16,   DIMENSION(Morb) :: ang_mom_orb
      COMPLEX*16,   DIMENSION(Morb) :: ang_mom_orb_2
      COMPLEX*16,   DIMENSION(Morb) :: ang_mom_orb_work
      logical :: exists,exists2,exists3,exists4,exists5

!      complex*16, allocatable :: Phi_Vec(:,:)
!      complex*16, allocatable :: LPhi(:,:)
      complex*16, allocatable :: temp(:,:)
      complex*16, allocatable :: v_temp(:,:)
      complex*16, allocatable :: deriv1(:,:)
      complex*16, allocatable :: deriv2(:,:)
      complex*16, allocatable :: v_deriv1(:,:)
      complex*16, allocatable :: v_deriv2(:,:)
      complex*16, allocatable :: v(:,:)

c!=========================================================================


        allocate(temp(NDX*NDY*NDZ,Morb),stat=ierr)
        if(ierr /= 0) write(*,*)"allocation error for temp"

        allocate(v_temp(NDX*NDY*NDZ,Morb),stat=ierr)
        if(ierr /= 0) write(*,*)"allocation error for v_temp"

        allocate(deriv1(NDX*NDY*NDZ,Morb),stat=ierr)
        if(ierr /= 0) write(*,*)"allocation error for deriv1"

        allocate(deriv2(NDX*NDY*NDZ,Morb),stat=ierr)
        if(ierr /= 0) write(*,*)"allocation error for deriv2"

        allocate(v_deriv1(NDX*NDY*NDZ,Morb),stat=ierr)
        if(ierr /= 0) write(*,*)"allocation error for v_deriv1"

        allocate(v_deriv2(NDX*NDY*NDZ,Morb),stat=ierr)
        if(ierr /= 0) write(*,*)"allocation error for v_deriv2"

        allocate(v(NDX*NDY*NDZ,Morb),stat=ierr)
        if(ierr /= 0) write(*,*)"allocation error for v"
         
         LPhi=cmplx(0.0,0.0,kind=8) 
         temp=cmplx(0.0,0.0,kind=8) 
         v_temp=cmplx(0.0,0.0,kind=8) 
         ang_mom=cmplx(0.0,0.0,kind=8) 
         ang_mom_orb=cmplx(0.0,0.0,kind=8) 
         ang_mom_orb_2=cmplx(0.0,0.0,kind=8) 
         ang_mom_orb_work=cmplx(0.0,0.0,kind=8) 
         v=cmplx(0.0,0.0,kind=8) 
          


!==== NATURAL ORBITALS
      do ind=1, NDX*NDY*NDZ
      Do k1=1,Morb
          Do l1=1,Morb
!================Natural orbitals are
              v(ind,k1)=v(ind,k1)+Conjg(NatVec(l1,k1))*phi_vec(ind,l1) ! Natural orbitals are ORG
! They are computed according to: \psi^{NO}_i= \sum_j U^{-1}_ij \phi_j = \sum_j U^*_ji \phi_j
! rho(x,x')= U U^-1 rho U U^-1=>(\phi_1^*,\phi_2^*,...\phi_M^*) U rho_d U^-1 ((\phi_1^*,\phi_2^*,...\phi_M^*)^T
          END DO
!       write(*,*) Nocc(k1)
      END DO
      end do

c=========  Get D/Dx,D/Dy,D/Dz
         do Iorb=1, Morb

            deriv1(:,Iorb) = phi_vec(:,Iorb)
            deriv2(:,Iorb) = phi_vec(:,Iorb)
            v_deriv1(:,Iorb) = v(:,Iorb)
            v_deriv2(:,Iorb) = v(:,Iorb)

        
          IF(dir.eq.1) then

            call Get_d_PSI(deriv1(:,Iorb),3)
            call Get_d_PSI(deriv2(:,Iorb),2)
            call Get_d_PSI(v_deriv1(:,Iorb),3)
            call Get_d_PSI(v_deriv2(:,Iorb),2)
 

          End IF

          IF(dir.eq.2) then

            call Get_d_PSI(deriv1(:,Iorb),1)
            call Get_d_PSI(deriv2(:,Iorb),3)
            call Get_d_PSI(v_deriv1(:,Iorb),1)
            call Get_d_PSI(v_deriv2(:,Iorb),3)

          end if


          IF(dir.eq.3) then

            call Get_d_PSI(deriv1(:,Iorb),2)
            call Get_d_PSI(v_deriv1(:,Iorb),2)


            call Get_d_PSI(deriv2(:,Iorb),1)
            call Get_d_PSI(v_deriv2(:,Iorb),1)

          end if 


c========= Compute L_xyz(i,j)

          ind=1
         Do K=1,NDZ
         Do J=1,NDY
         Do I=1,NDX

c        r=DSQRT(ort_X(i)**2+ort_Y(j)**2)
            x=ort_X(i)
            y=ort_y(j)
            z=ort_z(k)
          

            IF(dir.eq.1) then
              temp(ind,Iorb)=ort_Y(J)*deriv1(ind,Iorb)-ort_Z(k)
     .*deriv2(ind,Iorb)
              v_temp(ind,Iorb)=ort_Y(J)*v_deriv1(ind,Iorb)-ort_Z(k)
     .*v_deriv2(ind,Iorb)
            end if

            IF(dir.eq.2) then
              temp(ind,Iorb)=ort_Z(k)*deriv1(ind,Iorb)-ort_x(i)
     .*deriv2(ind,Iorb)
              v_temp(ind,Iorb)=ort_Z(k)*v_deriv1(ind,Iorb)-ort_x(i)
     .*v_deriv2(ind,Iorb)
            end if

            IF(dir.eq.3) then
              temp(ind,Iorb)=ort_x(i)*deriv1(ind,Iorb)-ort_y(j)
     .*deriv2(ind,Iorb)
              v_temp(ind,Iorb)=ort_x(i)*v_deriv1(ind,Iorb)-ort_y(j)
     .*v_deriv2(ind,Iorb)
            end if


           ind=ind+1
         Enddo
         Enddo
         Enddo

          end do



       INQUIRE(FILE="Lz_orb.dat", EXIST=exists2)
       INQUIRE(FILE="Lz_orb_work.dat", EXIST=exists3)
       INQUIRE(FILE="Lz_orb_nat2.dat", EXIST=exists4)
       INQUIRE(FILE="NatVec.dat", EXIST=exists5)
 
 
       if(exists2) then
             open(unit=125, file='Lz_orb.dat', status='old',
     .position='append', action='readwrite', iostat=status_open)
           if(status_open /= 0) write(*,*)'Lz_orb.dat error'
       else
            open(unit=125, file='Lz_orb.dat',
     .status='new', action='readwrite', iostat=status_open)
              if(status_open /= 0) write(*,*)'Lz_orb.dat error!'
       end if
 
      if(exists3) then
            open(unit=126, file='Lz_orb_work.dat', status='old',
     .position='append', action='readwrite', iostat=status_open)
          if(status_open /= 0) write(*,*)'Lz_orb_work.dat error'
      else
            open(unit=126, file='Lz_orb_work.dat',
     .status='new', action='readwrite', iostat=status_open)
          if(status_open /= 0) write(*,*)'Lz_orb_work.dat error!'
      end if


      if(exists4) then
            open(unit=127, file='Lz_orb_nat2.dat', status='old',
     .position='append', action='readwrite', iostat=status_open)
          if(status_open /= 0) write(*,*)'Lz_orb_nat2.dat error'
       else
             open(unit=127, file='Lz_orb_nat2.dat',
     .status='new', action='readwrite', iostat=status_open)
           if(status_open /= 0) write(*,*)'Lz_orb_nat2.dat error!'
       end if


      if(exists5) then
             open(unit=128, file='NatVec.dat', status='old',
     .position='append', action='readwrite', iostat=status_open)
          if(status_open /= 0) write(*,*)'NatVec.dat error'
        else
              open(unit=128, file='NatVec.dat',
     .status='new', action='readwrite', iostat=status_open)
           if(status_open /= 0) write(*,*)'NatVec.dat error!'
        end if

 

       do Iorb=1, Morb
        do Jorb=1,Morb

            LPhi(Iorb,Jorb)=SUM(conjg(phi_vec(:,Iorb))*(-1.0)*(0.d0,1.d0
     .)*temp(:,Jorb)) 
            
            ang_mom=ang_mom+AllZRij(Iorb,Jorb)*LPHI(Iorb,Jorb)
      
          end do

            ang_mom_orb(Iorb)=SUM(conjg(v(:,Iorb))*(-1.0)*(0.d0,1.d0
     .)*v_temp(:,Iorb))

            ang_mom_orb_work(Iorb)=LPhi(Iorb,Iorb)

!            ang_mom_orb(Iorb)=ang_mom_orb(Iorb)/Nocc(Iorb) 

         end do

      do i=1, Morb
      do j=1, Morb
        do k=1, Morb

            LPhi(j,k)=SUM(conjg(phi_vec(:,j))*(-1.0)*(0.d0,1.d0
     .)*temp(:,k)) 


            ang_mom_orb_2(i)=ang_mom_orb_2(i)+LPHI(j,k)*
     .conjg(NatVec(j,i))*NatVec(i,k)

        end do
      end do
      end do


      write(125,526) t,(real(ang_mom_orb(jj)),jj=1,Morb)
      write(126,526) t,(real(ang_mom_orb_work(jj)),jj=1,Morb)
      write(127,526) t,(real(ang_mom_orb_work(jj)),jj=1,Morb)
      write(128,526) t,NatVec(1,1),NatVec(2,2),NatVec(3,3),NatVec(4,4)
!     .jj=1,Morb)
       close(125)
       close(126)
       close(127)
       close(128)
 



          write(*,*) "Ang-Mom: ", ang_mom


!      if(abs(t_next-t_print).le.1.0d-8.or.mod(t,0.1).eq.0.0) then
!       if(ABS(t-t_print).le.1.0d-8) then
 
       INQUIRE(FILE="Lz.dat", EXIST=exists)
 
 
      if(exists) then
            open(unit=124, file='Lz.dat', status='old',
     .position='append', action='readwrite', iostat=status_open)
          if(status_open /= 0) write(*,*)'Lz.dat error'
      else
           open(unit=124, file='Lz.dat',
     .status='new', action='readwrite', iostat=status_open)
             if(status_open /= 0) write(*,*)'Lz.dat error!'
      end if
 
       write(124,*) t,real(ang_mom)
 
       close(124)
 
!       end if


      deallocate(v)
      deallocate(temp)
      deallocate(v_temp)
      deallocate(deriv1)
      deallocate(deriv2)
      deallocate(v_deriv1)
      deallocate(v_deriv2)

526    format(12F26.16)

       END SUBROUTINE get_ang_mom

