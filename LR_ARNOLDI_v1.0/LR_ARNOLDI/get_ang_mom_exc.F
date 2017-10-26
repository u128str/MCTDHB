c#  This routine gives you the elements L_xyz(Iorb,Jorb)
c#  to compute total angular momentum L=Sum_i,j rho_ij*L_xyz(i,j)        


      SUBROUTINE  get_ang_mom_exc(u_vecE,v_vecE,dir,conv_count)
      USE   SHARED_DIMS
      USE   W_INTERPARTICLE
      USE   DVR_ALL
      USE   rR_hW
      USE USR_PAR
      IMPLICIT NONE
      integer            ::   ND,i,j,k,kk,ind,status_open,conv_count
      integer                    ::      ierr,dir,Iorb,Jorb,jj,k1,l1
      real*8                     ::      x,y,z,t,t_next,t_print
      complex*16 compt,Z1
      COMPLEX*16, DIMENSION(NDX*NDY*NDZ,conv_count,Morb) :: u_vecE
      COMPLEX*16, DIMENSION(NDX*NDY*NDZ,conv_count,Morb) :: v_vecE
      logical :: exists,exists2,exists3,exists4,exists5

      complex*16, allocatable :: tempE(:)
      complex*16, allocatable :: deriv1E(:)
      complex*16, allocatable :: deriv2E(:)
      complex*16, allocatable :: phiE(:)
      complex*16, allocatable :: phiE_new(:)
      complex*16, allocatable :: ang_momE(:,:,:)

c!=========================================================================

        ND=NDX*NDY*NDZ
       
        write(*,*) "HERE"
        if(.not.allocated(phiE)) 
     .    allocate(phiE(ND))
        write(*,*) "HERE"


        if(.not.allocated(phiE_new)) 
     .    allocate(phiE_new(ND))
        write(*,*) "HERE"

c        if(.not.allocated(u_vecE)) then
c         allocate(u_vecE(ND,conv_count,Morb))
c        end if
c        write(*,*) "HERE"

c        if(.not.allocated(v_vecE))
c     .   allocate(v_vecE(ND,conv_count,Morb))
c        write(*,*) "HERE"


        if(.not.allocated(ang_momE))
     .    allocate(ang_momE(2,conv_count,Morb))
        write(*,*) "HERE"

        if(.not.allocated(tempE))
     .    allocate(tempE(NDX*NDY*NDZ),stat=ierr)
        if(ierr /= 0) write(*,*)"allocation error for tempE"
        write(*,*) "HERE"


        if(.not.allocated(deriv1E))
     .    allocate(deriv1E(NDX*NDY*NDZ),stat=ierr)
        if(ierr /= 0) write(*,*)"allocation error for deriv1E"
        write(*,*) "HERE"

        if(.not.allocated(deriv2E))
     .   allocate(deriv2E(NDX*NDY*NDZ),stat=ierr)
        if(ierr /= 0) write(*,*)"allocation error for deriv2E"
        write(*,*) "HERE"

         
         phiE=cmplx(0.0,0.0,kind=8) 
         tempE=cmplx(0.0,0.0,kind=8) 
         ang_momE=cmplx(0.0,0.0,kind=8) 
          

 
 
            open(unit=125, file='get_ang_mom_exc.dat', status='replace',
     . action='readwrite', iostat=status_open)
           if(status_open /= 0) write(*,*)'get_ang_mom_exc.dat error'


          write(125,'(100(A15))') 
     .      "# root","   ",
     .               "(L_z(u_i),i=1..M)","   ",
     .            "(L_z(u_i),i=1..M)"
          write(125,*) '  '


!=========  Angmom for u_i
c=========  Get D/Dx,D/Dy,D/Dz
        do kk=1, conv_count 
         do Iorb=1, Morb

            phiE(:)=u_vecE(:,kk,Iorb)
     .           /sqrt(SUM(conjg(u_vecE(:,kk,Iorb))*u_vecE(:,kk,Iorb))) 

!            ind=1
!            do i=1,NDY
!             do j=1,NDX
!               phiE_new(ind)=phiE((j-1)*NDY+i)
!               ind=ind+1
!             end do
!            end do      
!            phiE=phiE_new
            
!
!            if(kk.eq.2) then
!        write(*,*) sqrt(SUM(conjg(u_vecE(:,kk,Iorb))*u_vecE(:,kk,Iorb)))
!            stop 
!            end if
            deriv1E = phiE
            deriv2E = phiE

        
          IF(dir.eq.1) then

            call Get_d_PSI(deriv1E,3)
            call Get_d_PSI(deriv2E,2)
 

          End IF

          IF(dir.eq.2) then

            call Get_d_PSI(deriv1E,1)
            call Get_d_PSI(deriv2E,3)

          end if


          IF(dir.eq.3) then

            call Get_d_PSI(deriv1E,2)
            call Get_d_PSI(deriv2E,1)

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
              tempE(ind)=ort_Y(J)*deriv1E(ind)-ort_Z(k)
     .*deriv2E(ind)
            end if

            IF(dir.eq.2) then
              tempE(ind)=ort_Z(k)*deriv1E(ind)-ort_x(i)
     .*deriv2E(ind)
            end if

            IF(dir.eq.3) then
              tempE(ind)=ort_x(i)*deriv1E(ind)-ort_y(j)
     .*deriv2E(ind)
!              tempE(ind)=ort_x(i)*deriv2E(ind)-ort_y(j)
!     .*deriv1E(ind)
            end if


           ind=ind+1
         Enddo
         Enddo
         Enddo


            ang_momE(1,kk,Iorb)=
     .           SUM(conjg(phiE(:))
     .         *(-1.0)*(0.d0,1.d0)*tempE(:))/SUM(conjg(phiE(:))*phiE(:))


         end do
        end do




!========  Angmom for v_i
        do kk=1, conv_count 
         do Iorb=1, Morb

            phiE(:)=v_vecE(:,kk,Iorb)
            deriv1E = phiE
            deriv2E = phiE

        
          IF(dir.eq.1) then

            call Get_d_PSI(deriv1E,3)
            call Get_d_PSI(deriv2E,2)
 

          End IF

          IF(dir.eq.2) then

            call Get_d_PSI(deriv1E,1)
            call Get_d_PSI(deriv2E,3)

          end if


          IF(dir.eq.3) then

            call Get_d_PSI(deriv1E,2)
            call Get_d_PSI(deriv2E,1)

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
              tempE(ind)=ort_Y(J)*deriv1E(ind)-ort_Z(k)
     .*deriv2E(ind)
            end if

            IF(dir.eq.2) then
              tempE(ind)=ort_Z(k)*deriv1E(ind)-ort_x(i)
     .*deriv2E(ind)
            end if

            IF(dir.eq.3) then
              tempE(ind)=ort_x(i)*deriv1E(ind)-ort_y(j)
     .*deriv2E(ind)
            end if


           ind=ind+1
         Enddo
         Enddo
         Enddo


            ang_momE(2,kk,Iorb)=SUM(conjg(phiE(:))*(-1.0)*(0.d0,1.d0)
     .                   *tempE(:)) 


         end do
        end do


 
!====== OUTPUT
       do kk=1, conv_count

        ind=1
!       do kk=1, NDY
!        do j=1, NDX
         

         write(125,'(I5,10x,100(F20.10,10x))') kk,
     .                 ((Real(ang_momE(j,kk,i)),i=1,Morb,1),j=1,2,1)
!         write(125,'(3(F20.10,10x))') ort_X(j),ort_Y(kk),
!     .                 Dimag(u_vecE(ind,2,1))
         ind=ind+1
 
!        end do
         write(125,*) '  '
       end do 

       close(125)


      deallocate(tempE)
      deallocate(deriv1E)
      deallocate(deriv2E)


       END SUBROUTINE get_ang_mom_exc

