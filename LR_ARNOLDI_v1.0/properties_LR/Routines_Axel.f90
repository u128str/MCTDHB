MODULE ANALYZER
use Prop_MB
USE W_INTERPARTICLE
USE CORRELATIONFUNCTIONS

IMPLICIT NONE
CONTAINS


!###############################################################
!         This routine integrates the one body density 
!    from xi to xf and writes the output to the file nonescape
!###############################################################
subroutine density_nonescape(xi,xf,rho_jk,time,PSI)

use CI_ALL,ONLY:NPar
use DVR_ALL
USE W_INTERPARTICLE

  IMPLICIT NONE

  COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb)  :: PSI
  complex*16,DIMENSION(Morb,Morb),intent(in) :: rho_jk
  real*8, intent(in) :: xi,xf
  real*8 :: nonescape,time,w
  
  complex*16, DIMENSION(NDX*NDY*NDZ)  :: DENSITY    

  integer :: i,j,k,l,m,n,NOGRIDPOINTS


  w=DSQRT((xf-xi)/NDX)

  if(DIM_MCTDHB.gt.1) then
    write(*,*) 'For Dimension>1 you have to code that yourself ;)'
  endif

  NOGRIDPOINTS=NDX*NDY*NDZ
! get the density from the one-body-matrix elements 
  DENSITY   = (0.d0,0.d0)
  do n=1,NOGRIDPOINTS
     do j = 1,Morb  
        do k=1,MOrb
           DENSITY(n) = DENSITY(n) + rho_jk(j,k) * DCONJG(PSI(n,j))*PSI(n,k)
        end do
     end do
  end do

! integrate from DENSITY from xi to xf
  nonescape=0.d0
  do n=1,NOGRIDPOINTS
     if ((Ort_X(n).lt.xf).and.(Ort_X(n).gt.xi)) then
        nonescape=nonescape +dble(DENSITY(n))
     endif
  end do
! open/write/close output
  open(unit=542,file='nonescape',form='formatted',ACCESS='APPEND')

  write(542,6545) time,' ',nonescape 

6545  format(2(F21.16,a1))

  close(542)
end subroutine density_nonescape

!computes the expectation values for N=2 particle loss operators when partitioning 
!the two particle Hilbert space at 'border'

subroutine lossops_two_bosons(time,PSI,rho_ijkl,border)
use CI_ALL,ONLY:NPar
use DVR_ALL
USE W_INTERPARTICLE

IMPLICIT NONE

COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb)  :: PSI
complex*16,DIMENSION(Morb,Morb,Morb,Morb),intent(in) :: rho_ijkl
complex*16 :: G2m_n
real*8, intent(in) :: time,border
real*8 :: int_2in,int_2out,int_1in_1out,int_1out_1in,check

integer            :: NOGRIDPOINTS,LORB,KORB,JORB,IORB,m,n,i_m,i_n,j_m,j_n,k_m,k_n
character(len=10)   :: borderAsString
character(len=100)   :: fullName


int_2in=0.d0
int_2out=0.d0
int_1in_1out=0.d0
int_1out_1in=0.d0

if (NPar.ne.2) then 
   write(*,*) 'this is so far only implemented for N=2'
   return
endif

if (DIM_MCTDHB.ne.1) then
   write(*,*) 'this is so far only implemented for 1D systems'
   return
endif
NOGRIDPOINTS = NDX*NDY*NDZ

  do n=1,NOGRIDPOINTS
     do m=1,NOGRIDPOINTS
        call get_ijk_from_m(m,NDX,NDY,i_m,j_m,k_m)
        call get_ijk_from_m(n,NDX,NDY,i_n,j_n,k_n)
        if (((ort_x(i_m).gt.border).and.(ort_x(i_n).gt.border)).eqv..FALSE.) then
        G2m_n=(0.d0,0.d0)
        do Lorb=1,Morb
           do Korb=1,Morb
              do Jorb=1,Morb
                 do Iorb=1,Morb
                    G2m_n = G2m_n + rho_ijkl(Iorb,Jorb,Korb,Lorb)* &
                              DCONJG(PSI(m,Iorb))*DCONJG(PSI(n,Jorb))*PSI(m,Korb)*PSI(n,Lorb)
  
                 end do
              end do
           end do
        end do


        if ((ort_x(i_m).lt.border).and.(ort_x(i_n).lt.border)) int_2in=int_2in+Dble(G2m_n)
        if ((ort_x(i_m).ge.border).and.(ort_x(i_n).lt.border)) int_1in_1out=int_1in_1out+2*Dble(G2m_n)
        
!        if ((ort_x(i_m).ge.border).and.(ort_x(i_n).ge.border)) int_2out=int_2out+Dble(G2m_n)
!        if ((ort_x(i_n).ge.border).and.(ort_x(i_m).lt.border)) int_1out_1in=int_1out_1in+Dble(G2m_n)
        endif
        
     end do
  end do
  int_2out=2-int_2in-int_1in_1out
  call getDoubleAsString(border,borderAsString) 
  fullName='lossops_N2_'//trim(borderAsString)//'.dat'
  
  open(unit=12,file=trim(fullName),form='formatted')

  write(12,8765) time,int_2in,int_1in_1out,int_2out
 
 8765 FORMAT(5E25.16)

return
  
end subroutine lossops_two_bosons

subroutine get_corr_one_restricted(time,PSI,mode,dilation,rho_jk,rho_ijkl,&
                                   xini,xfin,xpts,kxini,kxfin,kpts) 
!--------------------------------------------------------------------------------------------
! computes the correlation functions 
! G1(r,r), G1(r,r'), G1(r',r') where r=(x), 
! and writes them to a file. The first 2 columns contain
! x,x'. 
!
! time          : current time 
! PSI           : orbital array 
! mode          : mode=1: x-space (input PSI=PSI) 
!                 mode=2: k-space (input PSI=FT[PSI])
! rho_jk        : the 1-body density matrix where rho_jk(j,k)   = rho_jk 
! rho_ijkl      : the 2-body density matrix where rho_ijkl(j,k) = rho_ijkl 
! xini,xfin     : spatial region, for which data is written for xpts
!--------------------------------------------------------------------------------------------
!  {{{
use CI_ALL,ONLY:Npar
use DVR_ALL
  implicit none
  integer  :: NOGRIDPOINTS,last_i_m,last_i_n
  integer,intent(in) :: mode,dilation,xpts,kpts
  real*8,intent(in)  :: time
  real*8     :: w,dk,mom_X(NDX),mom_Y(NDY),mom_Z(NDZ)
  COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb)  :: PSI
  complex*16,intent(in) :: rho_ijkl(Morb,Morb,Morb,Morb)
  complex*16,intent(in) :: rho_jk(Morb,Morb)
  character(len=2)   :: aString
  character(len=200) :: fullName
  character(len=10)  :: timeAsString
  character(len=10)  :: MorbAsString
  character(len=1)   :: NparAsString1
  character(len=2)   :: NparAsString2
  character(len=3)   :: NparAsString3
  character(len=4)   :: NparAsString4
  character(len=5)   :: NparAsString5
  character(len=6)   :: NparAsString6
  character(len=7)   :: NparAsString7
  character(len=8)   :: NparAsString8

  integer            :: m,n,Iorb,Jorb,Korb,Lorb,i_m,j_m,k_m,i_n,j_n,k_n
  complex*16         :: G1m_m,G1m_n,G1n_n,G2m_n

  complex*16         :: maximum,integ
  real*8,parameter   :: TOL=0.00001d0 
  real*8     :: threshold
  real*8     :: outG1m_m,outReG1m_n,outImG1m_n,outG1n_n,outG2m_n,increment
  real*8             :: xi,xf,yi,yf,zi,zf,xini,xfin,kxini,kxfin,onebdens(NDX)

  integ=0.d0
  NOGRIDPOINTS=NDX*NDY*NDZ
   xi=Time_xint
   xf=Time_xfnl
   yi=Time_yint
   yf=Time_yfnl
   zi=Time_zint
   zf=Time_zfnl
  if (mode.eq.1) then
     increment= (xfin-xini)/(xpts*1.d0+2.d0)
  else if (mode.eq.2) then 
     increment=(kxfin-kxini)/(kpts*1.d0+2.d0)
  endif

! some exceptions
  if (Time_DVRMETHODX.eq.1) then
     write(*,*) 'The correlations in restricted space work only aequidistant grids...'
     return
  else if ((mode.eq.1).and.(increment.le.ort_x(2)-ort_x(1))) then
     write(*,*) 'CORR1-RESTRICTED:The number of points you asked implies a denser grid than the output has,&
          I will use the full number - I decrease the number of points accordingly. Mode',mode
       increment=0.d0
  else if ((mode.eq.2).and.(increment.le.abs(mom_x(2)-mom_x(1)))) then 
       write(*,*) 'CORR1-RESTRICTED:The number of points you asked implies a denser grid than the output has,&
           I will use the full number - I decrease the number of points accordingly. Mode',mode
       increment=0.d0
  elseif ((xi.gt.xini).or.(xfin.gt.xf)) then
       write(*,*) 'CORR1-RESTRICTED:The asked grid is bigger than the one in the output - i decrease the size accordingly!'
       if ((xi.gt.xini)) xini = xi
       if ((xfin.gt.xf)) xfin = xf
  endif
  if (mode.eq.1) then 
    IF (TIME_DVRMETHODX.ne.3) then
      w= DSQRT((xf-xi)/(NDX*1.d0))
    else if (TIME_DVRMETHODX.eq.3) then
      w= DSQRT((xf-xi)/(NDX*1.d0+1.d0))
    endif   
  else if (mode.eq.2) then
     call get_mom(xi,xf,NDX,mom_X)
     dk     = 2.d0*PI/(xf-xi)
     w      = DSQRT( dk )
  else if ((mode.ne.1).and.(mode.ne.2)) then
        write(*,*) "wrong mode"
  end if

!  write (*,*) 'the grid spacing for the correlations is:', increment 

  if (MORB.lt.10) then
     write (MorbAsString, '(I1)') Morb
  else if ((MORB.ge.10).and.(MORB.lt.100)) then
     write (MorbAsString, '(I2)') Morb
  else 
     write(*,*) 'Bigger orbital number than 100 not implemented!'
  endif
  
  call getDoubleAsString(time,timeAsString) 

  if (mode.eq.1) then
     aString = 'x-'
  else if (mode.eq.2) then 
     aString = 'k-'
  else 
     write(*,*)"ERROR in get_correlations"
     stop
  endif   

  if ((NPar.gt.0) .and. (NPar.lt.10)) then 
    write (NParAsString1, '(I1)') NPar
    fullName=timeasString//'N'//trim(NParAsString1)//'M'//trim(MorbAsString)//aString//'corr1restr.dat'
  else  if ((NPar.ge.10) .and. (NPar.lt.100)) then 
    write (NParAsString2, '(I2)') NPar
    fullName=timeasString//'N'//trim(NParAsString2)//'M'//trim(MorbAsString)//aString//'corr1restr.dat'
  else if ((NPar.ge.100) .and. (NPar.lt.1000)) then 
    write (NParAsString3, '(I3)') NPar
    fullName=timeasString//'N'//trim(NParAsString3)//'M'//trim(MorbAsString)//aString//'corr1restr.dat'
  else if ((NPar.ge.1000) .and. (NPar.lt.10000)) then 
    write (NParAsString4, '(I4)') NPar
    fullName=timeasString//'N'//trim(NParAsString4)//'M'//trim(MorbAsString)//aString//'corr1restr.dat'
  else if ((NPar.ge.10000) .and. (NPar.lt.100000)) then 
    write (NParAsString5, '(I5)') NPar
    fullName=timeasString//'N'//trim(NParAsString5)//'M'//trim(MorbAsString)//aString//'corr1restr.dat'
  else if ((NPar.ge.100000) .and. (NPar.lt.1000000)) then 
    write (NParAsString6, '(I6)') NPar
    fullName=timeasString//'N'//trim(NParAsString6)//'M'//trim(MorbAsString)//aString//'corr1restr.dat'
  else if ((NPar.ge.1000000) .and. (NPar.lt.10000000)) then 
    write (NParAsString7, '(I7)') NPar
    fullName=timeasString//'N'//trim(NParAsString7)//'M'//trim(MorbAsString)//aString//'corr1restr.dat'
  else if ((NPar.ge.10000000) .and. (NPar.lt.100000000)) then 
    write (NParAsString8, '(I8)') NPar
    fullName=timeasString//'N'//trim(NParAsString8)//'M'//trim(MorbAsString)//aString//'corr1restr.dat'

  else 
    NParAsString4='XXXX'
    fullName=timeasString//'N'//trim(NParAsString4)//'M'//trim(MorbAsString)//aString//'corr1restr.dat'
  end if
  
  open(unit=12,file=trim(fullName),form='formatted')

!compute onnebody density
        onebdens = (0.d0,0.d0) 
  do m=1,NDX
        G1m_m = (0.d0,0.d0) 
        do Jorb = 1,Morb  
           do Korb=1,Morb
              G1m_m = G1m_m + rho_jk(Jorb,Korb) * DCONJG(PSI(m,Jorb))*PSI(m,Korb)
           end do
        end do
        onebdens(m)=dble(G1m_m)/w**2
  end do
  integ=0.d0
   
!compute G1n_n
  last_i_m=NDX/2
  last_i_n=NDX/2
  i_n=1
  i_m=1
!---------------
  do n=1,NOGRIDPOINTS
    call get_ijk_from_m(n,NDX,NDY,i_n,j_n,k_n)
    if ((((mode.eq.1).and.&
         abs(ort_x(last_i_n)-ort_x(i_n)).ge.increment).and.&
         (Ort_X(i_n).lt.xfin).and.(Ort_x(i_n).gt.xini)).or. &
         ((mode.eq.2).and.&
         (abs(mom_x(last_i_n)-mom_x(i_n)).ge.increment).and.&
         (mom_X(i_n).lt.kxfin).and.(mom_X(i_n).gt.kxini)))&
         then
     do m=1,NOGRIDPOINTS
        call get_ijk_from_m(m,NDX,NDY,i_m,j_m,k_m)
        if (((mode.eq.1).and.&
           (Ort_X(i_m).lt.xfin).and.&
           (Ort_x(i_m).gt.xini).and.&
           (abs(Ort_X(i_m)-Ort_X(last_i_m)).ge.increment)).or.&
           ((mode.eq.2).and.&
           (mom_X(i_m).lt.kxfin).and.&
           (mom_X(i_m).gt.kxini).and.&
           (abs(mom_X(i_m)-mom_X(last_i_m)).ge.increment)))&
           then
!compute G1m_m
        
!compute G1m_n
        G1m_n = (0.d0,0.d0) 
        do Jorb = 1,Morb  
           do Korb=1,Morb
              G1m_n = G1m_n + rho_jk(Jorb,Korb) * DCONJG(PSI(m,Jorb))*PSI(n,Korb)
           end do
        end do
        
        G1m_n=G1m_n*1.d0/w**2 
        outReG1m_n=DBLE(G1m_n)   
        outImG1m_n=DIMAG(G1m_n)
        if (m.eq.n) integ=integ+g1m_n
!----------------
!#########################################################################################
!####### output only in restricted space 
!#########################################################################################
           if (mode.eq.1) then
                write(12,6669) ort_X(i_m),&
                      ort_X(i_n),outREG1m_n,outIMG1m_n,&
                      onebdens(m),onebdens(n)
           else if (mode.eq.2) then
                write(12,6669) mom_X(i_m),&
                      mom_X(i_n),outREG1m_n,outIMG1m_n,&
                      onebdens(m),onebdens(n)
           endif
      else !if conditions for first coordinate are not met
            goto 9877
      endif
      last_i_m=i_m
 9877 continue
      if(MOD(i_m,NDX).eq.0) write(12,*)"                                          " 
     end do
     else !if conditions for second coordinate are not met 
       goto 9878
     endif
 9878 continue
  end do 
  close(12)
 
 
6669  format(2(F8.3),4(E25.16))
!---}}}
end subroutine get_corr_one_restricted

subroutine get_corr_two_restricted(time,PSI,mode,dilation,rho_jk,rho_ijkl,&
              xini,xfin,xpts,kxini,kxfin,kpts)
!--------------------------------------------------------------------------------------------
! computes the correlation functions 
! G1(r,r), G1(r,r'), G1(r',r') and G2(r,r',r',r), where r=(x), 
! and writes them to a file. The first 2 columns contain
! x,x'. Currently only 1D has been tested.
!
! time          : current time 
! PSI           : orbital array 
! mode          : mode=1: x-space (input PSI=PSI) 
!                 mode=2: k-space (input PSI=FT[PSI])
! rho_jk        : the 1-body density matrix where rho_jk(j,k)   = rho_jk 
! rho_ijkl      : the 2-body density matrix where rho_ijkl(j,k) = rho_ijkl 
! xini,xfin         : spatial region, for which data is written
! kxini,kxfin         : spatial region, for which data is written
! xpts,kpts                : number of points for x/k space output
!--------------------------------------------------------------------------------------------
!  {{{
use CI_ALL,ONLY:Npar
use DVR_ALL
  implicit none
  integer  :: NOGRIDPOINTS,last_i_m,last_i_n,L,I,J
  integer,intent(in) :: mode,dilation,xpts,kpts
  real*8,intent(in)  :: time
  real*8     :: w,dk,mom_X(NDX),mom_Y(NDY),mom_Z(NDZ)
   COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb)  :: PSI
  complex*16,intent(in) :: rho_ijkl(Morb,Morb,Morb,Morb)
  complex*16,intent(in) :: rho_jk(Morb,Morb)
  character(len=2)   :: aString
  character(len=200) :: fullName
  character(len=10)  :: timeAsString
  character(len=10)  :: MorbAsString
  character(len=1)   :: NparAsString1
  character(len=2)   :: NparAsString2
  character(len=3)   :: NparAsString3
  character(len=4)   :: NparAsString4
  character(len=5)   :: NparAsString5
  character(len=6)   :: NparAsString6
  character(len=7)   :: NparAsString7
  character(len=8)   :: NparAsString8

  integer            :: m,n,Iorb,Jorb,Korb,Lorb,i_m,j_m,k_m,i_n,j_n,k_n
  complex*16         :: G1m_m,G1m_n,G1n_n,G2m_n,corr2,denscheck,integ

  real*8         :: onebdens(NDX)
  real*8,parameter   :: TOL=0.00001d0 
  real*8     :: threshold
  real*8     :: outG1m_m,outReG1m_n,outImG1m_n,outG1n_n,outG2m_n
  real*8             :: xi,xf,yi,yf,zi,zf,xini,xfin,kxini,kxfin
  real*8     ::   increment


    IF (TIME_DVRMETHODX.eq.1) then
       write(*,*) 'Please code the analysis 2-body correlations on restricted subspaces for HO DVR!!!'
       return
    endif

   

  NOGRIDPOINTS=NDX*NDY*NDZ
   xi=Time_xint
   xf=Time_xfnl
   yi=Time_yint
   yf=Time_yfnl
   zi=Time_zint
   zf=Time_zfnl
   integ=0.d0
  if (mode.eq.1) then 
    IF (TIME_DVRMETHODX.ne.3) then
      w= DSQRT((xf-xi)/(NDX*1.d0))
    else if (TIME_DVRMETHODX.eq.3) then
      w= DSQRT((xf-xi)/(NDX*1.d0+1.d0))
    endif   
  else if (mode.eq.2) then
     call get_mom(xi,xf,NDX,mom_X)
     dk     = 2.d0*PI/(xf-xi)
     w      = DSQRT( dk )

  else if ((mode.ne.1).and.(mode.ne.2)) then
        write(*,*) "wrong mode"
  end if

  if (mode.eq.1) then
     increment= (xfin-xini)/(xpts*1.d0+2.d0)
  else if (mode.eq.2) then 
     increment=(kxfin-kxini)/(kpts*1.d0+2.d0)
  endif

! some exceptions
  if (Time_DVRMETHODX.eq.1) then
     write(*,*) 'The correlations in restricted space work only aequidistant grids...'
     return
  else if ((mode.eq.1).and.(increment.le.abs(ort_x(2)-ort_x(1)))) then
       write(*,*) 'CORR2-RESTRICTED:The number of points you asked implies a denser grid than the output has,&
           I will use the full number - I decrease the number of points accordingly. Mode',mode
       increment=0.d0
  else if ((mode.eq.2).and.(increment.le.abs(mom_x(2)-mom_x(1)))) then 
       write(*,*) 'CORR2-RESTRICTED:The number of points you asked implies a denser grid than the output has,&
           I will use the full number - I decrease the number of points accordingly. Mode',mode
       increment=0.d0
  elseif ((xi.gt.xini).or.(xfin.gt.xf)) then
       write(*,*) 'The asked grid is bigger than the one in the output - i decrease the size accordingly!'
       if ((xi.gt.xini)) xini = xi
       if ((xfin.gt.xf)) xfin = xf
  endif

!  if (mode.eq.2) then
!     write(*,*) 'Increment:',increment,'grid-spacing:',mom_X(2)-mom_X(1)
!  else if (mode.eq.1) then
!     write(*,*) 'Increment:',increment,'grid-spacing:',ort_X(2)-ort_X(1)
!  endif

! write (*,*) 'the grid spacing for the correlations is:', increment 
  if (MORB.lt.10) then
     write (MorbAsString, '(I1)') Morb
  else if ((MORB.ge.10).and.(MORB.lt.100)) then
     write (MorbAsString, '(I2)') Morb
  else 
     write(*,*) 'Bigger orbital number than 100 not implemented!'
  endif


  call getDoubleAsString(time,timeAsString) 

  if (mode.eq.1) then
     aString = 'x-'
  else if (mode.eq.2) then 
     aString = 'k-'
  else 
     write(*,*)"ERROR in get_correlations"
     stop
  endif   

!!!!!!!!!!!!!!!!!!! IF YOU WANT THE SUM OF two-body MATRIX ELEMENTS TO BE DISPLAYED
!  write(*,*) 'ANALYZER:', sum(abs(rho_ijkl)),   &
!  sum(rho_ijkl)
!!!!!!!!!!!!!!!!!!! IF YOU WANT THE SUM OF two-body MATRIX ELEMENTS TO BE DISPLAYED
  
  if ((NPar.gt.0) .and. (NPar.lt.10)) then 
    write (NParAsString1, '(I1)') NPar
    fullName=timeasString//'N'//trim(NParAsString1)//'M'//trim(MorbAsString)//aString//'corr2restr.dat'
  else  if ((NPar.ge.10) .and. (NPar.lt.100)) then 
    write (NParAsString2, '(I2)') NPar
    fullName=timeasString//'N'//trim(NParAsString2)//'M'//trim(MorbAsString)//aString//'corr2restr.dat'
  else if ((NPar.ge.100) .and. (NPar.lt.1000)) then 
    write (NParAsString3, '(I3)') NPar
    fullName=timeasString//'N'//trim(NParAsString3)//'M'//trim(MorbAsString)//aString//'corr2restr.dat'
  else if ((NPar.ge.1000) .and. (NPar.lt.10000)) then 
    write (NParAsString4, '(I4)') NPar
    fullName=timeasString//'N'//trim(NParAsString4)//'M'//trim(MorbAsString)//aString//'corr2restr.dat'
  else if ((NPar.ge.10000) .and. (NPar.lt.100000)) then 
    write (NParAsString5, '(I5)') NPar
    fullName=timeasString//'N'//trim(NParAsString5)//'M'//trim(MorbAsString)//aString//'corr2restr.dat'
  else if ((NPar.ge.100000) .and. (NPar.lt.1000000)) then 
    write (NParAsString6, '(I6)') NPar
    fullName=timeasString//'N'//trim(NParAsString6)//'M'//trim(MorbAsString)//aString//'corr2restr.dat'
  else if ((NPar.ge.1000000) .and. (NPar.lt.10000000)) then 
    write (NParAsString7, '(I7)') NPar
    fullName=timeasString//'N'//trim(NParAsString7)//'M'//trim(MorbAsString)//aString//'corr2restr.dat'
  else if ((NPar.ge.10000000) .and. (NPar.lt.100000000)) then 
    write (NParAsString8, '(I8)') NPar
    fullName=timeasString//'N'//trim(NParAsString8)//'M'//trim(MorbAsString)//aString//'corr2restr.dat'
 
  else 
    NParAsString4='XXXX'
    fullName=timeasString//'N'//trim(NParAsString4)//'M'//trim(MorbAsString)//aString//'corr2restr.dat'
  end if

  open(unit=12,file=trim(fullName),form='formatted')

!compute onnebody density
        onebdens = (0.d0,0.d0) 
  do m=1,NDX
        G1m_m = (0.d0,0.d0) 
        do Jorb = 1,Morb  
           do Korb=1,Morb
              G1m_m = G1m_m + rho_jk(Jorb,Korb) * DCONJG(PSI(m,Jorb))*PSI(m,Korb)
           end do
         end do
        onebdens(m)=dble(G1m_m*1.d0/(w**2))
  end do
!  write(*,*) 'N=', SUM(onebdens*w**2),'in corr2 restr'
!  write(*,*) 'computed weight', w,'in corr2 restr'

  last_i_m=NDX/2
  last_i_n=NDX/2
  i_n=1
  i_m=1
  integ=0.d0
  do n=1,NOGRIDPOINTS
    call get_ijk_from_m(n,NDX,NDY,i_n,j_n,k_n)
    if ((((mode.eq.1).and.&
         abs(ort_x(last_i_n)-ort_x(i_n)).ge.increment).and.&
         (Ort_X(i_n).lt.xfin).and.(Ort_x(i_n).gt.xini)).or. &
         ((mode.eq.2).and.&
         (abs(mom_x(last_i_n)-mom_x(i_n)).ge.increment).and.&
         (mom_X(i_n).lt.kxfin).and.(mom_X(i_n).gt.kxini)))&
         then
     do m=1,NOGRIDPOINTS
        call get_ijk_from_m(m,NDX,NDY,i_m,j_m,k_m)
        if (((mode.eq.1).and.&
           (Ort_X(i_m).lt.xfin).and.&
           (Ort_x(i_m).gt.xini).and.&
           (abs(Ort_X(i_m)-Ort_X(last_i_m)).ge.increment)).or.&
           ((mode.eq.2).and.&
           (mom_X(i_m).lt.kxfin).and.&
           (mom_X(i_m).gt.kxini).and.&
           (abs(mom_X(i_m)-mom_X(last_i_m)).ge.increment)))&
         then
!compute G2m_n 
        G2m_n=(0.d0,0.d0)
        do Lorb=1,Morb
           do Korb=1,Morb
              do Jorb=1,Morb
                 do Iorb=1,Morb
!                       RhoAll(Iorb,Jorb,Korb,Lorb)
                    G2m_n = G2m_n + rho_ijkl(Iorb,Jorb,Korb,Lorb)* &
                              DCONJG(PSI(m,Iorb))*DCONJG(PSI(n,Jorb))*PSI(m,Korb)*PSI(n,Lorb)
                 end do
              end do
           end do
        end do
        G2m_n=G2m_n*(1.d0/w**4)
           if (mode.eq.1) then
                write(12,6666) ort_X(i_m),&
                      ort_X(i_n),G2m_n,&
                      onebdens(m),onebdens(n)
           else if (mode.eq.2) then
                write(12,6666) mom_X(i_m),&
                      mom_X(i_n),G2m_n,&
                      onebdens(m),onebdens(n)
           endif
        else !if conditions for first coordinate are not met
            goto 9879
        endif
!        if (mode.eq.2) then
!             write(*,*) 'two body restricted momentum correlation function not implemented so far'
!        end if
        last_i_m=i_m
 9879 continue
        if(MOD(i_m,NDX).eq.0) write(12,*)"                                 " 
     end do
     else !if conditions for second coordinate are not met 
       goto 9880
     endif
     last_i_n=i_n                        
 9880 continue
  end do  

  close(12)
 
 
 6666  format(2(F8.3),4(E25.16))
!---}}}
end subroutine get_corr_two_restricted


subroutine get_corr_slice(time,PSI,mode,rho_jk,rho_ijkl,&
                       x1const,x1slice,y1const,y1slice,x2const,x2slice,y2const,y2slice,dilation) 
!--------------------------------------------------------------------------------------------
! computes cuts of the correlation functions 
! G1(r,r), G1(r,r'), G1(r',r') and G2(r,r',r',r), where r=(x,y), 
! and writes them to a file.  
! 
! the user has to choose two of the coordinates to stay constant, to obtain a 3D plot.
! x1,y1,x2,y2. Currently only 1D has been tested.
!
!
! time          : current time 
! PSI           : orbital array 
! mode          : mode=1: x-space (input PSI=PSI) 
!                 mode=2: k-space (input PSI=FT[PSI])
! dilation      : the integer dilation gives the factor by which the momentum intervall should
!                 be shrunk or equivalently the factor by which the sampling rate is increased.
!                 The momentum distribution at the outer ends is padded with zeros then. 
!                 Normally dilation is computed by get_dilation. If unsure put
!                 dilation = 1, but the momentum grid will then usually have a bad
!                 resolution.
! rho_jk        : the 1-body density matrix where rho_jk(j,k)   = rho_jk 
! rho_ijkl      : the 2-body density matrix where rho_ijkl(j,k) = rho_ijkl 
!--------------------------------------------------------------------------------------------
!  {{{
use CI_ALL,ONLY:Npar
use DVR_ALL
  implicit none
!  integer,parameter  :: NOGRIDPOINTS=NDX*NDY*NDZ
  integer  :: NOGRIDPOINTS
  integer,intent(in) :: mode
  integer,intent(in) :: dilation
  real*8,intent(in)  :: time
  real*8,intent(in)  :: x1slice,x2slice,y1slice,y2slice
  logical,intent(in)  :: x1const,y1const,x2const,y2const
  logical  :: x1done,y1done,x2done,y2done

  real*8     :: w,dk,mom_X(NDX),mom_Y(NDY),mom_Z(NDZ)
!  complex*16,DIMENSION(NOGRIDPOINTS,Morb),intent(in) :: PSI 
   COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb)  :: PSI
  complex*16,intent(in) :: rho_ijkl(Morb,Morb,Morb,Morb)
  complex*16,intent(in) :: rho_jk(Morb,Morb)
  character(len=2)   :: aString
  character(len=200) :: fullName   
  character(len=10)  :: timeAsString
  character(len=7)  :: slice1AsString
  character(len=7)  :: slice2AsString
  real*8            :: slice1,slice2
  integer           :: slice1pt,slice2pt
  character(len=10)  :: MorbAsString
  character(len=1)   :: NparAsString1
  character(len=2)   :: NparAsString2
  character(len=3)   :: NparAsString3
  character(len=4)   :: NparAsString4
  character(len=5)   :: NparAsString5
  character(len=6)   :: NparAsString6
  character(len=7)   :: NparAsString7
  character(len=8)   :: NparAsString8
  integer            :: m,n,Iorb,Jorb,Korb,Lorb,i_m,j_m,k_m,i_n,j_n,k_n,k_k
  complex*16         :: G1m_m,G1m_n,G1n_n,G2m_n
  Real*8             :: onebdens(NDX*NDY)

  complex*16         :: maximum,integ,integ2
  real*8,parameter   :: TOL=0.00001d0 
  real*8     :: threshold
  real*8     :: outG1m_m,outReG1m_n,outImG1m_n,outG1n_n,outG2m_n
  real*8             :: xi,xf,yi,yf,zi,zf

  if (NDY.eq.1) then
     write(*,*) 'The sliced correlations are implemented for >2D systems only!!'
     return
  endif

 
  integ=0.d0
  integ2=0.d0


  NOGRIDPOINTS=NDX*NDY*NDZ
  xi=Time_xint
  xf=Time_xfnl
  yi=Time_xint
  yf=Time_xfnl
  if(mode.eq.1) then
      w=weight(1)
  end if
  if (mode.eq.2) then
     call get_mom(xi,xf,NDX,mom_X)
     call get_mom(yi,yf,NDY,mom_Y)
     mom_X=mom_X/dilation
     mom_Y=mom_Y/dilation
     dk  = 2.d0*PI/((xf-xi)*dilation)*2.d0*PI/((yf-yi)*dilation)
     w   = sqrt( dk )


  end if
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! some exceptions 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (mode.eq.1) then

  if ((x1const.eqv..true.).and.((x1slice.lt.xi).or.(x1slice.gt.xf))) then
     write(*,*) 'The wanted x1 slice of the correlations is outside the grid!!!'
     return 
  endif
  if ((x2const.eqv..true.).and.((x2slice.lt.xi).or.(x2slice.gt.xf))) then
     write(*,*) 'The wanted x1 slice of the correlations is outside the grid!!!'
     return 
  endif
  if ((y1const.eqv..true.).and.((y1slice.lt.yi).or.(y1slice.gt.yf))) then
     write(*,*) 'The wanted x1 slice of the correlations is outside the grid!!!'
     return 
  endif
  if ((y2const.eqv..true.).and.((y2slice.lt.xi).or.(y2slice.gt.xf))) then
     write(*,*) 'The wanted x1 slice of the correlations is outside the grid!!!'
     return 
  endif

  else if (mode.eq.2) then

  if ((x1const.eqv..true.).and.((x1slice.gt.mom_x(NDX)).or.(x1slice.lt.mom_x(1)))) then
     write(*,*) 'The wanted x1 slice of the correlations is outside the grid!!!'
     return 
  endif
  if ((x2const.eqv..true.).and.((x2slice.gt.mom_x(ndx)).or.(x2slice.lt.mom_x(1)))) then
     write(*,*) 'The wanted x1 slice of the correlations is outside the grid!!!'
     return 
  endif
  if ((y1const.eqv..true.).and.((y1slice.gt.mom_y(ndx)).or.(y1slice.lt.mom_y(1)))) then
     write(*,*) 'The wanted x1 slice of the correlations is outside the grid!!!'
     return 
  endif
  if ((y2const.eqv..true.).and.((y2slice.gt.mom_y(ndx)).or.(y2slice.lt.mom_y(1)))) then
     write(*,*) 'The wanted x1 slice of the correlations is outside the grid!!!'
     return 
  endif

  endif 
  m=0


  if (x1const.eqv..true.) then
      m=m+1
      call getsliceasstring(1,x1slice,slice1AsString)
      slice1=x1slice
  endif
  if (y1const.eqv..true.) then
      m=m+1
      if (m.eq.1) then
         call getsliceasstring(2,y1slice,slice1AsString)
         slice1=y1slice
      endif
      if (m.eq.2) then
         call getsliceasstring(2,y1slice,slice2AsString)
         slice2=y1slice
      endif
  endif
  if (x2const.eqv..true.) then 
      m=m+1
      if (m.eq.1) then 
         call getsliceasstring(3,x2slice,slice1AsString)
         slice1=x2slice
      endif
      if (m.eq.2) then 
        call getsliceasstring(3,x2slice,slice2AsString)
        slice2=x2slice
      endif
  endif
  if (y2const.eqv..true.) then 
     m=m+1
     if (m.eq.1) then
        call getsliceasstring(4,y2slice,slice1AsString)
        slice1=y2slice
      endif
     if (m.eq.2) then
        call getsliceasstring(4,y2slice,slice2AsString)
        slice2=y2slice
      endif
  endif

  if (m.ne.2) then
      write(*,*) 'You must keep two of the four coordinates of the correlations constant!!!'
      return
  endif

 
  if (MORB.lt.10) then
     write (MorbAsString, '(I1)') Morb
  else if ((MORB.ge.10).and.(MORB.lt.100)) then
     write (MorbAsString, '(I2)') Morb
  else 
     write(*,*) 'Bigger orbital number than 100 not implemented!'
  endif

  call getDoubleAsString(time,timeAsString) 

  if (mode.eq.1) then
     aString = 'x-'
  else if (mode.eq.2) then 
     aString = 'k-'
  else 
     write(*,*)"ERROR in get_correlations"
     stop
  endif   
!   assemble filename
  if ((NPar.gt.0) .and. (NPar.lt.10)) then 
    write (NParAsString1, '(I1)') NPar
    fullName=timeasString//'N'//trim(NParAsString1)//'M'//trim(MorbAsString)
    fullName=fullName//aString//slice1Asstring//'-'//slice2Asstring//'-'//'correlations.dat'
  else  if ((NPar.ge.10) .and. (NPar.lt.100)) then 
    write (NParAsString2, '(I2)') NPar
    fullName=timeasString//'N'//trim(NParAsString2)//'M'//trim(MorbAsString)
    fullName=fullName//aString//slice1Asstring//'-'//slice2Asstring//'-'//'correlations.dat'
  else if ((NPar.ge.100) .and. (NPar.lt.1000)) then 
    write (NParAsString3, '(I3)') NPar
    fullName=timeasString//'N'//trim(NParAsString3)//'M'//trim(MorbAsString)
    fullName=fullName//aString//slice1Asstring//'-'//slice2Asstring//'-'//'correlations.dat'
  else if ((NPar.ge.1000) .and. (NPar.lt.10000)) then 
    write (NParAsString4, '(I4)') NPar
    fullName=timeasString//'N'//trim(NParAsString4)//'M'//trim(MorbAsString)
     fullName=fullName//aString//slice1Asstring//'-'//slice2Asstring//'-'//'correlations.dat'
  else if ((NPar.ge.10000) .and. (NPar.lt.100000)) then 
    write (NParAsString5, '(I5)') NPar
    fullName=timeasString//'N'//trim(NParAsString5)//'M'//trim(MorbAsString)
    fullName=fullName//aString//slice1Asstring//'-'//slice2Asstring//'-'//'correlations.dat'

  else if ((NPar.ge.100000) .and. (NPar.lt.1000000)) then 
    write (NParAsString6, '(I6)') NPar
    fullName=timeasString//'N'//trim(NParAsString6)//'M'//trim(MorbAsString)
    fullName=fullName//aString//slice1Asstring//'-'//slice2Asstring//'-'//'correlations.dat'

  else if ((NPar.ge.1000000) .and. (NPar.lt.10000000)) then 
    write (NParAsString7, '(I7)') NPar
    fullName=timeasString//'N'//trim(NParAsString7)//'M'//trim(MorbAsString)//aString
    fullName=fullName//slice1Asstring//'-'//slice2Asstring//'-'//'correlations.dat'

  else if ((NPar.ge.10000000) .and. (NPar.lt.100000000)) then 
    write (NParAsString8, '(I8)') NPar
    fullName=timeasString//'N'//trim(NParAsString8)//'M'//trim(MorbAsString)//aString
     fullName=fullName//slice1Asstring//'-'//slice2Asstring//'-'//'correlations.dat'

  else 
    NParAsString4='XXXX'
    fullName=timeasString//'N'//trim(NParAsString4)//'M'//trim(MorbAsString)//aString
   fullName=fullName//slice1Asstring//'-'//slice2Asstring//'-'//'correlations.dat'
  end if
  call get_maximum(PSI,rho_jk,maximum)
  threshold = Kdip * DBLE(maximum)/w**2 ! set the minimal density to be plotted
  
  open(unit=12,file=trim(fullName),form='formatted')
  do n=1,NOGRIDPOINTS
!compute one body density onebdens 
     G1n_n = (0.d0,0.d0)
     do Jorb = 1,Morb  
        do Korb=1,Morb
           G1n_n = G1n_n + rho_jk(Jorb,Korb) * DCONJG(PSI(n,Jorb))*PSI(n,Korb)
        end do
     end do
     G1n_n=G1n_n*1.d0/w**2 
     onebdens(n)=DBLE(G1n_n)
!---------------
   end do

! determine which points are to be plotted    
   x1done=.FALSE.
   y1done=.FALSE.
   x2done=.FALSE.
   y2done=.FALSE.


   do n=1,NOGRIDPOINTS
      call get_ijk_from_m(n,NDX,NDY,i_n,j_n,k_n)

 IF (x1done.eqv..FALSE.) then
!     x1slice
      if ((mode.eq.1).and.((abs(ort_X(i_n)-x1slice).lt.ort_X(2)-ort_X(1)).and.&
           (x1const.eqv..true.))&
     .or.((mode.eq.2).and.(abs(mom_x(i_n)-x1slice).lt.&
          mom_X(2)-mom_X(1)).and.(x1const.eqv..true.))) then
          call get_ijk_from_m(n+1,NDX,NDY,i_m,j_m,k_m) 
          if ((abs(ort_X(i_n)-x1slice)).lt.(abs(ort_X(i_m)-x1slice))) then
             slice1pt=i_n
                write(*,*) 'I assigned 1 slice1',i_n
             x1done=.TRUE.
          else
             slice1pt=i_m
                write(*,*) 'I assigned 2 slice1',i_m
             x1done=.TRUE.
          endif
      endif
 ENDIF


 IF (y1done.eqv..FALSE.) then
!     y1slice

      if (((mode.eq.1).and.(abs(ort_Y(j_n)-y1slice).lt.&
        ort_Y(2)-ort_Y(1)).and.(y1const.eqv..true.)).or.&
      ((mode.eq.2).and.(abs(mom_Y(j_n)-y1slice).lt.&
        mom_Y(2)-mom_Y(1)).and.(y1const.eqv..true.))) then
          call get_ijk_from_m(n+1,NDX,NDY,i_m,j_m,k_m) 
          if ((abs(ort_Y(j_n)-y1slice)).lt.(abs(ort_Y(j_m)-y1slice))) then
             if (x1const.eqv..false.)  then
                slice1pt=j_n
                write(*,*) 'I assigned 3 slice1',j_n
                y1done=.TRUE.
             endif
          else
             if (x1const.eqv..false.) slice1pt=j_m
             if  (x1const.eqv..true.) then 
                slice2pt=j_m
                write(*,*) 'I assigned 4 slice2',j_m
                y1done=.TRUE.
                goto 8855
             endif
          endif
      endif
 ENDIF
 IF (x2done.eqv..FALSE.) then
!     x2slice
      if (((mode.eq.1).and.(abs(ort_X(i_n)-x2slice).lt.&
         ort_x(2)-ort_x(1)).and.(x2const.eqv..true.)).or.&
       ((mode.eq.2).and.(abs(mom_X(i_n)-x2slice).lt.mom_x(2)-mom_x(1)).and.&
            (x2const.eqv..true.))) then
          call get_ijk_from_m(n+1,NDX,NDY,i_m,j_m,k_m) 
          if ((abs(ort_x(i_n)-x2slice)).lt.(abs(ort_X(i_m)-x2slice))) then
             if ((x1const.eqv..false.).and.(y1const.eqv..false.)) then
                    slice1pt=i_n
                    write(*,*) 'I assigned 5 slice2=',i_n
             endif
             if ((x1const.eqv..true.).or.(y1const.eqv..true.)) then
                 slice2pt=i_n
                 write(*,*) 'I assigned 6 slice2=',i_n
                 x2done=.TRUE.
                 goto 8855
             endif
          else
             if ((x1const.eqv..false.).and.(y1const.eqv..false.)) then 
                 slice1pt=i_m
                 write(*,*) 'I assigned 7 slice2=',i_m
             endif
             if  ((x1const.eqv..true.).or.(y1const.eqv..true.)) then 
                 slice2pt=i_m
                 write(*,*) 'I assigned 8 slice2=',i_m
                 x2done=.TRUE.
                goto 8855
             endif
          endif
      endif
 ENDIF


 IF (y2done.eqv..FALSE.) then
!     y2slice
      if (((mode.eq.1).and.(abs(ort_Y(j_n)-y2slice).lt.&
        ort_Y(2)-ort_Y(1)).and.(y2const.eqv..true.)).or.&
      ((mode.eq.2).and.(abs(mom_Y(j_n)-y2slice).lt.&
        mom_Y(2)-mom_Y(1)).and.(y2const.eqv..true.))) then 
          call get_ijk_from_m(n+1,NDX,NDY,i_m,j_m,k_m) 
          if ((abs(ort_Y(j_n)-y2slice)).lt.(abs(ort_Y(j_m)-y2slice))) then
             if ((x1const.eqv..true.).or.(y1const.eqv..true.).or.((x2const.eqv..true.))) then
                 slice2pt=j_n
                 y2done=.TRUE.
                 write(*,*) 'I assigned 9 slice2=',j_n
                 goto 8855
             else 
                 write(*,*) 'Something went wrong with the assignment of the slices!!!'
                 stop
             endif
          else
             if ((x1const.eqv..true.).or.(y1const.eqv..true.).or.((x2const.eqv..true.))) then
                slice2pt=j_m
                 write(*,*) 'I assigned 10 slice2=',j_m
                y2done=.TRUE.
                goto 8855
             else 
                 write(*,*) 'Something went wrong with the assignment of the slices!!!'
                 stop
             endif
          endif
       endif
 ENDIF
 8855 continue
   end do
! 8855 continue
 do n=1,NOGRIDPOINTS
   call get_ijk_from_m(n,NDX,NDY,i_n,j_n,k_n)
!  x1const and y1const
   if (((x1const.eqv..true.).and.(i_n.eq.slice1pt)).and.&
      ((y1const.eqv..true.).and.(j_n.eq.slice2pt))) then
   do m=1,NOGRIDPOINTS
      call get_ijk_from_m(m,NDX,NDY,i_m,j_m,k_m)

!compute G1m_n
        G1m_n = (0.d0,0.d0) 
        do Jorb = 1,Morb  
           do Korb=1,Morb
              G1m_n = G1m_n + rho_jk(Jorb,Korb) * DCONJG(PSI(m,Jorb))*PSI(n,Korb)
           end do
        end do
        G1m_n=G1m_n*1.d0/(w**2) 
        outReG1m_n=DBLE(G1m_n)   
        outImG1m_n=DIMAG(G1m_n)
        if (m.eq.n) integ=integ+G1m_n
!----------------


!compute G2m_n 
        G2m_n=(0.d0,0.d0)
        do Lorb=1,Morb
           do Korb=1,Morb
              do Jorb=1,Morb
                 do Iorb=1,Morb
                    G2m_n = G2m_n + rho_ijkl(Iorb,Jorb,Korb,Lorb)* &
                              DCONJG(PSI(m,Iorb))*DCONJG(PSI(n,Jorb))*PSI(m,Korb)*PSI(n,Lorb)
  
                 end do
              end do
           end do
        end do
        G2m_n=G2m_n*1.d0/(w*w*w*w) 
        integ2=integ2+G2m_n
        outG2m_n=DBLE(G2m_n)
!----------------

        if (mode.eq.1) then
           write(12,6666) ort_X(i_m),ort_Y(j_m),ort_X(i_n),ort_Y(j_n), &
                       onebdens(m),onebdens(n),outReG1m_n,outImG1m_n,outG2m_n
        else if (mode.eq.2) then
           write(12,6666) mom_X(i_m),mom_Y(j_m),mom_X(i_n),mom_Y(j_n), &
                       onebdens(m),onebdens(n),outReG1m_n,outImG1m_n,outG2m_n
        end if
     if (mod(m,ndx).eq.0)     write(12,*)"                                                       " 

     end do
! if x1const and x2const or y2const
     else if ((x1const.eqv..true.).and.((x2const.eqv..true.).or.(y2const.eqv..true.))&
              .and.(i_n.eq.slice1pt)) then 
     do m=1,NOGRIDPOINTS

        call get_ijk_from_m(m,NDX,NDY,i_m,j_m,k_m)
! if the point is in the slice
        if (((x2const.eqv..true.).and.(i_m.eq.slice2pt))&
         .or.((y2const.eqv..true.).and.(j_m.eq.slice2pt))) then
!compute G1m_n
        G1m_n = (0.d0,0.d0) 
        do Jorb = 1,Morb  
           do Korb=1,Morb
              G1m_n = G1m_n + rho_jk(Jorb,Korb) * DCONJG(PSI(m,Jorb))*PSI(n,Korb)
           end do
        end do
        G1m_n=G1m_n*1.d0/(w**2) 
        outReG1m_n=DBLE(G1m_n)   
        outImG1m_n=DIMAG(G1m_n)
        if (m.eq.n) integ=integ+G1m_n
!----------------


!compute G2m_n 
        G2m_n=(0.d0,0.d0)
        do Lorb=1,Morb
           do Korb=1,Morb
              do Jorb=1,Morb
                 do Iorb=1,Morb
                    G2m_n = G2m_n + rho_ijkl(Iorb,Jorb,Korb,Lorb)* &
                              DCONJG(PSI(m,Iorb))*DCONJG(PSI(n,Jorb))*PSI(m,Korb)*PSI(n,Lorb)
  
                 end do
              end do
           end do
        end do
        G2m_n=G2m_n*1.d0/(w*w*w*w) 
        integ2=integ2+G2m_n
        outG2m_n=DBLE(G2m_n)
!----------------

        if (mode.eq.1) then
           write(12,6666) ort_X(i_m),ort_Y(j_m),ort_X(i_n),ort_Y(j_n), &
                       onebdens(m),onebdens(n),outReG1m_n,outImG1m_n,outG2m_n
        else if (mode.eq.2) then
           write(12,6666) mom_X(i_m),mom_Y(j_m),mom_X(i_n),mom_Y(j_n), &
                       onebdens(m),onebdens(n),outReG1m_n,outImG1m_n,outG2m_n
        end if

        endif ! loops only executed if point is in the slice
     enddo
     write(12,*) "                                                " 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  if x1const false, y1const true and x2/y2const true
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     elseif ((x1const.eqv..false.).and.(y1const.eqv..true.).and.(j_n.eq.slice1pt)) then

     do m=1,NOGRIDPOINTS
        call get_ijk_from_m(m,NDX,NDY,i_m,j_m,k_m)
! if the point is in the slice
        if (((x2const.eqv..true.).and.(i_m.eq.slice2pt))&
         .or.((y2const.eqv..true.).and.(j_m.eq.slice2pt))) then
!compute G1m_n
        G1m_n = (0.d0,0.d0) 
        do Jorb = 1,Morb  
           do Korb=1,Morb
              G1m_n = G1m_n + rho_jk(Jorb,Korb) * DCONJG(PSI(m,Jorb))*PSI(n,Korb)
           end do
        end do
        G1m_n=G1m_n*1.d0/(w**2) 
        outReG1m_n=DBLE(G1m_n)   
        outImG1m_n=DIMAG(G1m_n)
        if (m.eq.n) integ=integ+G1m_n
!----------------


!compute G2m_n 
        G2m_n=(0.d0,0.d0)
        do Lorb=1,Morb
           do Korb=1,Morb
              do Jorb=1,Morb
                 do Iorb=1,Morb
                    G2m_n = G2m_n + rho_ijkl(Iorb,Jorb,Korb,Lorb)* &
                              DCONJG(PSI(m,Iorb))*DCONJG(PSI(n,Jorb))*PSI(m,Korb)*PSI(n,Lorb)
  
                 end do
              end do
           end do
        end do
        G2m_n=G2m_n*1.d0/(w*w*w*w) 
        integ2=integ2+G2m_n
        outG2m_n=DBLE(G2m_n)
!----------------

        if (mode.eq.1) then
           write(12,6666) ort_X(i_m),ort_Y(j_m),ort_X(i_n),ort_Y(j_n), &
                       onebdens(m),onebdens(n),outReG1m_n,outImG1m_n,outG2m_n
        else if (mode.eq.2) then
           write(12,6666) mom_X(i_m),mom_Y(j_m),mom_X(i_n),mom_Y(j_n), &
                       onebdens(m),onebdens(n),outReG1m_n,outImG1m_n,outG2m_n
        end if

        endif ! loops only executed if point is in the slice
     enddo
      write(12,*)"                                                 " 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  if x1const false, y1const false and x2const and y2const true
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     elseif ((x1const.eqv..false.).and.(y1const.eqv..false.)) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do m=1,NOGRIDPOINTS
        call get_ijk_from_m(m,NDX,NDY,i_m,j_m,k_m)
! if the point is in the slice
        if (((x2const.eqv..true.).and.(i_m.eq.slice1pt))&
         .and.((y2const.eqv..true.).and.(j_m.eq.slice2pt))) then
!compute G1m_n
        G1m_n = (0.d0,0.d0) 
        do Jorb = 1,Morb  
           do Korb=1,Morb
              G1m_n = G1m_n + rho_jk(Jorb,Korb) * DCONJG(PSI(m,Jorb))*PSI(n,Korb)
           end do
        end do
        G1m_n=G1m_n*1.d0/(w**2) 
        outReG1m_n=DBLE(G1m_n)   
        outImG1m_n=DIMAG(G1m_n)
        if (m.eq.n) integ=integ+G1m_n
!----------------


!compute G2m_n 
        G2m_n=(0.d0,0.d0)
        do Lorb=1,Morb
           do Korb=1,Morb
              do Jorb=1,Morb
                 do Iorb=1,Morb
                    G2m_n = G2m_n + rho_ijkl(Iorb,Jorb,Korb,Lorb)* &
                              DCONJG(PSI(m,Iorb))*DCONJG(PSI(n,Jorb))*PSI(m,Korb)*PSI(n,Lorb)
  
                 end do
              end do
           end do
        end do
        G2m_n=G2m_n*1.d0/(w*w*w*w) 
        integ2=integ2+G2m_n
        outG2m_n=DBLE(G2m_n)
!----------------
        if (mode.eq.1) then
           write(12,6666) ort_X(i_m),ort_Y(j_m),ort_X(i_n),ort_Y(j_n), &
                       onebdens(m),onebdens(n),outReG1m_n,outImG1m_n,outG2m_n
           if (mod(i_n,NDX).eq.0) then 
              write(12,*) "                                                                   "
           endif
        else if (mode.eq.2) then

           write(12,6666) mom_X(i_m),mom_Y(j_m),mom_X(i_n),mom_Y(j_n), &
                       onebdens(m),onebdens(n),outReG1m_n,outImG1m_n,outG2m_n
           if (mod(j_n,NDY).eq.0) then 
              write(12,*) "                                                                   "
           endif


        end if

        endif ! loops only executed if point is in the slice
     enddo
     endif !loop m only if n is in the slice
  end do         

 6666 FORMAT(9(E20.10))

  close(12)
 
 
!---}}}
end subroutine get_corr_slice

!subroutine MKL_FT(PSI,FTPSI) 
!!---{{{
!USE SHARED_DIMS, ONLY : DIM_MCTDHB
!use DVR_ALL
!use W_INTERPARTICLE
!Use MKL_DFTI
!  implicit none
!  integer                  :: Iorb,m,n,i,j,k
!  integer                  :: jmin,jzero,jminone,jmax 
!  integer                  :: imin,izero,iminone,imax 
!
!  complex*16,DIMENSION(NDX*NDY*NDZ,Morb),intent(in)  :: PSI 
!  complex*16,DIMENSION(NDX*NDY*NDZ,Morb),intent(out) :: FTPSI
!  complex*16,DIMENSION(NDX*NDY*NDZ) :: ORBTMP 
! 
!  real*8,dimension(NDX*NDY*NDZ) :: kvecs
!  real*8                        :: xall,xfirst,xlast,dk,dx,kmin,norm,xi,xf
!
!  type(DFTI_DESCRIPTOR), POINTER :: Desc_Handle_Dim1
!  
!  Integer :: Status,STRIDE(2),L(2)
!
!  if (NDZ.gt.1) then
!     write(*,*) "more than 2-D FFT has not been implemented yet"
!     return
!  end if
!  NORM=1.d0
!!##########################################################  
!!##########################################################  
!!     1D -CASE
!    if (DIM_MCTDHB.eq.1) then
!!##########################################################  
!!##########################################################  
!   FTPSI=PSI
!      
!!  do the mkl fft voodoo
!
!   Status = DftiCreateDescriptor(Desc_Handle_Dim1, DFTI_DOUBLE,&
!    DFTI_COMPLEX, 1, NDX )
!   Status = DftiSetValue( Desc_Handle_Dim1,&
!    DFTI_NUMBER_OF_TRANSFORMS, NDY*NDZ )
!   Status = DftiSetValue( Desc_Handle_Dim1, DFTI_INPUT_DISTANCE, NDX)
!   Status = DftiSetValue( Desc_Handle_Dim1, DFTI_OUTPUT_DISTANCE,NDX)
!   do Iorb=1,Morb
!         ORBTMP=(FTPSI(:,IORB))
!         Status = DftiCommitDescriptor(Desc_Handle_Dim1)
!         Status = DftiComputeForward(Desc_Handle_Dim1,ORBTMP)
!
!         if (mod(NDX,2).eq.0) then
!           do m=(NDX/2+1),NDX
!              FTPSI(m-(NDX/2),IORB)=ORBTMP(m)
!              FTPSI(m,IORB)=ORBTMP(m-(NDX/2))
!           enddo
!         else
!           do m=(NDX/2+1),NDX
!              FTPSI(m-(NDX/2),IORB)=ORBTMP(m)
!              FTPSI(m,IORB)=ORBTMP(m-(NDX/2))
!           enddo
!         endif
!
!     call normvxz(FTPSI(:,Iorb),norm,NDX)
!!         write(6,*)"Norm ",norm
!     call xvixdzo(norm,FTPSI(:,Iorb),NDX)

 !        do m=1,NDX
 !           FTPSI(m,IORB)=FTPSI(m,IORB)*(-1.d0)**(m-1)
 !        end do 
!  end do
!!##########################################################  
!!##########################################################  
!!     1D -CASE
!    endif 
!##########################################################  
!##########################################################  

!##########################################################  
!##########################################################  
!     2D -CASE
!    if (DIM_MCTDHB.eq.2) then
!!##########################################################  
!!##########################################################  
!     FTPSI=PSI
!     L(1)=NDX
!     L(2)=NDY
!     Status = DftiCreateDescriptor(Desc_Handle_Dim1, DFTI_DOUBLE,&
!           DFTI_COMPLEX, 2, L )
!
!     do Iorb=1,Morb
!
!         ORBTMP=(FTPSI(:,IORB))
!      
!         Status = DftiCommitDescriptor(Desc_Handle_Dim1)
!         Status = DftiComputeForward(Desc_Handle_Dim1,ORBTMP)
! 
!         call normvxz(ORBTMP,norm,NDX*NDY)
!         call xvixdzo(norm,ORBTMP,NDX*NDY)
!         
!         call sort_FFT_to_ascending_2d(ORBTMP,NDX,NDY)
!
!         FTPSI(:,IORB)=ORBTMP
!     end do
!!##########################################################  
!!##########################################################  
!!     2D -CASE
!    endif
!!##########################################################  
!!##########################################################  
!!---}}}
!end subroutine MKL_FT
!

subroutine sort_FFT_to_ascending_2d(FFT,NDX,NDY)
! sort an array representing a 2D object given in FFT order, compact storage 
! to be in strictly ascending order (needed for gnuplots pm3d), compact storage

IMPLICIT NONE
integer NDX,NDY,i,j,k
COMPLEX*16 FFT(NDX*NDY),ASC(NDX*NDY)

         ASC=FFT
 
! copy the first half of each vector in X direction to the corresponding last half

   IF (MOD(NDX,2).eq.0) THEN
         do k=1,NDY
            do i=1,NDX
             if (i.le.((NDX)/2)) then
                FFT((k-1)*NDX+i)=ASC((k-1)*NDX+i+(NDX)/2)
             else if (i.gt.((NDX)/2)) then
                FFT((k-1)*NDX+i)=ASC((k-1)*NDX+i-((NDX)/2))
             endif
            end do
         end do 
   ENDIF


! copy the first half of each vector in X direction to the corresponding last half
   IF (MOD(NDX,2).ne.0) then
         do k=1,NDY
            do i=1,NDX
             if (i.le.((NDX+1)/2-1)) then
                FFT((k-1)*NDX+i)=ASC((k-1)*NDX+i+(NDX+1)/2)
             else if (i.gt.((NDX+1)/2)-1) then
                FFT((k-1)*NDX+i)=ASC((k-1)*NDX+i-((NDX+1)/2-1))
             endif
            end do
         end do 
   ENDIF

   IF (MOD(NDY,2).eq.0) then
! copy y direction from 0 -> kmax to last half

         do k=1,NDX*NDY
            if (k.le.(NDX*(NDY/2))) then
               ASC(k)=FFT(NDX*(NDY/2)+k)
            elseif (k.gt.NDX*((NDY)/2)) then
               ASC(k)=FFT(k-NDX*(NDY/2))
            endif
         end do
   ENDIF
   IF (MOD(NDY,2).ne.0) then
! copy y direction from 0 -> kmax to last half
         do k=1,NDX*NDY
            if (k.le.(NDX*(NDY/2))) then
               ASC(k)=FFT(NDX*(NDY/2+1)+k)
            elseif (k.gt.NDX*(NDY/2)) then
               ASC(k)=FFT(k-NDX*((NDY/2)))
            endif
         end do
   ENDIF
   FFT=ASC


end subroutine sort_FFT_to_ascending_2d

subroutine get_1BodyDiag(time,PSI,rho_jk,mode,dilation,pathPRN) 
!--------------------------------------------------------------------------------------------
! computes G1(x,x) if PSI is provided or G1(k,k) if FT(PSI) is provided and
! writes it to a file. 
!
! time          : current time 
! PSI           : orbital array 
! rho_jk        : the 1-body density matrix where rho_jk(j,k)=rho_jk 
! mode          : mode=1: G1(x,x) 
!                 mode=2: G1(k,k) is computed
! dilation      : the integer dilation gives the factor by which the momentum intervall should
!                 be shrunk or equivalently the factor by which the sampling rate is increased.
!                 The momentum distribution at the outer ends is padded with zeros then. 
!                 Normally dilation is computed by get_dilation. If unsure put
!                 dilation = 1, but the momentum grid will then usually have a bad
!                 resolution. 
! pathPRN       :  character*10 path,pathCI,pathORB,pathORBK 
!              pathORB='DATA/orb_R'
!              pathCI='DATA/CIcnf'
!              pathORBK='DATA/orb_K'
!{ {{
!--------------------------------------------------------------------------------------------
use CI_ALL,ONLY:NPar
use DVR_ALL
USE W_INTERPARTICLE

  implicit none
   character*10 pathPRN
!  integer,parameter  :: NOGRIDPOINTS=NDX*NDY*NDZ
  integer  :: NOGRIDPOINTS
  integer,intent(in) :: mode,dilation
  INTEGER :: dil
  real*8             :: time,w,dk,mom_X(NDX),mom_Y(NDY),mom_Z(NDZ)
!  complex*16,DIMENSION(NOGRIDPOINTS,Morb),intent(in) :: PSI 
  COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb)  :: PSI
  COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb)  :: FTPSI
  complex*16,DIMENSION(Morb,Morb),intent(in) :: rho_jk
  character(len=2)   :: aString
  character(len=200) :: fullName
  character(len=10)  :: timeAsString
  character(len=10)  :: MorbAsString
  character(len=1)   :: NparAsString1
  character(len=2)   :: NparAsString2
  character(len=3)   :: NparAsString3
  character(len=4)   :: NparAsString4
  character(len=5)   :: NparAsString5
  character(len=6)   :: NparAsString6
  character(len=7)   :: NparAsString7
  character(len=8)   :: NparAsString8
  integer            :: n,Jorb,Korb,i_n,j_n,k_n
  complex*16         :: G1n_n
  complex*16         :: maximum
  real*8, parameter  :: TOL=0.00001d0 
  real*8             :: threshold,wdil,dkdil

  real*8             :: outG1n_n,integ2
  real*8             :: xi,xf,yi,yf,zi,zf


   xi=Time_xint
   xf=Time_xfnl
   yi=Time_yint
   yf=Time_yfnl
   zi=Time_zint
   zf=Time_zfnl
  NOGRIDPOINTS=NDX*NDY*NDZ
!  if((NDY.gt.1).or. (NDZ.gt.1)) then
!     write(*,*)"only 1d implemented so far"
!     stop
!  end if

   

  if(mode.eq.1) then
     w=DSQRT((xf-xi)/NDX)
  else if (mode.eq.2) then ! when 'Pedestrian FT' was used
     dk     = 2.d0*PI/(xf-xi)/dilation
     w      = DSQRT( dk )
     if (DIM_MCTDHB.ge.2) then 
        write(*,*) 'Pedestrian_FT only works in 1D (onebody density)'
        return
     endif
     CALL get_mom(xi,xf,NDX,mom_X) !returns momentum grid in strictly ascending order
     CALL get_mom(yi,yf,NDY,mom_Y)
     CALL get_mom(zi,zf,NDZ,mom_Z)

     mom_X=mom_X/dilation
     mom_Y=mom_Y/dilation
     mom_Z=mom_Z/dilation

  else if (mode.eq.3) then ! when MKL_FT was used
     CALL get_mom(xi,xf,NDX,mom_X) !returns momentum grid in strictly ascending order
     CALL get_mom(yi,yf,NDY,mom_Y) !returns momentum grid in strictly ascending order
     CALL get_mom(zi,zf,NDZ,mom_Z) !returns momentum grid in strictly ascending order
     if (DIM_MCTDHB.eq.1) then
       dk     = mom_X(2)-mom_X(1) 
     elseif (DIM_MCTDHB.eq.2) then
       dk     = (mom_Y(2)-mom_Y(1))*(mom_X(2)-mom_X(1))
     else if (DIM_MCTDHB.eq.3) then
       dk     =(mom_z(2)-mom_z(1))*(mom_Y(2)-mom_Y(1))*(mom_X(2)-mom_X(1))
!      write(*,*) 'one body momentum density is not implemented for 3D so far' 
!      return
     endif
     w      = DSQRT( dk )
  else 
        write(*,*) "wrong mode"
  end if
  
  if (MORB.lt.10) then
     write (MorbAsString, '(I1)') Morb
  else if ((MORB.ge.10).and.(MORB.lt.100)) then
     write (MorbAsString, '(I2)') Morb
  else 
     write(*,*) 'Bigger orbital number than 100 not implemented!'
  endif
  call getDoubleAsString(time,timeAsString) 

  if (mode.eq.1) then
     aString = 'x-'
  else if ((mode.eq.2).or.(mode.eq.3)) then 
     aString = 'k-'
  else 
     write(*,*)"ERROR in get_correlations"
     stop
  endif   
  if ((NPar.gt.0) .and. (NPar.lt.10)) then 
    write (NParAsString1, '(I1)') NPar
    fullName=timeasString//'N'//trim(NParAsString1)//'M'//trim(MorbAsString)//aString//'density.dat'
  else  if ((NPar.ge.10) .and. (NPar.lt.100)) then 
    write (NParAsString2, '(I2)') NPar
    fullName=timeasString//'N'//trim(NParAsString2)//'M'//trim(MorbAsString)//aString//'density.dat'
  else if ((NPar.ge.100) .and. (NPar.lt.1000)) then 
    write (NParAsString3, '(I3)') NPar
    fullName=timeasString//'N'//trim(NParAsString3)//'M'//trim(MorbAsString)//aString//'density.dat'
  else if ((NPar.ge.1000) .and. (NPar.lt.10000)) then 
    write (NParAsString4, '(I4)') NPar
    fullName=timeasString//'N'//trim(NParAsString4)//'M'//trim(MorbAsString)//aString//'density.dat'
  else if ((NPar.ge.10000) .and. (NPar.lt.100000)) then 
    write (NParAsString5, '(I5)') NPar
    fullName=timeasString//'N'//trim(NParAsString5)//'M'//trim(MorbAsString)//aString//'density.dat'
  else if ((NPar.ge.100000) .and. (NPar.lt.1000000)) then 
    write (NParAsString6, '(I6)') NPar
    fullName=timeasString//'N'//trim(NParAsString6)//'M'//trim(MorbAsString)//aString//'density.dat'
  else if ((NPar.ge.1000000) .and. (NPar.lt.10000000)) then 
    write (NParAsString7, '(I7)') NPar
    fullName=timeasString//'N'//trim(NParAsString7)//'M'//trim(MorbAsString)//aString//'density.dat'
  else if ((NPar.ge.10000000) .and. (NPar.lt.100000000)) then 
    write (NParAsString8, '(I8)') NPar
    fullName=timeasString//'N'//trim(NParAsString8)//'M'//trim(MorbAsString)//aString//'density.dat'
  else 
    NParAsString4='XXXX'
    fullName=timeasString//'N'//trim(NParAsString4)//'M'//trim(MorbAsString)//aString//'density.dat'
  end if
    

    fullName=timeasString//aString//'density.dat'
              integ2=0d0
  open(unit=12,file=trim(pathPRN)//'/'//trim(fullName),form='formatted')
  do n=1,NOGRIDPOINTS

!compute G1n_n
     G1n_n = (0.d0,0.d0)
     do Jorb = 1,Morb  
        do Korb=1,Morb

           G1n_n = G1n_n + rho_jk(Jorb,Korb) * DCONJG(PSI(n,Jorb))*PSI(n,Korb)

        end do
     end do
     call get_ijk_from_m(n,NDX,NDY,i_n,j_n,k_n)
     if (mode.eq.1)  w=weight(n)
     G1n_n=G1n_n*1.d0/w**(2)
     outG1n_n=DBLE(G1n_n)/Npar
      integ2=integ2+outG1n_n*w**2
     if (mode.eq.1) then

         write(12,6667) ort_X(i_n),ort_Y(j_n),ort_Z(k_n),outG1n_n,time

     else if ((mode.eq.2).or.(mode.eq.3)) then
    
!         write(12,6667) mom_X(i_n),mom_Y(j_n),mom_Z(k_n),outG1n_n,time
         write(12,6668) mom_X(i_n),mom_Y(j_n),mom_Z(k_n),outG1n_n,time,(DBLE(DCONJG(PSI(n,Jorb))*PSI(n,Jorb)/w/w), Jorb=1,Morb)

     end if
     if (mod(i_n,NDX).eq.0) write(12,*) '                             '
  end do
                             
  close(12)
  write(6,'(a28,F12.6,a4,i2,a7,F8.4)') 'rho: Intg DNS(R[mode=1]|K)=',integ2,' mode ',mode,' T=',time


6667  format(3(F12.6),2(F18.9))
6668  format(3(F12.6),2(F18.9),25(F18.9))
!---}}}
end subroutine get_1BodyDiag

subroutine get_2D_1BodyDiag_dilated(time,rho_jk,PSI,dil)

  use CI_ALL,ONLY:NPar
  use DVR_ALL
  USE W_INTERPARTICLE

  INTEGER,INTENT(IN) ::  dil
  real*8,INTENT(IN)  ::  time
  COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb),INTENT(IN)  :: PSI
  complex*16,DIMENSION(Morb,Morb),intent(in) :: rho_jk
  COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb)             :: FTPSI
  character(len=2)   :: aString
  character(len=200) :: fullName
  character(len=10)  :: timeAsString
  character(len=10)  :: MorbAsString
  character(len=1)   :: NparAsString1
  character(len=2)   :: NparAsString2
  character(len=3)   :: NparAsString3
  character(len=4)   :: NparAsString4
  character(len=5)   :: NparAsString5
  character(len=6)   :: NparAsString6
  character(len=7)   :: NparAsString7
  character(len=8)   :: NparAsString8

  real*8             :: dilmom_X(NDX),dilmom_Y(NDY),dkdil,wdil,outG1n_n
  real*8             :: xi,xf,yi,yf,zi,zf

  complex*16         :: G1n_n
  INTEGER            :: i_n,j_n,k_n,JORB,KORB,n
   xi=Time_xint
   xf=Time_xfnl
   yi=Time_yint
   yf=Time_yfnl


  if (MORB.lt.10) then
     write (MorbAsString, '(I1)') Morb
  else if ((MORB.ge.10).and.(MORB.lt.100)) then
     write (MorbAsString, '(I2)') Morb
  else 
     write(*,*) 'Bigger orbital number than 100 not implemented!'
  endif
  call getDoubleAsString(time,timeAsString) 

  aString = 'k-'
  write(*,*)"ERROR in get_correlations"
  if ((NPar.gt.0) .and. (NPar.lt.10)) then 
    write (NParAsString1, '(I1)') NPar
    fullName=timeasString//'N'//trim(NParAsString1)//'M'//trim(MorbAsString)//aString//'density.dat'
  else  if ((NPar.ge.10) .and. (NPar.lt.100)) then 
    write (NParAsString2, '(I2)') NPar
    fullName=timeasString//'N'//trim(NParAsString2)//'M'//trim(MorbAsString)//aString//'density.dat'
  else if ((NPar.ge.100) .and. (NPar.lt.1000)) then 
    write (NParAsString3, '(I3)') NPar
    fullName=timeasString//'N'//trim(NParAsString3)//'M'//trim(MorbAsString)//aString//'density.dat'
  else if ((NPar.ge.1000) .and. (NPar.lt.10000)) then 
    write (NParAsString4, '(I4)') NPar
    fullName=timeasString//'N'//trim(NParAsString4)//'M'//trim(MorbAsString)//aString//'density.dat'
  else if ((NPar.ge.10000) .and. (NPar.lt.100000)) then 
    write (NParAsString5, '(I5)') NPar
    fullName=timeasString//'N'//trim(NParAsString5)//'M'//trim(MorbAsString)//aString//'density.dat'
  else if ((NPar.ge.100000) .and. (NPar.lt.1000000)) then 
    write (NParAsString6, '(I6)') NPar
    fullName=timeasString//'N'//trim(NParAsString6)//'M'//trim(MorbAsString)//aString//'density.dat'
  else if ((NPar.ge.1000000) .and. (NPar.lt.10000000)) then 
    write (NParAsString7, '(I7)') NPar
    fullName=timeasString//'N'//trim(NParAsString7)//'M'//trim(MorbAsString)//aString//'density.dat'
  else if ((NPar.ge.10000000) .and. (NPar.lt.100000000)) then 
    write (NParAsString8, '(I8)') NPar
    fullName=timeasString//'N'//trim(NParAsString8)//'M'//trim(MorbAsString)//aString//'density.dat'
  else 
    NParAsString4='XXXX'
    fullName=timeasString//'N'//trim(NParAsString4)//'M'//trim(MorbAsString)//aString//'density.dat'
  end if

  open(unit=12,file=trim(fullName),form='formatted')

  CALL get_mom(xi,xf,NDX,dilmom_X) !returns momentum grid in strictly ascending order
  CALL get_mom(yi,yf,NDY,dilmom_Y)

  dilmom_X=dilmom_X/dil
  dilmom_Y=dilmom_Y/dil
  dkdil     = 2.d0*PI/((xf-xi)*dil)*2.d0*PI/((yf-yi)*dil)
  wdil      = sqrt( dkdil )

!  CALL twoD_zero_padding(PSI,FTPSI,dil)

  do n=1,NDX*NDY
!compute G1n_n
     G1n_n = (0.d0,0.d0)
       do Jorb = 1,Morb  
          do Korb=1,Morb

             G1n_n = G1n_n + rho_jk(Jorb,Korb) * DCONJG(FTPSI(n,Jorb))*FTPSI(n,Korb)
  
          end do
       end do
       G1n_n=G1n_n*1.d0/wdil**2
       outG1n_n=DBLE(G1n_n)
       call get_ijk_from_m(n,NDX,NDY,i_n,j_n,k_n)
       write(12,'(E25.16,E25.16,E25.16)') dilmom_X(i_n),dilmom_Y(j_n),outG1n_n
       if (mod(i_n,(NDX)).eq.0) write(12,*) '                             '
       
    end do   
!       return
    close(12)
end subroutine get_2D_1BodyDiag_dilated


subroutine getsliceAsString(d,x,doubleString)
!       {{{ writes a double in 4 digit format,plus 
!           usefull  for naming files without blanks       
IMPLICIT NONE
integer       :: d
integer,parameter         :: DBL=8
real(kind=DBL)            :: x
CHARACTER(len=7)         :: doubleString
CHARACTER(len=2)          :: slc

if (d.eq.1) then 
    write(slc,'(A2)') 'c1'
elseif (d.eq.2) then 
    write(slc,'(A2)') 'c2'
elseif (d.eq.3)  then 
    write(slc,'(A2)') 'c3'
elseif (d.eq.4) then 
    write(slc,'(A2)') 'c4'
else 
  write(*,*) 'This slice cannot be right!'
  stop
endif

IF (x.ge.0.d0) then
    if (x .LE. 999.999999d0) WRITE(doubleString,'(A1,F4.0,A2)') 'P',x,slc
    if (x .LE. 99.9999999d0) WRITE(doubleString,'(A1,F4.1,A2)') 'P',x,slc
    if (x .LE. 9.99999999d0) WRITE(doubleString,'(A1,F4.2,A2)') 'P',x,slc
ELSE IF ((x.lt.0.d0).and.(abs(x).lt.1000.d0)) then
    if (abs(x) .LE. 999.999999d0) WRITE(doubleString,'(A1,F4.0,A2)') 'M',abs(x),slc
    if (abs(x) .LE. 99.9999999d0) WRITE(doubleString,'(A1,F4.1,A2)') 'M',abs(x),slc
    if (abs(x) .LE. 9.99999999d0) WRITE(doubleString,'(A1,F4.2,A2)') 'M',abs(x),slc
ELSE
  WRITE(6,*)"get slice as string: x not implemented!"
  WRITE(doubleString,'(A7)') 'XXXXXXX' 
END IF

!}}}
end subroutine getsliceasString

SUBROUTINE Get_VtrapProjection(PSI,Vtrap,Rho_IJ,time,dir)
  USE SHARED_DIMS
  USE CI_ALL
  USE DVR_ALL
  USE W_INTERPARTICLE
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb)       :: PSI
  COMPLEX*16, DIMENSION(NDX*NDY*NDZ)  ::  Vtrap
  COMPLEX*16, ALLOCATABLE ::  Vtrap_NormX(:,:,:)
  COMPLEX*16, ALLOCATABLE ::  Vtrap_NormY(:,:,:)
  COMPLEX*16, DIMENSION(MOrb,MOrb) ::  Rho_IJ
  COMPLEX*16, ALLOCATABLE ::  VX(:,:,:)
  REAL*8, ALLOCATABLE  ::  VX_EFF(:)
  COMPLEX*16, ALLOCATABLE ::  VY(:,:,:)
  REAL*8, ALLOCATABLE  ::  VY_EFF(:)
  REAL*8  :: time 
  INTEGER :: I,J,K,L,M,N,O,ierr


  character(len=200) :: fullName
  character(len=10)  :: timeAsString
  character(len=10)  :: MorbAsString
  character(len=1)   :: dir 
  character(len=1)   :: NparAsString1
  character(len=2)   :: NparAsString2
  character(len=3)   :: NparAsString3
  character(len=4)   :: NparAsString4
  character(len=5)   :: NparAsString5
  character(len=6)   :: NparAsString6
  character(len=7)   :: NparAsString7
  character(len=8)   :: NparAsString8

  IF (DIR.eq.'B') then
     N=1
     O=2
  ELSE IF (DIR.eq.'X') then
     N=1
     O=1
  ELSE IF (DIR.eq.'Y') then
     N=2
     O=2
  ENDIF

  do M=N,O 

    IF (M.eq.1) DIR='X'

    IF (M.eq.2) DIR='Y'

 
    IF (DIR.eq.'X') then
 
          allocate(VX_EFF(NDX),stat=ierr)
          if(ierr /= 0) then
             write(*,*)"allocation error for VX_EFF"
          endif
          allocate(Vtrap_NormX(NDX,MORB,MORB),stat=ierr)
          if(ierr /= 0) then
             write(*,*)"allocation error for Vtrap_norm"
          endif
          allocate(VX(NDX,MORB,MORB),stat=ierr)
          if(ierr /= 0) then
             write(*,*)"allocation error for VX"
          endif
  
    ELSEIF (DIR.eq.'Y') then 
 
          allocate(VY_EFF(NDY),stat=ierr)
          if(ierr /= 0) then
             write(*,*)"allocation error for VX_EFF"
          endif
 
          allocate(Vtrap_NormY(NDY,MORB,MORB),stat=ierr)
          if(ierr /= 0) then
             write(*,*)"allocation error for Vtrap_norm"
          endif
          allocate(VY(NDY,MORB,MORB),stat=ierr)
          if(ierr /= 0) then
             write(*,*)"allocation error for VX"
          endif
    ENDIF
 
    write(6,*) "Vtrap projections "
    write(6,*) "V_ij(x)= int(psi*(x,y) V(x,y) psi(x,y))dy for 2D Only" 
    if(DIM_MCTDHB.ne.2) then
      write(6,*)"It does not work in 1D and 3D "
      return         
    endif
 
    VX=(0d0,0d0)
    Vtrap_NormX=(0.d0,0.d0)
    VX_EFF=(0.d0,0.d0)
    VY=(0d0,0d0)
    Vtrap_NormY=(0.d0,0.d0)
    VY_EFF=(0.d0,0.d0)
    if (dir.eq.'X') then
  
    DO L=1,MORB
      DO K=1,MORB
 
        DO J=1,NDY
          DO I=1,NDX
 
             VX(I,L,K)=VX(I,L,K)+Rho_IJ(L,K)*&
                                 Conjg(PSI(I+NDX*(J-1),L))*VTRAP_EXT(I+NDX*(J-1))&
                                 *PSI(I+NDX*(J-1),K)&
                                 /weight_Y(J)/weight_Y(J)!/weight(I+NDX*(J-1))/weight(I+NDX*(J-1))
             Vtrap_NormX(I,L,K)=Vtrap_NormX(I,L,K)+Conjg(PSI(I+NDX*(J-1),L))&
                                                *PSI(I+NDX*(J-1),K)/weight_Y(J)/weight_Y(J)

          ENDDO
 
        ENDDO
 
      ENDDO
    ENDDO
 
       
         Do I=1,NDX
         Do L=1,Morb
         Do K=1,Morb
         if(ABS(Vtrap_NormX(I,L,K)).ge.1.0d-16) then
            VX(I,L,K)=VX(I,L,K)/Vtrap_NormX(I,L,K)
         else
            VX(I,L,K)=VX(I,L,K)/1.0d-16
         endif
         ENDDO
         ENDDO
         ENDDO

     DO L=1,MORB
      DO K=1,MORB
 
        DO I=1,NDX
    
          VX_EFF(I)=VX_EFF(I)+ABS(VX(I,L,K))
        ENDDO
 
      ENDDO
    ENDDO
 
 
    else if (dir.eq.'Y') then
    DO L=1,MORB
      DO K=1,MORB
 
         DO J=1,NDY
         DO I=1,NDX

             VY(J,L,K)=VY(J,L,K)+Rho_IJ(L,K)*&
                                 Conjg(PSI(I+NDX*(J-1),L))*VTRAP_EXT(I+NDX*(J-1))&
                                 *PSI(I+NDX*(J-1),K)&
                                 /weight_X(I)/weight_X(I)!/weight(I+NDX*(J-1))/weight(I+NDX*(J-1))
             Vtrap_NormY(J,L,K)=Vtrap_NormY(J,L,K)+Conjg(PSI(I+NDX*(J-1),L))&
                                                *PSI(I+NDX*(J-1),K)/weight_X(I)/weight_X(I)
          ENDDO
 
        ENDDO
 
      ENDDO
    ENDDO
         

!    VY=VY/Vtrap_NormY
           Do I=1,NDY
           Do L=1,Morb
           Do K=1,Morb
           if(ABS(Vtrap_NormY(I,L,K)).ge.1.0d-16) then
              VY(I,L,K)=VY(I,L,K)/Vtrap_NormY(I,L,K)
           else
              VY(I,L,K)=1.0d+16
           endif
           ENDDO
           ENDDO
           ENDDO

     DO L=1,MORB
      DO K=1,MORB
 
        DO I=1,NDY
    
          VY_EFF(I)=VY_EFF(I)+ABS(VY(I,L,K))
        ENDDO
 
      ENDDO
    ENDDO
 
    endif
 
    if (MORB.lt.10) then
       write (MorbAsString, '(I1)') Morb
    else if ((MORB.ge.10).and.(MORB.lt.100)) then
       write (MorbAsString, '(I2)') Morb
    else 
       write(*,*) 'Bigger orbital number than 100 not implemented!'
    endif
    call getDoubleAsString(time,timeAsString) 
 
    if ((NPar.gt.0) .and. (NPar.lt.10)) then 
      write (NParAsString1, '(I1)') NPar
      fullName=timeasString//'N'//trim(NParAsString1)//'M'//trim(MorbAsString)//'Vtrap_proj_'//trim(dir)//'.dat'
    else  if ((NPar.ge.10) .and. (NPar.lt.100)) then 
      write (NParAsString2, '(I2)') NPar
      fullName=timeasString//'N'//trim(NParAsString2)//'M'//trim(MorbAsString)//'Vtrap_proj_'//trim(dir)//'.dat'
    else if ((NPar.ge.100) .and. (NPar.lt.1000)) then 
      write (NParAsString3, '(I3)') NPar
      fullName=timeasString//'N'//trim(NParAsString3)//'M'//trim(MorbAsString)//'Vtrap_proj_'//trim(dir)//'.dat'
    else if ((NPar.ge.1000) .and. (NPar.lt.10000)) then 
      write (NParAsString4, '(I4)') NPar
      fullName=timeasString//'N'//trim(NParAsString4)//'M'//trim(MorbAsString)//'Vtrap_proj_'//trim(dir)//'.dat'
    else if ((NPar.ge.10000) .and. (NPar.lt.100000)) then 
      write (NParAsString5, '(I5)') NPar
      fullName=timeasString//'N'//trim(NParAsString5)//'M'//trim(MorbAsString)//'Vtrap_proj_'//trim(dir)//'.dat'
    else if ((NPar.ge.100000) .and. (NPar.lt.1000000)) then 
      write (NParAsString6, '(I6)') NPar
      fullName=timeasString//'N'//trim(NParAsString6)//'M'//trim(MorbAsString)//'Vtrap_proj_'//trim(dir)//'.dat'
    else if ((NPar.ge.1000000) .and. (NPar.lt.10000000)) then 
      write (NParAsString7, '(I7)') NPar
      fullName=timeasString//'N'//trim(NParAsString7)//'M'//trim(MorbAsString)//'Vtrap_proj_'//trim(dir)//'.dat'
    else if ((NPar.ge.10000000) .and. (NPar.lt.100000000)) then 
      write (NParAsString8, '(I8)') NPar
      fullName=timeasString//'N'//trim(NParAsString8)//'M'//trim(MorbAsString)//'Vtrap_proj_'//trim(dir)//'.dat'
    else 
      NParAsString4='XXXX'
      fullName=timeasString//'N'//trim(NParAsString4)//'M'//trim(MorbAsString)//'Vtrap_proj_'//trim(dir)//'.dat'
    end if
 
 
    open(unit=120,file=trim(fullName),form='formatted')
    if (DIR.eq.'X') then 
    DO I=1,NDX
      write(120,'(99F22.16)') ort_x(i),ABS(VX_EFF(i)),((ABS(VX(i,K,L)),K=1,MORB),L=1,MORB)
    ENDDO
    elseif (DIR.eq.'Y') then 
    DO I=1,NDY
      write(120,'(99F22.16)') ort_Y(i),ABS(VY_EFF(i)),((ABS(VY(i,K,L)),K=1,MORB),L=1,MORB)
    ENDDO
    endif
    close(120)
  end do

    IF(ALLOCATED(VX_EFF)) deallocate(VX_EFF)
    IF(ALLOCATED(VY_EFF)) deallocate(VY_EFF)
    IF(ALLOCATED(VX)) deallocate(VX)
    IF(ALLOCATED(VY)) deallocate(VY)
    IF(ALLOCATED(Vtrap_normX)) deallocate(Vtrap_normX)
    IF(ALLOCATED(Vtrap_normY)) deallocate(Vtrap_normY)
     
    return         
END  SUBROUTINE Get_VtrapProjection


SUBROUTINE Get_LZ(PSI,NOcc,NatVec,time,Rho_ij)

  USE SHARED_DIMS
  USE DVR_ALL
  USE W_INTERPARTICLE

  REAL*8,INTENT(IN) :: time 
  REAL*8 :: norm
  COMPLEX*16 :: Total_L_z,Orbital_L_z(MOrb,MOrb)
  COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb), INTENT(IN)       :: PSI
  COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb)                   :: PSIsave
  REAL*8, DIMENSION(MOrb),             INTENT(IN)       :: NOcc
  COMPLEX*16, DIMENSION(MOrb,MOrb),        INTENT(IN)       :: NatVec
  COMPLEX*16, DIMENSION(MOrb,MOrb),        INTENT(IN)       :: Rho_ij 
  COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb)       :: PSI2
  INTEGER :: K,Ind,J,I,l
 CALL GET_NOs(PSI,PSI2,NatVec)
! PSI2=PSI
 PSIsave=PSI2
 DO I=1,MOrb
!AIS 22JUL14 to cmpile with FFTW    CALL Get_Op_PSI_L_FFT_MKL(psi2(:,I),3)
    call normvxz(psi2(:,I),norm,NDX*NDY)
    do K=1,MORB  
      Orbital_L_z(I,K)=sum(dconjg(Psisave(:,K))*Psi2(:,I))
    end do
!    call xvixdzo(norm,psi(:,I),NDX*NDY)
 END DO
 write(*,*) 'ANGULAR MOMENTUM MATRIX ELEMENTS'
 write(*,'(3(2F16.8,1X))') (((Orbital_L_z(I,K)),I=1,Morb),K=1,MorB)

 Total_L_z=(0.d0,0.d0)
! DO K=1,MOrb
    DO I=1,MOrb
       Total_L_z=Total_L_z+(Nocc(I))*Orbital_L_z(I,I)
    END DO
! END DO
! Total_L_z=SUM(rho_ij*Orbital_L_z)

 write(6,*) 'The determined total angular momentum of this state at t=',time,'is <L_z>=', Total_L_z
 open(unit=122,file='L_z_expectation.dat', ACCESS = 'APPEND')
 WRITE(122,'(3F22.16)') time,Total_L_z
 write(122,'(A99)') 'ORBITAL ANGULAR MOMENTUM MATRIX ELEMENTS'
 DO I=1,MORB
    WRITE(122,'(99F22.16)') ((Orbital_L_z(I,K)),K=1,MORB)
 END DO
 close(122)
!       ind=1
!       do K=1,NDZ
!       do J=1,NDY
!       do I=1,NDX
!      write(112,2222) ort_X(I),"  ",ort_Y(J),"  ",ort_Z(K),"  "&
!     ,(DREAL(PSIsave(ind,l)),"  ",DIMAG(PSIsave(ind,l)),"  ",l=1,Morb)
!      write(113,2222) ort_X(I),"  ",ort_Y(J),"  ",ort_Z(K),"  "&
!     ,(DREAL(PSI2(ind,l)),"  ",DIMAG(PSI2(ind,l)),"  ",l=1,Morb)
!                 ind=ind+1
!       enddo
!       write(112,*)"    " 
!       write(113,*)"    " 
!       enddo
!       enddo

! 2222  format((138(F26.16,a3)))

END SUBROUTINE Get_LZ

SUBROUTINE GET_NOs(PSI,NOs,Natvec)
COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb), INTENT(IN)       :: PSI
COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb), INTENT(OUT)      :: NOs
COMPLEX*16, DIMENSION(MOrb,MOrb), INTENT(IN)              :: Natvec
COMPLEX*16                                                :: xx
INTEGER :: ind,K,J,I,k1,l1

   ind=1
   do K=1,NDZ
     do J=1,NDY
       do I=1,NDX
                 xx=Zero
        Do k1=1,Morb
            NOs(ind,k1)=Zero
          Do l1=1,Morb
!================Natural orbitals are
            NOs(ind,k1)=NOs(ind,k1)+Conjg(NatVec(l1,k1))*PSI(ind,l1) ! Natural orbitals 
! They are computed according to: \psi^{NO}_i= \sum_j U^{-1}_ij \phi_j = \sum_j U^*_ji \phi_j
! rho(x,x')= U U^-1 rho U U^-1=>(\phi_1^*,\phi_2^*,...\phi_M^*) U rho_d U^-1 ((\phi_1^*,\phi_2^*,...\phi_M^*)^T
          END DO
        END DO
        ind=ind+1
       enddo
     enddo
   enddo


END SUBROUTINE GET_NOs



END MODULE ANALYZER
