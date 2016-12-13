! vim:fdm=marker:

MODULE CORRELATIONFUNCTIONS
use Prop_MB
USE W_INTERPARTICLE

IMPLICIT NONE
CONTAINS

subroutine getDoubleAsString(x,doubleString)
!       {{{ writes a double in ten digit format, 
!           usefull  for naming files without blanks       
IMPLICIT NONE
integer,parameter         :: DBL=8
real(kind=DBL)            :: x
CHARACTER(len=10)         :: doubleString

if (x .LE. 999999.999d0) WRITE(doubleString,'(F10.3)') x
if (x .LE. 99999.9999d0) WRITE(doubleString,'(F10.4)') x
if (x .LE. 9999.99999d0) WRITE(doubleString,'(F10.5)') x
if (x .LE. 999.999999d0) WRITE(doubleString,'(F10.6)') x
if (x .LE. 99.9999999d0) WRITE(doubleString,'(F10.7)') x
if (x .LE. 9.99999999d0) WRITE(doubleString,'(F10.8)') x

IF ((x .LT. 0.d0) .OR. (x .GE. 1000000.0d0)) THEN
  WRITE(6,*)"get double as string: x not implemented"
  STOP
END IF

!}}}
end subroutine getDoubleAsString

subroutine get_maximum(PSI,rho_jk,maximum)
!---{{{
use DVR_ALL

  implicit none
!ORG  integer,parameter  :: NOGRIDPOINTS=NDX*NDY*NDZ
  integer  :: NOGRIDPOINTS
!  complex*16,DIMENSION(NOGRIDPOINTS,Morb),intent(in) :: PSI 
    COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb)  :: PSI
  complex*16,DIMENSION(Morb,Morb),intent(in) :: rho_jk
  complex*16,intent(out)  :: maximum
  integer                 :: m,Jorb,Korb
  complex*16              :: G1

  NOGRIDPOINTS=NDX*NDY*NDZ

  maximum = (0.d0,0.d0)
  do m=1,NOGRIDPOINTS
 
 !compute G1
     G1 = (0.d0,0.d0)
     do Korb = 1,Morb  
        do Jorb=1,Morb

           G1 = G1 + rho_jk(Jorb,Korb) * DCONJG(PSI(m,Jorb))*PSI(m,Korb)

        end do
     end do

     if (REAL(CDABS(G1)) .gt. REAL(CDABS(maximum))) maximum = G1
 
  end do
!---}}}
end subroutine get_maximum


subroutine get_correlations(time,PSI,mode,dilation,rho_jk,rho_ijkl,pathPRN) 
!--------------------------------------------------------------------------------------------
! computes the correlation functions 
! G1(r,r), G1(r,r'), G1(r',r') and G2(r,r',r',r), where r=(x,y,z), 
! and writes them to a file. The first six columns contain
! x,y,z,x',y',z'. Currently only 1D has been tested.
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
! pathPRN       :  character*10 path,pathCI,pathORB,pathORBK 
!              ='DATA/g1_RR'
!              ='DATA/g1_KK'
!--------------------------------------------------------------------------------------------
!  {{{
use CI_ALL,ONLY:Npar
use DVR_ALL
  implicit none
      character*10 pathPRN
!  integer,parameter  :: NOGRIDPOINTS=NDX*NDY*NDZ
  integer  :: NOGRIDPOINTS
  integer,intent(in) :: mode,dilation
  real*8,intent(in)  :: time
  real*8     :: w,w1,dk,mom_X(NDX),mom_Y(NDY),mom_Z(NDZ)
  real*8     :: x1,y1,z1,x2,y2,z2
!  complex*16,DIMENSION(NOGRIDPOINTS,Morb),intent(in) :: PSI 
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

  complex*16         :: maximum,integ,integ2
  real*8,parameter   :: TOL=0.00001d0 
  real*8     :: threshold
  real*8     :: outG1m_m,outReG1m_n,outImG1m_n,outG1n_n,outG2m_n
  real*8             :: xi,xf,yi,yf,zi,zf

  integ=0.d0
  integ2=0.d0


  NOGRIDPOINTS=NDX*NDY*NDZ
   xi=Time_xint
   xf=Time_xfnl

  if(mode.eq.1) then
     w=DSQRT((xf-xi)/NDX)
  else if (mode.eq.2) then
     call get_mom(xi,xf,NDX,mom_X)
     dk     = mom_X(2)-mom_X(1) 
     w      = DSQRT( dk )
     mom_Y=0d0
     mom_Z=0d0
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
  else if (mode.eq.2) then 
     aString = 'k-'
  else 
     write(*,*)"ERROR in get_correlations"
     stop
  endif   

  if ((NPar.gt.0) .and. (NPar.lt.10)) then 
    write (NParAsString1, '(I1)') NPar
    fullName=timeasString//'N'//trim(NParAsString1)//'M'//trim(MorbAsString)//aString//'correlations.dat'
  else  if ((NPar.ge.10) .and. (NPar.lt.100)) then 
    write (NParAsString2, '(I2)') NPar
    fullName=timeasString//'N'//trim(NParAsString2)//'M'//trim(MorbAsString)//aString//'correlations.dat'
  else if ((NPar.ge.100) .and. (NPar.lt.1000)) then 
    write (NParAsString3, '(I3)') NPar
    fullName=timeasString//'N'//trim(NParAsString3)//'M'//trim(MorbAsString)//aString//'correlations.dat'
  else if ((NPar.ge.1000) .and. (NPar.lt.10000)) then 
    write (NParAsString4, '(I4)') NPar
    fullName=timeasString//'N'//trim(NParAsString4)//'M'//trim(MorbAsString)//aString//'correlations.dat'
  else  if ((NPar.ge.10000) .and. (NPar.lt.100000)) then 
    write (NParAsString5, '(I5)') NPar
    fullName=timeasString//'N'//trim(NParAsString5)//'M'//trim(MorbAsString)//aString//'correlations.dat'
  else if ((NPar.ge.100000) .and. (NPar.lt.1000000)) then 
    write (NParAsString6, '(I6)') NPar
    fullName=timeasString//'N'//trim(NParAsString6)//'M'//trim(MorbAsString)//aString//'correlations.dat'
  else if ((NPar.ge.1000000) .and. (NPar.lt.10000000)) then 
    write (NParAsString7, '(I7)') NPar
    fullName=timeasString//'N'//trim(NParAsString7)//'M'//trim(MorbAsString)//aString//'correlations.dat'
  else if ((NPar.ge.10000000) .and. (NPar.lt.100000000)) then 
    write (NParAsString8, '(I8)') NPar
    fullName=timeasString//'N'//trim(NParAsString8)//'M'//trim(MorbAsString)//aString//'correlations.dat'
  else 
    NParAsString4='XXXX'
    fullName=timeasString//'N'//trim(NParAsString4)//'M'//trim(MorbAsString)//aString//'correlations.dat'
  end if
 ! needs to be edited for other dvrs with nonconstant weights! 
  call get_maximum(PSI,rho_jk,maximum)
  threshold = Kdip * DBLE(maximum)/w**2 ! set the minimal density to be plotted
  
    fullName=timeasString//aString//'correlations.dat'
  open(unit=12,file=trim(pathPRN)//'/'//trim(fullName),form='formatted') !STR 2014 corrected for DIR path
  do n=1,NOGRIDPOINTS

        w=weight(n)
!compute G1n_n
     G1n_n = (0.d0,0.d0)
     do Jorb = 1,Morb  
        do Korb=1,Morb

           G1n_n = G1n_n + rho_jk(Jorb,Korb) * DCONJG(PSI(n,Jorb))*PSI(n,Korb)

        end do
     end do
     G1n_n=G1n_n*1.d0/w**2/Npar
     outG1n_n=DBLE(G1n_n)
!---------------
     do m=1,NOGRIDPOINTS         

        w1=weight(m)
!compute G1m_m
        G1m_m = (0.d0,0.d0) 
        do Jorb = 1,Morb  
           do Korb=1,Morb

              G1m_m = G1m_m + rho_jk(Jorb,Korb) * DCONJG(PSI(m,Jorb))*PSI(m,Korb)

           end do
        end do
        G1m_m=G1m_m*1.d0/w1**2/Npar
        outG1m_m=DBLE(G1m_m)   
!----------------
        
!compute G1m_n
        G1m_n = (0.d0,0.d0) 
        do Jorb = 1,Morb  
           do Korb=1,Morb

              G1m_n = G1m_n + rho_jk(Jorb,Korb) * DCONJG(PSI(m,Jorb))*PSI(n,Korb)

           end do
        end do
        G1m_n=G1m_n*1.d0/(w*w1)/Npar 
        outReG1m_n=DBLE(G1m_n)   
        outImG1m_n=DIMAG(G1m_n)
        if (m.eq.n) integ=integ+G1m_n*w*w1
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
        G2m_n=G2m_n*1.d0/(w*w*w1*w1)/Npar/(Npar-1)
        integ2=integ2+G2m_n*w*w*w1*w1
!        write(*,*) G2m_n,integ2
        outG2m_n=DBLE(G2m_n)
!----------------
!        IF((DABS(outG1m_m).LT. MAX(TOL,threshold))     .OR.   &
!           (DABS(outG1n_n).LT. MAX(TOL,threshold)) )   THEN 
!           outReG1m_n =  -600.d0    
!           outImG1m_n =  -600.d0    
!           outG2m_n   = -1.d0    
!        END IF

        call get_ijk_from_m(m,NDX,NDY,i_m,j_m,k_m)
        call get_ijk_from_m(n,NDX,NDY,i_n,j_n,k_n)
        if (mode.eq.1) then
           write(12,6666) ort_X(i_m)," ",ort_Y(j_m)," ",ort_Z(k_m)," ",ort_X(i_n)," ",ort_Y(j_n)," ",ort_Z(k_n)," ", &
                       outG1m_m," ",outReG1m_n," ",outImG1m_n," ",outG1n_n," ",outG2m_n
        else if (mode.eq.2) then

           write(12,6666) mom_X(i_m)," ",mom_Y(j_m)," ",mom_Z(k_m)," ",mom_X(i_n)," ",mom_Y(j_n)," ",mom_Z(k_n)," ", &
                       outG1m_m," ",outReG1m_n," ",outImG1m_n," ",outG1n_n," ",outG2m_n
        end if

        if(MOD(m,NDX).eq.0) write(12,*)"                                            " 

     end do
                             
  end do         
  close(12)
 
 
6666  format(6(F8.3,a1),6(F16.9,a1))
!---}}}
!     write(*,*) 'int of rho - old', integ, 'mode', mode
!     write(*,*) 'int of rho2 - old', integ2, 'mode', mode
  write(6,'(a28,F12.6,a4,i2,a7,F8.4)') 'g1: Intg DNS1(R[mode=1]|K)=',ABS(integ),' mode ',mode,' T=',time
  write(6,'(a28,F12.6,a4,i2,a7,F8.4)') 'g1: Intg DNS2(R[mode=1]|K)=',ABS(integ2),' mode ',mode,' T=',time
!     write(*,*) 'norm of orb',sum(conjg(PSI(:,1))*PSI(:,1)),'sum of 2b-ele',sum(rho_ijkl)
end subroutine get_correlations

subroutine get_ijk_from_m(m,gdimx,gdimy,i,j,k)
!the routine returns i,j,k given m---{{{
!of a vector that stores the ijk-th element in
!the position m=i+gdimx(j-1)+gdimy*gdimx*(k-1)
![(ijk) start with (111)]
implicit none
integer,intent(IN)      :: m,gdimx,gdimy
integer,intent(OUT)     :: i,j,k
integer                 :: n

  n = m
  k = (n-1)/(gdimx*gdimy) + 1
  n =  n-(k-1)*gdimx*gdimy
  j = (n-1)/gdimx + 1
  n =  n-(j-1)*gdimx
  i =  n

!---}}}
end subroutine get_ijk_from_m

subroutine get_FFTmom(xi,xf,gdim,mom)
USE SHARED_DIMS, ONLY: PI
 !---{{{ returns momentum grid in fft order
! xi   - left border of box
! xf   - right border of box ( =last grid point + dx)
! gdim - number of grid points
! mom  - array of momentum space grid points in FFT order
implicit none

integer,intent(IN)              :: gdim
!double precision,intent(IN)     :: xi,xf
double precision,intent(OUT)    :: mom(gdim)

double precision                :: dp
integer                         :: i
  real*8             :: xi,xf

   xi=Time_xint
   xf=Time_xfnl
mom(1)=0.d0

if (gdim.gt.1) then 
   dp             =     2.d0*PI/(xf-xi)
end if

do i=2,gdim/2
  mom(i) = (i-1)*dp
enddo

do i=gdim/2+1,gdim !put momenta in FFT order (0,dp,...,N/2*dp,(-N/2+1)*dp,...,-dp)
  mom(i)=1.d0*(i-gdim)*dp
end do

!---}}}
end subroutine get_FFTmom

subroutine get_mom(xl,xu,gdim,mom)
USE SHARED_DIMS, ONLY: PI
!---{{{ returns momentum grid in ascending order
implicit none

integer,intent(IN)              :: gdim
double precision,intent(IN)     :: xl,xu
double precision,intent(OUT)    :: mom(gdim)

double precision                :: dp
integer                         :: i

if (gdim.gt.1) then 
   dp             =     2.d0*PI/(xu-xl)

   if (MOD(gdim,2).eq.0) then !gdim even

      do i=1,gdim
         mom(i)=(i-gdim/2-1)*dp  ! e.g.: (-dp,0,dp,2dp)
      end do

   else ! gdim uneven

      do i=1,gdim
         mom(i)=(i-gdim/2-1)*dp ! e.g.: (-2dp,-dp,0,dp,2dp)
      end do

   end if 

else if (gdim.eq.1) then

   mom(1) = 0.d0

else

   write(*,*)"ERROR in dimension"
   stop

end if

!---}}}
end subroutine get_mom

subroutine get_dilation(dilation,Psi,rho_jk)
use DVR_ALL
!       {{{ computes the optimal k-space intervall for the FT
!       if unwanted, set dilation = 1 everywhere
  implicit none  
  integer,intent(inout)     :: dilation
!  complex*16,DIMENSION(NDX*NDY*NDZ,Morb),intent(in) :: PSI 
  COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb)  :: PSI
  complex*16,DIMENSION(Morb,Morb),intent(in) :: rho_jk

  integer                   :: m,iLower,iUpper
  integer                   :: Jorb,Korb
  complex*16                :: cmax,G1m_m
  real*8                    :: G1max,THRESHOLD,G1(NDX*NDY*NDZ)
  real*8                    :: imin,imax

!  if((NDY.gt.1) .or. (NDZ.gt.1)) then
!     write(*,*)"only 1d implemented"
!     stop
!  end if 

  THRESHOLD =       0.01d0       
  dilation  =       1
  G1max     =       0.d0

  call get_maximum(PSI,rho_jk,cmax)
  G1max = DBLE(cmax)

!---compute G1m_m

     do m=1,NDX*NDY*NDZ 

        G1m_m = (0.d0,0.d0) 
        do Jorb = 1,Morb  
           do Korb=1,Morb

              G1m_m = G1m_m + rho_jk(Jorb,Korb) * DCONJG(PSI(m,Jorb))*PSI(m,Korb)

           end do
        end do

        G1(m)=DBLE(G1m_m)   

     end do 

!----------------
!----find boundaries

      do m=1,NDX/2
        if(G1(m)/G1max .GT. THRESHOLD) THEN
          iLower = m  
          iUpper = NDX-m+1
          exit
        end if
        if(G1(NDX-m+1)/G1max .GT. THRESHOLD) THEN
          iLower = m  
          iUpper = NDX-m+1
          exit
        end if  
      end do 

        imax     = 1.d0*(iUpper-(NDX)/2 - 1)
        imin     = 1.d0*(iLower-(NDX)/2 - 1)
        imax     = MAX(DABS(imax),DABS(imin))    
        dilation = MAX(INT(NDX/2.d0 /(2.d0*imax)),1)

!!}}}      
end  subroutine get_dilation

subroutine pedestrian_FT(PSI,FTPSI,dilation) 
!---{{{
use DVR_ALL
use W_INTERPARTICLE
  implicit none
  integer                  :: Iorb,m,n

  complex*16,DIMENSION(NDX*NDY*NDZ,Morb),intent(in)  :: PSI 
  complex*16,DIMENSION(NDX*NDY*NDZ,Morb),intent(out) :: FTPSI
 
  real*8,dimension(NDX*NDY*NDZ) :: kvecs
  real*8                        :: xall,xfirst,xlast,dk,dx,kmin,k,norm,xi,xf
  integer,intent(in)                    :: dilation
  integer :: DVRMETHOD
!---{{{ check for consistency
      DVRMETHOD=Time_DVRMETHODX

  if ((NDY.gt.1) .or. (NDZ.gt.1)) then
     write(*,*) "more than 1-D has not been implemented yet"
     return
  end if

  if ((DVRMETHOD .ne. 5).and.(DVRMETHOD .ne. 4))  then
     write(*,*) "other than exp-dvr or FFT has not been implemented yet"
     return
  end if
!---}}}

   xi=Time_xint
   xf=Time_xfnl
  xfirst = xi
  xall   = xf-xi
  dx     = xall/NDX 
  dk     = 2.d0*PI/(xall*dilation) 
  kmin   = -(NDX/2) * dk
  
  FTPSI=(0.d0,0.d0)

  do Iorb=1,Morb


     do m=1,NDX   

        k = kmin+(m-1)*dk

        do n=1,NDX
           FTPSI(m,Iorb) = FTPSI(m,Iorb) + dx * PSI(n,Iorb)/DSQRT(dx*2.d0*PI) * CDEXP( (0.d0,-1.d0)*k*ort_X(n) )
        end do 

     end do 

     call normvxz(FTPSI(:,Iorb),norm,NDX) ! STR 2014 to activate normalization in K space
     call xvixdzo(norm,FTPSI(:,Iorb),NDX) ! STR 2014

  end do


!---}}}
end subroutine pedestrian_FT


END MODULE CORRELATIONFUNCTIONS
