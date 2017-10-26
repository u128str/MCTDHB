       subroutine print_rho(time,PSI,path,Rho_PR)
       USE   SHARED_DIMS
       USE   rR_hW
       USE   CI_All
       USE   W_INTERPARTICLE
       USE   DVR_ALL

c       implicit real*8 (a-h,o-z)
c       implicit integer (i-n)
       COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb)  :: PSI
       COMPLEX*16, DIMENSION(NDX*NDY*NDZ)  :: Rho_PR
       COMPLEX*16, DIMENSION(Morb)  :: v !,NO_Expectation_x(Morb)
       complex*16 xvec1,xvec2,rhoPR
       complex*16 xx,rij11,rij12,rij21,rij22,xN,yy
       real*8 time,t_total,ww
       integer :: jj,Ndim,ind,k,j,i,k1,l1
       character*18 lname
       character*10 path
c=================================================
        ndim=NDX*NDY*NDZ
       t_total=time
       if(t_total.lt.99999.99999999D0) then
       Write( lname, '(F10.4,a8)') t_total,"time.dat"
       endif
       if(t_total.lt.9999.99999999D0) then
       Write( lname, '(F10.5,a8)') t_total,"time.dat"
       endif
       if(t_total.lt.999.99999999D0) then
       Write( lname, '(F10.6,a8)') t_total,"time.dat"
       endif
       if(t_total.lt. 99.99999999D0) then
       Write( lname, '(F10.7,a8)') t_total,"time.dat"
       endif
       if(t_total.lt.  9.99999999D0) then
       Write( lname, '(F10.8,a8)') t_total,"time.dat"
       endif

c       NO_Expectation_x=Zero
      open(unit=112,file=trim(path)//'/'//lname,form='formatted')
      write(112,*)"          "
      write(112,'(a1,a28,999(I26,a3))')"#","",(jj,"  ",jj=2,10+5*Morb)
       
       xN=1.0d0*Npar

       ind=1
       do K=1,NDZ
       do J=1,NDY
       do I=1,NDX
                 xx=Zero
                 yy=Zero
                 ww=weight(ind)
                 Do k1=1,Morb
            v(k1)=Zero
                 Do l1=1,Morb
            xx=xx+dconjg(PSI(ind,l1))*AllZRIJ(l1,k1)*PSI(ind,k1)
!================Natural orbitals are
            v(k1)=v(k1)+Conjg(NatVec(l1,k1))*PSI(ind,l1) ! Natural orbitals are ORG
! They are computed according to: \psi^{NO}_i= \sum_j U^{-1}_ij \phi_j = \sum_j U^*_ji \phi_j
! rho(x,x')= U U^-1 rho U U^-1=>(\phi_1^*,\phi_2^*,...\phi_M^*) U rho_d U^-1 ((\phi_1^*,\phi_2^*,...\phi_M^*)^T
                 END DO
            yy=yy+Nocc(k1)*dconjg(v(k1))*v(k1)
                 END DO
                 rhoPR=xx/ww/ww
                 Rho_PR(ind)=rhoPR/xN !STR 1Aug2014
      write(112,2222) ort_X(I),"  ",ort_Y(J),"  ",ort_Z(K),"  ",ww,"  "
     $,DREAL(VTRAP_EXT(ind)),"  "
     $,DREAL(rhoPR/xN),"   ",DIMAG(rhoPR/xN),"   "
     $,DREAL(yy/ww/ww/xN),"   ",DIMAG(yy/ww/ww/xN),"   "
     $,(DREAL(PSI(ind,jj)),"  ",DIMAG(PSI(ind,jj)),"  ",jj=1,Morb)
     $,(DREAL(v(jj)),"  ",DIMAG(v(jj)),"   ",jj=1,Morb)
     $,(Nocc(jj)/Npar,"   ",jj=1,Morb)
     $,t_total
                 ind=ind+1
       enddo
       write(112,*)"    " 
       enddo
       enddo


       close(112)
       return
2222  format((999(F26.16,a3)))
       end SUBROUTINE print_rho

