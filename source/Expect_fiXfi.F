       subroutine Expect_fiXfi(PSI) ! evaluates <Psi_no^*| x | Psi_no>
       USE   SHARED_DIMS
       USE   rR_hW
       USE   CI_All
       USE   W_INTERPARTICLE
       USE   DVR_ALL

c       implicit real*8 (a-h,o-z)
c       implicit integer (i-n)
       COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb)  :: PSI
       COMPLEX*16, DIMENSION(Morb)  :: v !,NO_Expectation_x(Morb)
       complex*16 xvec1,xvec2,rhoPR
       complex*16 xx,rij11,rij12,rij21,rij22,xN,yy
       real*8 time,t_total,ww
       integer :: Ndim,ind,k,j,i,k1,l1,JJ
       character*18 lname
c=================================================
        ndim=NDX*NDY*NDZ

       NO_Expectation_x=Zero
       
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
            v(k1)=v(k1)+Conjg(NatVec(l1,k1))*PSI(ind,l1)
                 END DO
                 END DO
        Do jj=1,Morb
        NO_Expectation_x(jj)=NO_Expectation_x(jj)
     .                 +ort_X(I)*v(jj)*dconjg(v(jj))
        EnDDO
       ind=ind+1
       enddo
       enddo
       enddo


       return
       end SUBROUTINE  Expect_fiXfi

