       subroutine U_Transform (time,VIN)
       USE   SHARED_DIMS
       USE   rR_hW
       USE   CI_All
       USE   W_INTERPARTICLE
       USE   DVR_ALL

c       implicit real*8 (a-h,o-z)
c       implicit integer (i-n)
       COMPLEX*16, DIMENSION(Nconf) :: VIN,VOUT,VIN1 !,VOUT1
       COMPLEX*16 scl
       COMPLEX*16 scl1,scl2
       real*8 time,t_total,ww,prefactor
       character*18 lname
       double precision dlgama,LogFKTR
       COMPLEX*16, DIMENSION(Morb,Morb) :: U ! NatVecInv Added for Unitary rotations NatVec is given
       INTEGER, DIMENSION(Morb) :: IPIV,Ivec,Nvec
       COMPLEX*16, DIMENSION(Morb,Morb) :: WORK
       real*8  xnrm,theta,DZNRM2
       integer :: NDIM,N2,K1,K2,N,NN1,N1,M,JJ,IND
       external DZNRM2,ZGETRF,ZGETRI

c=================================================
            VIN1=VIN
            VIN1=VIN1/DZNRM2(Nconf,VIN1,1)
c        write(6,*)"U-transformed CI Vec NORM :",xnrm
c===============   U psi_new=psi_old 
           U=NatVec
c          U=Transpose(NatVec)
c          U=Conjg(NatVec)
c           VIN1=Zero
c           VIN(51)=1.0d0
c           theta=PI/4
c           theta=0.0d0
c           U(1,1)= cos(theta)
c           U(1,2)= sin(theta)
c           U(2,1)=-sin(theta)
c           U(2,2)= cos(theta)
c         write(6,*) U(2,1),U(2,2)
c             Do i=1,Morb
c             Do j=1,Morb
c             U(i,j)=Conjg(U(i,j))
c             Enddo
c             Enddo

c            U(:,1)=U(:,1)/DZNRM2(2,U(:,1),1)
c            U(:,2)=U(:,2)/DZNRM2(2,U(:,2),1)


          
c         write(6,*) U(1,1),U(1,2)
c         write(6,*) U(2,1),U(2,2)
c===============    Get Inverse U - vectors obtained by diagonalizing rho_ij
c        xnrm=DZNRM2(Nconf,VIN,1)
c        write(6,*)"U-transformed CI Vec NORM :",xnrm
c        scl=U(1,1)*Conjg(U(1,1))+U(2,1)*Conjg(U(2,1))
c        write(6,*)"U-Norm1:",scl
c        scl=U(1,2)*Conjg(U(1,2))+U(2,2)*Conjg(U(2,2))
c        write(6,*)"U-Norm2:",scl
c====================== Inverse U ===========================
c        CALL ZGETRF( Morb, Morb , U, Morb, IPIV, INFO )  ! LU Factorization
c        iF(INFO.ne.0) write(6,*)"Problem in Invertion of U 1",INFO
c        iF(INFO.ne.0) pause
c        LWORK=Morb*Morb
c        CALL ZGETRI( Morb, U, Morb, IPIV, WORK, LWORK, INFO )
c        iF(INFO.ne.0) write(6,*)"Problem in Invertion of U 2",INFO
c        iF(INFO.ne.0) pause
c        U=Zero
c        U(1,1)=Zoner
c        U(2,2)=Zoner
c=================================================
c              goto 11
c          ii=0
c1        IF(ii.eq.1) U=Transpose(NatVec)
c         write(6,*) U(1,1),U(1,2)
c         write(6,*) U(2,1),U(2,2)
       VOUT=Zero
c       VOUT1=Zero
       ndim=Nconf
c       do ind=1,ndim
       do n2=0,Npar
c         CALL GET_Ind_vs_ii(i+1,Npar,2,Ivec,Nvec)
c         write(6,*)i+1,"I Indx vs conf",Nvec
c        IF(ABS(VIN(n2+1)).ge.0.5d-2) THEN
        Do k1=0,Npar-n2
        Do k2=0,n2
c       scl=    log(U(1,2))*(Npar-n2-k1)
c       scl=scl+log(U(1,1))*k1
c       scl=scl+log(U(2,2))*(n2-k2)
c       scl=scl+log(U(2,1))*k2

       prefactor = (dlgama( 1.0d0*(Npar-n2+1))/2
     &           + dlgama ( 1.0d0*(n2     +1))/2
     &           + dlgama ( 1.0d0*(k1+k2+1))/2
     &           + dlgama ( 1.0d0*(Npar-k1-k2+1))/2
     &           - dlgama ( 1.0d0*(k1     +1)) 
     &           - dlgama ( 1.0d0*(Npar-n2-k1+1))
     &           - dlgama ( 1.0d0*(k2     +1)) 
     &           - dlgama ( 1.0d0*(n2-k2  +1)) )
c       prefactor = (LogFKTR( (Npar-n2))/2
c     &           +  LogFKTR( (n2     ))/2
c     &           +  LogFKTR ( (k1+k2))/2
c     &           +  LogFKTR ( (Npar-k1-k2))/2
c     &           -  LogFKTR ( (k1     )) 
c     &           -  LogFKTR ( (Npar-n2-k1))
c     &           -  LogFKTR ( (k2     )) 
c     &           -  LogFKTR ( (n2-k2  )) )


c       write(6,*)k1+k2,scl
       scl=    U(1,2)**((Npar-n2-k1))
       scl=scl*U(1,1)**(k1)
       scl=scl*U(2,2)**((n2-k2))
       scl=scl*U(2,1)**(k2)
c       scl=scl*prefactor
c       IF(k1+k2.le.51) write(6,*)exp(scl),exp(prefactor)
C       IF(k1+k2.le.51) write(6,*)exp(DREAL(scl+prefactor)),
C     $ exp(DIMAG(scl+prefactor))
c        scl1=exp(scl+prefactor)
        scl1=scl*exp(prefactor)
c       IF(k1+k2.le.51) write(6,*)exp(scl),exp(prefactor)
C     $ exp(DIMAG(scl+prefactor))

c       IF(exp(prefactor).ge.1.0d+15) scl1=Zero
c       IF(ABS(scl1).le.0.2d-15) scl1=Zero
c       IF(k1+k2.eq.51) write(6,*)scl1*VIN1(n2+1)

c        IF(ABS(scl1).le.0.2d-15) then
        VOUT(Npar-(k1+k2)+1)=VOUT(Npar-(k1+k2)+1)+scl1*VIN1(n2+1)
c        else
c        VOUT1(Npar-(k1+k2)+1)=VOUT1(Npar-(k1+k2)+1)+scl1*VIN1(n2+1)
c        endif

c       IF(k1+k2.eq.51) write(6,*)VOUT(Npar-(k1+k2)+1)
        Enddo
        Enddo
c       ENDIF
       enddo
c         VOUT=VOUT+VOUT1
c================Ofirs Version
              goto 12
11       continue
       N=Npar
       do NN1=0,N
       VOUT(N-NN1+1)=Zero
       do n1=0,N
       scl1=Zero
c         CALL GET_Ind_vs_ii(n1+1,N,2,Ivec,Nvec)
c         write(6,*)n1+1,"I Indx vs conf",Nvec
c        write(6,*)"From",MAX(NN1+n1-N,0),MIN(NN1,n1)
       do m=MAX(NN1+n1-N,0),MIN(NN1,n1)
       scl=    U(1,1)**(1.0d0*m)
       scl=scl*U(2,1)**(1.0d0*(NN1-m))
       scl=scl*U(1,2)**(1.0d0*(n1-m))
       scl=scl*U(2,2)**(1.0d0*(N-n1-NN1+m))
c       prefactor = Exp(dlgama ( 1.0d0*(n1+1))/2
c     &           + dlgama ( 1.0d0*(N-n1+1))/2
c     &           + dlgama ( 1.0d0*(NN1+1))/2
c     &           + dlgama ( 1.0d0*(N-NN1+1))/2
c     &           - dlgama ( 1.0d0*(m     +1)) 
c     &           - dlgama ( 1.0d0*(NN1-m+1))
c     &           - dlgama ( 1.0d0*(n1-m     +1)) 
c     &           - dlgama ( 1.0d0*(N-n1-NN1+m +1)) )

       scl2=U(2,1)**(1.0d0*NN1)*U(1,2)**(1.0d0*n1)
       scl2=scl2*U(2,2)**(1.0d0*(N-n1-NN1))
       scl2=scl2*Exp(dlgama ( 1.0d0*(n1+1))/2+dlgama(1.0d0*(N-n1+1))/2)

       prefactor = Exp(dlgama ( 1.0d0*(n1+1))/2
     &           + dlgama ( 1.0d0*(N-n1+1))/2
     &           + dlgama ( 1.0d0*(NN1+1))/2
     &           + dlgama(1.0d0*(N-NN1+1))/2
     &           - dlgama ( 1.0d0*(m     +1)) 
     &           - dlgama ( 1.0d0*(NN1-m+1))
     &           - dlgama ( 1.0d0*(n1-m     +1)) 
     &           - dlgama ( 1.0d0*(N-n1-NN1+m +1)) )



       scl=scl*prefactor
       scl1=scl1+scl
c       IF(n1.eq.51) write(6,*) scl
       Enddo
c       write(6,*) U(1,1)*U(2,2)/U(2,1)/U(1,2) 
       VOUT(N-NN1+1)=VOUT(N-NN1+1)+scl1*VIN1(N-n1+1)
       Enddo
       Enddo
12       continue



        xnrm=DZNRM2(Nconf,VOUT,1)
        write(6,*)"U-transformed CI Vec NORM :",xnrm
c           IF(ABS(1.0-xnrm).ge.0.000001d0) return
           IF(ABS(1.0-xnrm).ge.0.0001d0) return
c         VIN1=VOUT/xnrm
c          ii=ii+1
c         IF(ii.eq.1) goto 1
c=================================================
c========= Write_OUT results CI and ORBs
c=================================================
       t_total=time
       if(t_total.lt.99999.99999999D0) then
       Write( lname, '(F10.4,a8)') t_total,"coef.TRN"
       endif
       if(t_total.lt.9999.99999999D0) then
       Write( lname, '(F10.5,a8)') t_total,"coef.TRN"
       endif
       if(t_total.lt.999.99999999D0) then
       Write( lname, '(F10.6,a8)') t_total,"coef.TRN"
       endif
       if(t_total.lt. 99.99999999D0) then
       Write( lname, '(F10.7,a8)') t_total,"coef.TRN"
       endif
       if(t_total.lt.  9.99999999D0) then
       Write( lname, '(F10.8,a8)') t_total,"coef.TRN"
       endif

      open(unit=111,file=lname,form='formatted')
      write(111,*)"          "
      write(111,'(a1,a28,138(I26,a3))')"#","",(jj,"  ",jj=2,10)
       
       ndim=Nconf
       do ind=1,ndim
c        IF(ABS(VOUT(ind)).ge.0.5d-14) THEN
      write(111,2222) ind,"  "
     $,DREAL(VOUT(ind)),"   ",DIMAG(VOUT(ind)),"   "
     $,t_total
c        ENDIF
       enddo
       close(111)
       return
2222  format((i26,a3,10(F26.16,a3)))
       end SUBROUTINE U_Transform

