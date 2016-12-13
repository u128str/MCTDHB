c============= Fit of the interparticle potential - to make it separable
c============ F(x-x')=SUM F_i f^*_i(x)f_i(x')
c============     Wxx - VINTERP(ort(i),ort(j))
       subroutine get_POTFIT(ndim,ort,w,F,Fij) 
c       implicit real*8 (a-h,o-z)
c       implicit integer (i-n)
       integer NDIM,LWORK
       real*8  ort(ndim),RFij(ndim,ndim),w(ndim),RF(ndim)
       real*8  C(ndim,ndim),VINTERP
c       complex*16 Wxx(ndim,ndim)
       real*8  WORK(ndim*ndim+2*ndim)
       real*8 Fij(ndim,ndim),F(ndim)
       integer I,J,K,L,info
       real*8 XINTEGR,ALPHA,BETA,xdelta,xint,xnrmf,xtol
       EXTERNAL DSYEV

c================== Construction of Fij on the DVR set
       xintegr=0.0d0
       Do i=1,ndim
       Do j=1,ndim
       C(i,j) =VINTERP(ort(i),ort(j))
c       Wxx(i,j)=VINTERP(ort(i),ort(j))*(1.0d0,0d0)
       Enddo
       Enddo
       ALPHA=1.0d0
       BETA=0.0d0
c================evaluation of the potential density
       CALL DSYMM('L','L',ndim,ndim,ALPHA,C,ndim,C,ndim,BETA,RFij,ndim)

       LWORK=ndim*ndim+2*ndim
c      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
c================= deagonalization potential density to obtain natural potentials and their occupations
       CALL DSYEV('V','L',ndim,RFij,ndim,RF,WORK,LWORK,INFO)

c       write(6,*)"Diagonalization of interparticle interaction is ",INFo
c       write(6,*)"Eigenvalues are ",F
c================ Construction of the L^2-error Delta...
       xdelta=0.0d0
       xint=0.0d0
       xnrmf=0.0d0
       Do k=1,ndim
       Do l=1,ndim
        xint=0.0d0
        Do j=1,ndim
        xint=xint+DSQRT(DABS(RF(j)))*RFij(k,j)*RFij(l,j)
        Enddo
        xdelta=xdelta+(C(k,l)-xint)*(C(k,l)-xint)
c       write(6,*)"(V(x-x'))x=",ort(i)," x'=",ort(i),"= ",C(i+1,i)-xint
       xnrmf=xnrmf+C(k,l)*(w(k)*w(l))**2
       Enddo
       Enddo
       write(6,*)" Int (V(x-x')=",xnrmf
        xtol=1.0d-8
       write(6,*)" Delat SUM_ij (V(x-x')-Vappx)=",xdelta
       IF(xdelta.ge.xtol) then
       write(6,*)" Delta SUM_ij (V(x-x')-Vappx) >",xtol
       write(6,*)"Arroximation to interparticle interaction is BAD"
       pause
       ENDIF      
c        stop
       Do k=1,ndim
       Do l=1,ndim
       Fij(k,l)=RFij(k,l)
       Enddo
       F(k)=RF(k)
       Enddo
c       write(6,*)" F^0.5 out:"
c       write(6,*)(DSQRT(ABS(DREAL(RF(i)))),i=1,ndim)
c       write(6,*)(F(i),i=1,ndim)
c       write(6,*)(Fij(i,i),i=1,ndim)
c       pause
       return
       end subroutine get_POTFIT
c===================== Analitical expression for Interparticle interaction
       FUNCTION VINTERP(x,y)
       implicit real*8 (a-h,o-z)
       complex*16 compt
       real*8 PI,VINTERP,x,y,D
       parameter        (PI = 3.141592653589793238462643d0)
       compt=exp(-(x-y))
       compt=x**2+0.0d0*y**2
       D=0.05d0
c       D=0.5d0
       compt=Exp(-(x-y)**2/(2*D*D))/DSQRT(2*PI)/D
       compt=(Exp(-200*(x-y)**2))*DSQRT(200/PI)
c       w=1.0d0
c       w=3.25d0
c       w=10.0d0
c       w=3.2d0
c       w=100.0d0
c       compt=Exp(-w*(x-y)**2)
c       compt=Exp(-0.5*(x-y)**2)*(x-y)**2
c       compt=Exp(-0.5*(x-y)**2)-0.5*Exp(-2*(x-y)**2)
c       compt=1.0d0/DSQRT((x-y)**2+1.0d0)
c       e=1.0d0 
c       compt=e/(e**2+(x-y)**2)/Pi
       VINTERP=DREAL(compt)
       end 
