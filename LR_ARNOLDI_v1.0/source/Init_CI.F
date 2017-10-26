C=============== INtegral Hij Wijkl Evaluation 
       subroutine Init_CI
       USE SHARED_DIMS
       USE rR_hW
       USE CI_All
       implicit NONE
c=================== MPI 
c       INCLUDE 'mpif.h'
c       INTEGER ::  ierr,MYID,numprocs
c==========================================================
       INTEGER ::  I,J,L,cK,cL,CI,cJ,i2,P,K
c==========================================================
       REAL*8 :: CNK
c       EXTERNAL   CNK
c==========================================================
c==========================================================
            I=Npar+Morb-1
c            J=Npar
            J=Morb
       Nconf=NINT(CNK(I,Npar))
      write(6,*)" You task for ",Npar," particle and ", Morb,
     $" orbitals is  C(N+M-1,N) long:", Nconf,CNK(I,Npar)
c========================== MCNK(i,j) and TermIndex  ARE SENT directelly 
c========= Reconstruction of the Binomial Matrix MCNK ===================================
            Do i=0,Morb+Npar-1
            Do j=0,Morb-1
            MCNK(i,j)=NINT(CNK(i,j))
            EndDo
c            write(6,*)("C(",I,",",K,")=",MCNK(i,k),K=0,Morb)
            EndDo

c             pause
c==========================================================
c========== Reconstruction of the Term_Index_1B(L)=0+100*cJ+10000*cK+1000000*0
                L=0
               Do cK= 1,Morb
               Do cJ=cK,Morb
               L=L+1
c               P=cJ+100*cK
               P=cK+100*cJ
               TERM_INDEX_1B(L)=P
               TERM_REQ_1B(L)=L
c=========================== RESCALE(I) array
                     RESCALE_1B(L)= 1.0d0*ZONER+1.0d0*ZONEI
       IF(cJ.eq.cK)  RESCALE_1B(L)= 1.0d0*ZONER+0.0d0*ZONEI
c         write(6,*)P,"cJ-cK:",cJ,cK, RESCALE_1B(L)
C         write(6,*)P,"cJ-cK:",cK,cJ
c         write(6,*)P,"cJ-cK:",cJ,cK
c=========================== End Rescale
                EndDo
                EndDo
                MaxTrm1b=Morb*(Morb+1)/2
c========== Construction of the TERM_INDEX_2B(L)=cI+100*cJ+10000*cK+1000000*cL
                L=0
               Do cK= 1,Morb
               Do cL=cK,Morb
                  i2=cL
                  Do cI=cK,Morb
                  Do cJ=i2,Morb
          L=L+1
c          P=cI+100*cJ+10000*cK+1000000*cL
          P=cK+100*cL+10000*cI+1000000*cJ
          TERM_INDEX_2B(L)=P
          TERM_REQ_2B(L)=L
c=========================== RESCALE(I) array
C                                      RESCALE_2B(L)=4*ZONER+4*ZONEI
C        IF(cI.eq.cJ)                  RESCALE_2B(L)=2*ZONER+2*ZONEI
C        IF(cK.eq.cL)                  RESCALE_2B(L)=2*ZONER+2*ZONEI
C        IF((cI.eq.cK).and.(cJ.eq.cL)) RESCALE_2B(L)=2*ZONER+2*ZONEI
C        IF((cJ.eq.cI).and.(cL.eq.cK)) THEN
C                                      RESCALE_2B(L)= 1*ZONER+1*ZONEI
C        IF(cJ.eq.cK)                  RESCALE_2B(L)= 0.5*ZONER+0.5*ZONEI
C        END IF
                                      RESCALE_2B(L)=0.5*ZONER+0.5*ZONEI
       IF((cI.eq.cK).and.(cJ.eq.cL)) RESCALE_2B(L)=0.5*ZONER+0.0*ZONEI
c       IF(cI.eq.cJ)                  RESCALE_2B(L)=0.5*ZONER+0.5*ZONEI
c       IF(cK.eq.cL)                  RESCALE_2B(L)=0.5*ZONER+0.5*ZONEI
c       IF((cI.eq.cK).and.(cJ.eq.cL)) RESCALE_2B(L)=0.5*ZONER+0.5*ZONEI
c       IF((cJ.eq.cI).and.(cL.eq.cK)) THEN
c                                     RESCALE_2B(L)= 0.5*ZONER+0.5*ZONEI
c       IF(cJ.eq.cK)                  RESCALE_2B(L)= 0.5*ZONER+0.5*ZONEI
c        END IF



c         write(6,*)"cI-cJ-cK-cL:",cI,cJ,cK,cL, RESCALE_2B(L)
C      write(6,'(i8,a20,4i3,F5.2,F5.2)')L,"cIcJcKcL:",cK,cL,cI,cJ
C     &,DREAL(RESCALE_2B(L)),DIMAG(RESCALE_2B(L))
c         write(6,*)L,"cI-cJ-cK-cL:",cI,cJ,cK,cL
c=========================== End Rescale
                  EndDo
                  i2=cI+1
                  EndDo
               EndDo
               EndDo
c              RESCALE_2B=RESCALE_2B/2 ! Because of 1/2 Prefactor before Wijkl
c             
c       RESCALE_2B(1)= 0.5*ZONER+0.0*ZONEI
c       RESCALE_2B(2)= 0.5 *ZONER+0.5 *ZONEI
c       RESCALE_2B(3)= 0.25 *ZONER+0.25 *ZONEI
c       RESCALE_2B(4)= ZONER+ZONEI
c       RESCALE_2B(5)= 0.5 *ZONER+0.5 *ZONEI
c       RESCALE_2B(6)= 0.25*ZONER+0.25*ZONEI
c
c       RESCALE_2B(1)= 0.5*ZONER+0.0*ZONEI
c       RESCALE_2B(2)= 0.5 *ZONER+0.5 *ZONEI
c       RESCALE_2B(3)= 0.5 *ZONER+0.5 *ZONEI
c       RESCALE_2B(4)= 0.5 *ZONER+0.0 *ZONEI
c       RESCALE_2B(5)= 0.5 *ZONER+0.5 *ZONEI
c       RESCALE_2B(6)= 0.5*ZONER+0.0*ZONEI
               MaxTrm2b=(Morb*(Morb+1)/2+1)*(Morb*(Morb+1)/2)/2
c          pause
c================================================================         
c==================== MAP Constructor
c        call MAP_CONSTRUCTOR
c================================================================         
      end subroutine Init_CI
