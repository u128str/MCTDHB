
SUBROUTINE  Get_RhoIJKL(Rho_IJKL)
USE   PASS_ARG
USE   SHARED_DIMS
USE   rR_hW
USE   CI_All
USE   W_INTERPARTICLE
USE   DVR_ALL
USE   Parallel_Orb
IMPLICIT NONE
INTEGER :: cK,cJ,cL,cI,L2B,P,I
DOUBLE COMPLEX, Dimension(Morb,Morb,Morb,Morb) :: Rho_IJKL

Rho_ijkl=dcmplx(0.d0,0.d0)
L2B=Rdim*(Rdim+1)/2
DO I=1,L2B
   P=TERM_INDEX_2B(I)
   cL= INT(P/1000000)
   cK= INT((P-cL*1000000)/10000)
   cJ= INT((P-cL*1000000-cK*10000)/100)
   cI= P-cL*1000000-cK*10000-cJ*100 
   Rho_ijkl(cI,cJ,cK,cL)=ZRIJKL(I) 
   Rho_ijkl(cI,cJ,cL,cK)=ZRIJKL(I) 
   Rho_ijkl(cJ,cI,cK,cL)=ZRIJKL(I) 
   Rho_ijkl(cJ,cI,cL,cK)=ZRIJKL(I) 
   IF ((CI .NE. CK) .OR. (CJ .NE. CL)) THEN
      Rho_ijkl(cK,cL,cI,cJ)=Conjg(ZRIJKL(I))
      Rho_ijkl(cK,cL,cJ,cI)=Conjg(ZRIJKL(I))
      Rho_ijkl(cL,cK,cI,cJ)=Conjg(ZRIJKL(I))
      Rho_ijkl(cL,cK,cJ,cI)=Conjg(ZRIJKL(I))
   ENDIF    
END DO


END SUBROUTINE
