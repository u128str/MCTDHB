C=============== INtegral Hij Wijkl Evaluation for BOSE-HUBBARD MODELS
       subroutine Get_BH_integrals
       USE SHARED_DIMS
       USE rR_hW
       implicit NONE

       integer :: icntr, Iorb, Jorb, I
       integer :: P, cL, cK, cJ, cI
       real*8 :: rfactor, j_BH, U_BH, Delta_BH
       character*4 :: BH_type
c==========================================

        J_BH=1.d0
        U_BH=0.25d0
        Delta_BH=0.d0
 
        J_BH=0.0223340316036491d0
        U_BH=-0.0050537669574329d0
        Delta_BH=0.4954580796706878d0
 
c       J_BH=0.0000000000001d0
c       U_BH=2/380.0d0
c       Delta_BH=0.0d0
 
        Do I=1,MaxTrm1b
        P=TERM_INDEX_1B(I)
c================ Unpack cI cJ cK cL from P
        cK= INT(P/100)
        cJ= P-cK*100
        HIJ(I)=(0d0,0d0)
        IF(cK.eq.cJ+1) HIJ(I)=-J_BH           ! Comes from  -J ( C^*_2 C_1+ C^1 C_2) J==J_BH
        IF(cK.eq.cJ) HIJ(I)=Delta_BH
      write (6,'(a7,i2,a,i2,a3,2F16.10)')"1-B: H(",cJ,",",cK,")=",HIJ(I)
        EnDdo
c========== COMMENT the following line for BOSE_HUBBARD WITH OPEN BOUNDARY CONDITIONS
!        HIJ(Morb)=-J_BH ! uncomment it for even number of cites 2,4,6,, for odd 3,5,7, it does-not work 
       write (6,'(a7,i2,a,i2,a3,2F16.10)')
     . "1-B: H(",1,",",Morb,")=",HIJ(I)
        
        icntr=1
        DO Iorb=1,Rdim
        DO Jorb=Iorb,Rdim

        P=TERM_INDEX_2B(icntr)
c================ Unpack cI cJ cK cL from P
          cL= INT(P/1000000)
          cK= INT((P-cL*1000000)/10000)
          cJ= INT((P-cL*1000000-cK*10000)/100)
          cI= P-cL*1000000-cK*10000-cJ*100

        WIJKL(icntr)=(0d0,0d0)

ccccccccccc Remembert that the Hamiltonain is 1/2 WIJKL b^\dag_i b^\dag_j b_k b_l

        if (cI.eq.cJ.and.cJ.eq.cK.and.cK.eq.cL) WIJKL(icntr)=U_BH ! Comes from 1/2 U n1n1

       write (6,'(a7,i2,a,i2,a,i2,a,i2,a3,2F16.10)')
     . "2-B: W(",cI,",",cJ,",",cK,",",cL,")=",WIJKL(icntr)
c       write (6,*) ZRIJKL(icntr)
c       write (6,*)

        icntr=icntr+1
        EndDO
        EndDO



  
      end subroutine Get_BH_integrals
