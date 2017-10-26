!          Complex function  RhoAll returns the valuse desired element of R_IJKL  for
!          any given set of K,S,Q,L
        COMPLEX*16 Function RhoAll(K,S,Q,L) 
        USE SHARED_DIMS
        USE rR_hW
        implicit NONE
c        COMPLEX*16,INTENT(IN), DIMENSION(3000) :: ZRIJKL
        INTEGER,INTENT(IN) :: K,S,Q,L
        INTEGER:: K1,S1,Q1,L1,I,J,KSadr,QLadr,Nadr,DNSDim
        DNSDim=Morb*(Morb+1)/2
        K1=K
        S1=S
        Q1=Q
        L1=L
        IF(K-S.gt.0) THEN 
        K1=S
        S1=K 
        ENDIF
        KSadr=Morb*(K1-1)-(K1-1)*(K1-2)/2+(S1-K1)+1
        IF(Q-L.gt.0) THEN 
        Q1=L
        L1=Q
        ENDIF
        QLadr=Morb*(Q1-1)-(Q1-1)*(Q1-2)/2+(L1-Q1)+1

        IF(KSadr-QLadr.le.0) THEN 
        Nadr=DNSDim*(KSadr-1)-(KSadr-1)*(KSadr-2)/2+(QLadr-KSadr)+1
        Rhoall=ZRIJKL(Nadr)
c        write(6,*)"Normal RhoALL",K,S,Q,L
        ELSE
        Nadr=DNSDim*(QLadr-1)-(QLadr-1)*(QLadr-2)/2+(KSadr-QLadr)+1
        Rhoall=Conjg(ZRIJKL(Nadr))
c        write(6,*)"CONJGE RhoALL",K,S,Q,L
        ENDIF
        END FUNCTION RhoAll

        INTEGER Function Nadr(K,S,Q,L,Morb) 
        implicit NONE
        INTEGER,INTENT(IN) :: K,S,Q,L,Morb 
        INTEGER :: K1,S1,Q1,L1,I,J,KSadr,QLadr,DNSDim
c        write(6,*)"IN  KSQL",K,S,Q,L
        DNSDim=Morb*(Morb+1)/2
        K1=K
        S1=S
        Q1=Q
        L1=L
        IF(K-S.gt.0) THEN 
        K1=S
        S1=K 
        ENDIF
        KSadr=Morb*(K1-1)-(K1-1)*(K1-2)/2+(S1-K1)+1
        IF(Q-L.gt.0) THEN 
        Q1=L
        L1=Q
        ENDIF
        QLadr=Morb*(Q1-1)-(Q1-1)*(Q1-2)/2+(L1-Q1)+1
c        write(6,*)"KSQL",KSadr,QLadr

        IF(KSadr-QLadr.le.0) THEN 
        Nadr=DNSDim*(KSadr-1)-(KSadr-1)*(KSadr-2)/2+(QLadr-KSadr)+1
        ELSE
        Nadr=DNSDim*(QLadr-1)-(QLadr-1)*(QLadr-2)/2+(KSadr-QLadr)+1
        Nadr=-Nadr
        ENDIF
c        write(6,*)"OUT KSQL",K1,S1,Q1,L1
        END FUNCTION Nadr

        INTEGER Function NadrW(K,S,Q,L,Morb) 
        implicit NONE
        INTEGER,INTENT(IN) :: K,S,Q,L,Morb 
        INTEGER :: K1,S1,Q1,L1,I,J,KSadr,QLadr,WDNSDim
c        write(6,*)"IN  KSQL",K,S,Q,L
        WDNSDim=Morb**2*(Morb**2+1)/2
        K1=K
        S1=S
        Q1=Q
        L1=L
        KSadr=Morb*(K1-1)+S1
        QLadr=Morb*(Q1-1)+L1
c        write(6,*)"KSQL",KSadr,QLadr
        IF(KSadr-QLadr.le.0) THEN 
        NadrW=WDNSDim*(KSadr-1)+QLadr
        ELSE
        NadrW=WDNSDim*(QLadr-1)+KSadr
        NadrW=-NadrW
        ENDIF
c        write(6,*)"OUT KSQL",K1,S1,Q1,L1
        END FUNCTION NadrW
