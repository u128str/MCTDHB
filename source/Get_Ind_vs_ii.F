c=================================================
      SUBROUTINE GET_Ind_vs_ii(Ind0,N,M,Ivec,Nvec)
      USE CI_All
      implicit NONE
      real*8 :: CNK
      integer:: N,K,Ind,Ind0,j,M,MaxI
      integer::  i,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10
      integer, DIMENSION(100) :: Ivec
      integer, DIMENSION(100)   :: Nvec
      external CNK
c99          ind0=ind0-1
c         write(6,*)"Input Ind",IND0
c         write(6,*)"N",N
c         write(6,*)"M",M
         IND=Ind0
         Ivec=0 
         Nvec=0 
         MaxI=M+N-1 
              do i=1,M-1
      qq: Do j=M-i,MaxI
      Ivec(i)=j
c      write(6,*)"i=",i,"Ii",Ivec(i),"ind=",IND,"C()=",(MCNK(j-1,M-i))
      IF(Ind.le.MCNK(j-1,M-i)) then 
c      IF(Ind.le.NINT(CNK(j-1,M-i))) then 
      Ivec(i)=j-1
      goto 1
      endif
      End Do  qq
1      Ind=Ind-MCNK(Ivec(i)-1,M-i)
      MaxI=Ivec(i)
c      write(6,*)"i=",i,"Ii",Ivec(i),"New ind=",IND
c      if(Ind.eq.0)  goto 2
              End do

2             Nvec(1)=N+M-1-Ivec(1)
              Do i=2,M-1
              Nvec(i)=Ivec(i-1)-Ivec(i)-1
              Enddo
              Nvec(M)=Ivec(M-1)-1
              
          
c         write(6,*)INd0," Ivec", Ivec
c         write(6,*)INd0,"Nvec", Nvec
c           if (ind0.ne.1) goto 99
      end  subroutine  GET_Ind_vs_ii

       subroutine GetInd(N,M,Nvec,ind)
       USE CI_All
       implicit NONE
c      real*8 :: CNK
       integer:: N,K,Ind,Ind0,j,M,MaxI
       integer:: i,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10
       integer, DIMENSION(100)   :: Nvec
       integer, DIMENSION(100)   :: Ivec
c      external CNK
!              write(6,*)"In GetInd",N,M
!              write(6,*)"In GetInd",Nvec(1:M)
            Ivec(1)=N+M-1-Nvec(1)
            Do i=2,M-1
            Ivec(i)=Ivec(i-1)-Nvec(i)-1
            Enddo
            ind=1
            Do i=1,M-1
            Ind=Ind+MCNK(Ivec(i)-1,M-i)
            EndDo
      end  subroutine GetInd
