      subroutine  MAP_CONSTRUCTOR
      USE SHARED_DIMS
      USE CI_All 
      implicit NONE
c=========================================================
       integer::  i,j,k,l,ii
       integer::  i1,i2,i3,i4,i5,i6,i7,i8,i9,i10
       integer::  j1,j2,j3,j4,j5,j6,j7,j8,j9,j10
       integer::  n1,n2,n3,n4,n5,n6,n7,n8,n9,n10
       integer::  m1,m2,m3,m4,m5,m6,m7,m8,m9,m10
       integer::  Sh_m1,Sh_m2,Sh_m3,Sh_m4,Sh_m5
       integer::  Sh_m6,Sh_m7,Sh_m8,Sh_m9,Sh_m10
       integer::  Sh_p1,Sh_p2,Sh_p3,Sh_p4,Sh_p5
       integer::  Sh_p6,Sh_p7,Sh_p8,Sh_p9,Sh_p10
       integer::  I_current_term,Nterms,alloc_stat
       integer::  Sh_k,Sh_l
       integer::  Ncnf
c      INTEGER ::  MYID
c      COMPLEX*16, DIMENSION(Nconf) :: VIN
c       integer, ALLOCATABLE  :: MAPI(:),MAPII(:)
       integer, ALLOCATABLE  :: MAP(:)
c       integer, DIMENSION(27)  :: MAPII,MAPIII
       INTEGER, DIMENSION(Morb) :: nvecin
       
          real*8 :: CNK
c=========================================================
        ALLOCATE( MAPI(Nconf),stat=alloc_stat) 
        if(alloc_stat/=0) STOP "MEM Alloc MAPI"
        IF(ALLOCATED(MAPI).eqv..FALSE.)
     . write(6,*)" MEM for first MAPI allocation is NOT ok!",SIZE(MAPI)
        ALLOCATE( MAP(Nconf),stat=alloc_stat) 
        if(alloc_stat/=0) STOP "MEM Alloc MAPI"
        IF(ALLOCATED(MAP).eqv..FALSE.)
     . write(6,*)" MEM for first MAP allocation is NOT ok!",SIZE(MAP)
              MAP=0
              MAPI=0
c=========================================================
               Nmax=Npar
c               Nmax(1)=2
c               Nmax(1)=Npar
c               Nmax(2)=Npar
C               Nmax(3)=Npar/5
C               Nmax(4)=Npar/5
C               Nmax(3)=5
C               Nmax(4)=5
               Nmax(1)=Npar
               Nmax(2)=Npar
               Nmax(3)=Npar/2
               Nmax(4)=Npar/2
c               Nmax(3)=2
c               Nmax(4)=1
c          write(6,*)"1",MCNK(Npar+Morb-1-1,Morb-1)
c          write(6,*)"1",MCNK(Npar+Morb-1-Nmax(1)-1,Morb-1)
c          write(6,*)"2",MCNK(Npar+Morb-2-2,Morb-2)
c          write(6,*)"2",MCNK(Npar+Morb-2-Nmax(2)-2,Morb-2)
c          write(6,*)"3",MCNK(Npar+Morb-3-3,Morb-3)
c          write(6,*)"3",MCNK(Npar+Morb-3-Nmax(3)-3,Morb-3)
c          i1=Npar+Morb-1-(Nmax(1)+0)
c          i1=Npar+Morb-1-(Nmax(1)+0)
c          i2=Npar+Morb-1-(Nmax(1)+0)-1
c          i3=Npar+Morb-1-(Nmax(1)+0)-2
c          K=1+MCNK(i1-1,Morb-1)+MCNK(i2-1,Morb-2)+MCNK(i3-1,Morb-3)
c          write(6,*)"Nmax(1)000",K
c          i1=Npar+Morb-1-Npar
c          i2=Npar+Morb-1-Npar-1
c          i3=Npar+Morb-1-Npar-2
c          K=1+MCNK(i1-1,Morb-1)+MCNK(i2-1,Morb-2)+MCNK(i3-1,Morb-3)
c          write(6,*)"N000",K
c          i1=Npar+Morb-1
c          i2=Npar+Morb-1-(Nmax(2)+1)-1
c          i3=Npar+Morb-1-(Nmax(2)+1)-2
c          K=1+MCNK(i1-1,Morb-1)+MCNK(i2-1,Morb-2)+MCNK(i3-1,Morb-3)
c          write(6,*)"0Nmax(2)00",K
c          i1=Npar+Morb-1
c          i2=Npar+Morb-1-1
c          i3=Npar+Morb-1-1-1
c          K=1+MCNK(i1-1,Morb-1)+MCNK(i2-1,Morb-2)+MCNK(i3-1,Morb-3)
c          write(6,*)"last",K
c          i3=Nmax(4)+1
c          i3=Npar+Morb-1-Nmax(3)
c          K=1+MCNK(i1-1,Morb-1)+MCNK(i2-1,Morb-2)+MCNK(i3-1,Morb-3)
c          write(6,*)"last4",K



               L=0
c=========================================================
        MorbChoice: SELECT CASE (Morb)
                    CASE (1)
c                    CALL GetCIJKL1body_1(MYID,VIN)
                    CASE (2)
c                    CALL GetCIJKL1body_2(MYID,VIN)
                    CASE (3)
c                    CALL GetCIJKL1body_3(MYID,VIN)
                    CASE (4)
c=============================================================
                      m1=MAX0(Morb-1+Npar-Nmax(1),Morb-1) ! fix the max n on n3
       i1loop:    Do i1=Npar+Morb-1    ,m1  ,-1  
          nvecin(1)=Npar+Morb-1-i1
                     m2=MAX0(i1-1-Nmax(2),Morb-2)         ! fix the max n on n3
       i2loop:    Do i2=i1-1 ,m2 ,-1 
          nvecin(2)=i1-i2-1
                     m3=MAX0(i2-1-Nmax(3),Morb-3)         ! fix the max n on n3
                     m4=MIN0(i2-1,1+Nmax(4))              ! fix the max n on n4
       i3loop:    Do i3=m4 ,m3    ,-1
          nvecin(3)=i2-i3-1
          nvecin(4)=i3-1
          K=1+MCNK(i1-1,Morb-1)+MCNK(i2-1,Morb-2)+MCNK(i3-1,Morb-3)
           L=L+1
           MAPI(L)=K
           EndDo I3loop
           EndDo I2loop
           EndDo I1loop
           Nconf=L
c=============================================================
                    CASE (5)
c                    CALL GetCIJKL1body_5(MYID,VIN)
                    CASE (6)
c                    CALL GetCIJKL1body_6(MYID,VIN)
                    CASE (7)
c                    CALL GetCIJKL1body_7(MYID,VIN)
                    CASE (8)
c                    CALL GetCIJKL1body_8(MYID,VIN)
                    CASE (9:100)
       write(6,*)" Still NOT implemented !!!!!!!!!"
                    CASE DEFAULT
       write(6,*)" Something wrong in Number of orbitals!!!!!MAP_CNSTR"
                    END SELECT MorbChoice
c=============================================================
        write(6,*) "New Reduced CI size Nconf=",Nconf
        Ncnf=NINT(CNK(Npar+Morb-1,Npar))
        write(6,*)"Nconf",Nconf,"out of",Ncnf
        write(6,*)"Constrained Nmax used:",(Nmax(i),i=1,Morb)

c      if(SIZE(MAP).lt.Nconf) STOP 
c     . "MAP_CONSTRUCTOR: DIM(MAP) is not big enough- should be >=Nconf"
c                ALLOCATE( MAP(Nconf),stat=alloc_stat)
c                if(alloc_stat/=0) STOP "2" 
c                ALLOCATE( MAPII(Nconf),stat=alloc_stat)
c                if(alloc_stat/=0) STOP "3" 
c        IF(ALLOCATED(MAP).eqv..FALSE.)
c     . write(6,*)" MEM for second MAP allocation is NOT ok!",SIZE(MAP)
c        IF(ALLOCATED(MAPII).eqv..FALSE.)
c     . write(6,*)" MEM for second MAP allocation is NOT ok!",SIZE(MAPII)
c=============================================================
              L=1
              do i=Nconf,1,-1
              MAP(i)=MAPI(L)
              L=L+1
              enddo
              MAPI=0
c         MAPI(MAP)=(/(j,j=1,Nconf)/)
              Do j=1,Nconf
              i=MAP(j)
              MAPI(i)=j
              End Do

c              write(6,*)"MAP",(MAP(i),i=1,Nconf)
 
c              MAP=MAPII(MAPII)
c              MAP(MAP)=(/(i,i=1,Nconf)/)
c                 Do i=1,Nconf
c                MAP(MAP(i))=i
c                 Enddo
c             write(6,*)"MAP",(MAP(MAP(i)),i=1,Nconf)
c             write(6,*)"MAP"
c=============================================================
       deallocate(MAP)
c       deallocate(MAPII)
c       stop
      end subroutine MAP_CONSTRUCTOR
