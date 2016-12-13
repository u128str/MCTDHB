      subroutine loopy_main(Npar,Morb,Sh_m,Sh_p,
     .                      MCNK,CI,CJ,CK,
     .                      CL,nconf,I_current_term,nthreads,
     .                      TID,RESC,WIJ,VIN,
     .                      VOUT,RIJKL,RDIM1,zrho,
     .                      zrhocnjg)

      IMPLICIT NONE


      INTEGER :: CI,CJ,CK,CL,nconf,I_current_term,RDIM1,nthreads,TID
      INTEGER , DIMENSION(0:NPAR+MORB-1,MORB-1) :: MCNK
      integer lix(Morb-1), ix 
      INTEGER   , DIMENSION(MORB+1) :: Sh_m,Sh_p
      integer depth,FromN,ToN,Npar,Morb
      complex*16 :: zrho,zrhocnjg,RESC,WIJ
      complex*16 :: VIN(nconf),VOUT(nconf,nthreads),RIJKL(RDIM1) 

      lix = 0
      ix=1
!CCCCCCCCCCCCCCC  The number of nested loops is MORB-1 in DEPTH
      depth=Morb-1
      FromN= Npar+Morb-1-Sh_m(1)
      ToN=Morb-1
      zrhocnjg=dcmplx(0.d0,0.d0)
      zrho=dcmplx(0.d0,0.d0)
      if(depth==1) ToN=Morb-1 +Sh_m(2)
      call loopy(depth, FromN,ToN, lix, 
     .           ix, Sh_m,Sh_p,Morb,
     .           NPAR,MCNK,nconf,CI,
     .           CJ,CK,CL,I_current_term,
     .           zrho,zrhocnjg,nthreads,TID,
     .           RESC,WIJ,VIN,VOUT,
     .           RIJKL,RDIM1)
      end subroutine loopy_main

      recursive subroutine loopy(depth,FromN,ToN,lix,
     . ix,Sh_m,Sh_p,Morb,
     . NPAR,MCNK,nconf,CI,
     . CJ,CK,CL,I_current_term,
     . zrho,zrhocnjg,nthreads,TID,
     . RESC,WIJ,VIN,VOUT,
     . RIJKL,RDIM1)

      IMPLICIT NONE

      
      integer :: CI,CJ,CK,CL,nconf,I_current_term,RDIM1,l,nthreads,TID
      integer, intent(inout)     :: depth
      integer   :: FromN, ix, ToN,jj,ii,Morb,k,NPAR,IND1,IND2
      integer, intent(inout)  :: lix(depth)
      INTEGER   , DIMENSION(Morb+1) :: Sh_m,Sh_p
      INTEGER , DIMENSION(DEPTH) :: JVEC
      INTEGER  , DIMENSION(Morb) ::nvecin,nvecout
      INTEGER , DIMENSION(0:NPAR+MORB-1,0:MORB-1) :: MCNK

      REAL*8 :: xbprefac

      complex*16 :: zdir,zrho,zinv,zrhocnjg,RESC,WIJ
      complex*16 :: VIN(nconf),VOUT(nconf,nthreads),RIJKL(RDIM1) 

      do jj =  FromN,ToN,-1
!CCCCCCC lix contains the J indices
         lix(ix) = jj
         if (ix .eq. depth) then
!            write(*,*) 'WTF?!'
!            print *, "NEW:", ((lix(ii)), ii = 1, depth)
!CCCCCCCCCC assign the incoming fock vector
            nvecin(1)=Npar+Morb-1-lix(1)
            do k=2,MORB-1
               nvecin(k)=lix(k-1)-lix(k)-1
            end do
            nvecin(MORB)=lix(MORB-1)-1
!CCCCCCCCC assign the corresponding J indices
            JVEC(1)=NPAR+MORB-1-(nvecin(1)-Sh_m(1)+Sh_p(1))
            do k=2, MORB-1
               JVEC(k)=JVEC(k-1)-1-(nvecin(k)-Sh_m(k)+Sh_p(k))
            end do
!            print *, "NEWj:", ((JVEC(ii)), ii = 1, depth)
!CCCCCCCCC Finally calculate the indices of the CI-Vector elements
            IND1=1
            IND2=1
            do k=1,MORB-1
               IND1=IND1+MCNK(lix(k)-1,MORB-k)    ! this is K in the org
               IND2=IND2+MCNK(JVEC(k)-1,MORB-k)   ! this is J in the org
            end do
!            print *, "NEW CI:", IND1,IND2
            nvecout=nvecin
!CCCCCCCCCCCCC assign the the prefactor according to creators and
!CCCCCCCCCCCCC annihilators
            xbprefac=1.d0
            xbprefac=xbprefac*nvecout(cL)
            nvecout(cL)=nvecout(cL)-1
            xbprefac=xbprefac*nvecout(cK)
            nvecout(cK)=nvecout(cK)-1 
            nvecout(cJ)=nvecout(cJ)+1
            xbprefac=xbprefac*nvecout(cJ)
            nvecout(cI)=nvecout(cI)+1
            xbprefac=xbprefac*nvecout(cI)
!CCCCCCCCCCCCC ASSIGN the 2b-elements
       zdir=VIN(IND1)* SQRT(xbprefac)
       zrho=zrho+Conjg(VIN(IND2))*zdir
       RIJKL(I_current_term)=RIJKL(I_current_term)+Conjg(VIN(IND2))*zdir
!        write(6,*) 'NEW-RIJKL',RIJKL(I_current_term),IND1,IND2,zdir
!        pause
!        if (abs(zdir).gt.1000000.d0) stop
c=========================== Conjugate C^_K C^_L C_I C_J   rho_IJ_KL=rho^*_KL_IJ
        zinv=VIN(IND2)* SQRT(xbprefac)
!        write(*,*) 'NEWzrho',zrho,conjg(VIN(IND2)),xbprefac
!         write(*,*) 'New-NVECOUT:',(nvecout(l),l=1,MORB-1)
        zrhocnjg=zrhocnjg+Conjg(VIN(IND1))*zinv
        VOUT(IND2,TID)=VOUT(IND2,TID)
     .+zdir* WIJ
     .     *DREAL(RESC)
        VOUT(IND1,TID)=VOUT(IND1,TID)
     .+zinv* Conjg(WIJ)
     .     *DIMAG(RESC)
!        write(*,*) 'NEW-VOUT',VOUT(IND2,TID),VOUT(IND1,TID)
         else
            FromN= lix(ix)-1-Sh_m(ix+1)
            ToN=Morb-(ix+1)
         if(ix+1==depth) ToN=Morb-ix-1 +Sh_m(ix+2)        
      call loopy(depth,FromN,ToN,lix,
     .           ix+1,Sh_m,Sh_p,Morb,
     .           NPAR,MCNK,nconf,CI,
     .           CJ,CK,CL,I_current_term,
     .           zrho,zrhocnjg,nthreads,TID,
     .           RESC,WIJ,VIN,VOUT,
     .           RIJKL,RDIM1)
         end if
      end do
      end subroutine loopy
