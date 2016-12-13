C########################################################################
C                                                                       #
C                        SCHMIDTORTHO                                   #
C                                                                       #
C Library module making Schmidt-orthogonalization.                      #
C                                                                       #
C Contains:                                                             #
C       schmidtortho: Schmidt-orthogonalizes a set of complex vectors   #
C                     in column form.                                   #
C                                                                       #
C V7.0 MB                                                               #
C########################################################################


C-----------------------------------------------------------------------
C                         SCHMIDTORTHO
C
C Schmidt-orthogonalizes a set of complex vectors,
C (e.g. single-particle functions).
C The vectors are the column-vectors of psi(gdim,dim)
C-----------------------------------------------------------------------

      subroutine schmidtortho (psi,gdim,dim,ierr)

      implicit none

      integer    gdim,dim,e,e1,ierr
      complex*16 psi(gdim,dim),overlap
      real*8     norm

C-----------------------------------------------------------------------
C Schmidt-orthogonalize the set of functions.
C The orthonormalization is made twice to remove numerical inaccuracies.
C-----------------------------------------------------------------------
      ierr = 0
      do e=1,dim
         do e1=1,e-1
            call vvaxzz(psi(1,e1),psi(1,e),overlap,gdim)
            call xvxxzzs(overlap,psi(1,e1),psi(1,e),gdim)
         enddo
         call normvxz(psi(1,e),norm,gdim)
         if (norm .le. 1.0d-99 ) norm = 1.0d99
         call xvixdzo(norm,psi(1,e),gdim)
      enddo

      do e=1,dim
         do e1=1,e-1
            call vvaxzz(psi(1,e1),psi(1,e),overlap,gdim)
            call xvxxzzs(overlap,psi(1,e1),psi(1,e),gdim)
         enddo
         call normvxz(psi(1,e),norm,gdim)
         if( norm .le. 0.8d0 ) then
            ierr = e 
            return
         end if
         call xvixdzo(norm,psi(1,e),gdim)
      enddo

      return
      end


C-----------------------------------------------------------------------
C                            SCHMIDTORTHOD
C
C Schmidt-orthogonalizes (row-orthogonalization) a real quadratic
C matrix psi(dim,dim).
C
C SS 11/98
C-----------------------------------------------------------------------

      subroutine schmidtorthod (psi,dim)

      implicit none

      integer    dim,e,e1,i
      real*8 psi(dim,dim),overlap,norm
      call tranqxd(psi,dim)

C-----------------------------------------------------------------------
C Schmidt-orthogonalize the matrix.
C The orthonormalization is made twice to remove numerical inaccuracies.
C-----------------------------------------------------------------------
      do i=1,2
         do e=1,dim
            do e1=1,e-1
               call vvtxdd(psi(1,e1),psi(1,e),overlap,dim)
               call xvxxdds(overlap,psi(1,e1),psi(1,e),dim)
            enddo
            call normvxd(psi(1,e),norm,dim)
            call xvixddo(norm,psi(1,e),dim)
         enddo
      enddo

      call tranqxd(psi,dim)

      return
      end


C#######################################################################
