       subroutine print_cic(time,VIN,path)
       USE   SHARED_DIMS
       USE   rR_hW
       USE   CI_All
       USE   W_INTERPARTICLE
       USE   DVR_ALL

c       implicit real*8 (a-h,o-z)
c       implicit integer (i-n)
       COMPLEX*16, DIMENSION(Nconf) :: VIN
       COMPLEX*16, DIMENSION(1000)  :: v
       real*8 time,t_total,ww
       INTEGER :: jj,Ndim,ind
       character*18 lname
       character*10 path
       integer, DIMENSION(Morb)   :: Ivec
       integer, DIMENSION(Morb)   :: Nvec

c=================================================
        ndim=Nconf
       t_total=time
       if(t_total.lt.99999.99999999D0) then
       Write( lname, '(F10.4,a8)') t_total,"coef.dat"
       endif
       if(t_total.lt.9999.99999999D0) then
       Write( lname, '(F10.5,a8)') t_total,"coef.dat"
       endif
       if(t_total.lt.999.99999999D0) then
       Write( lname, '(F10.6,a8)') t_total,"coef.dat"
       endif
       if(t_total.lt. 99.99999999D0) then
       Write( lname, '(F10.7,a8)') t_total,"coef.dat"
       endif
       if(t_total.lt.  9.99999999D0) then
       Write( lname, '(F10.8,a8)') t_total,"coef.dat"
       endif

      open(unit=111,file=trim(path)//'/'//lname,form='formatted')
      write(111,*)"          "
!      write(111,'(a1,a28,138(I26,a3))')"#","",(jj,"  ",jj=2,10)
       
!       do ind=1,ndim
!        IF(ABS(VIN(ind)).ge.0.5d-14) THEN
!      write(111,2222) ind,"  "
!     $,DREAL(VIN(ind)),"   ",DIMAG(VIN(ind)),"   "
!     $,t_total
!        ENDIF
!       enddo

      write(111,2221)"n","   ","REAL(C(n))"
     $            ,"   ","IMAG(C(n))","   ",
     $            (jj,"  ",jj=1,Morb)

2221  format((a26,a3,2(a26,a3),14(i6,a3)))

      do ind=1,ndim
          IF(ABS(VIN(ind)).ge.0.5d-14) THEN
          Ivec=0
          Nvec=0
          call GET_Ind_vs_ii(Ind,Npar,Morb,Ivec,Nvec)
            write(111,2222) ind,"   "
     $           ,DREAL(VIN(ind)),"   ",DIMAG(VIN(ind)),"   "
     $           ,(Nvec(jj),"  ",jj=1,Morb)
         ENDIF
      enddo


      close(111)
      return
2222  format((i26,a3,2(F26.16,a3),14(i6,a3)))


!       close(111)
c       IF(MORB.eq.2)  CALL U_Transform(time,VIN)
!       return
!2222  format((i26,a3,10(F26.16,a3)))
       end SUBROUTINE print_cic

