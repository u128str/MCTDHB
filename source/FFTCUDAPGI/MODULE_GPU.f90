!====================== From example
! This module provides a simple facility for changing between single
! and double precision

module precision_m
  integer, parameter, public :: singlePrecision = kind(0.0) ! Single precision
  integer, parameter, public :: doublePrecision = kind(0.0d0) ! Double precision
  
  ! Comment out one of the lines below
  !integer, parameter, public :: fp_kind = singlePrecision
   integer, parameter, public :: fp_kind = doublePrecision
end module precision_m
!=======================================
!
! Define the INTERFACES to the NVIDIA CUFFT routines
!

module cufft_m

  integer, parameter, public :: CUFFT_FORWARD = -1
  integer, parameter, public :: CUFFT_INVERSE = 1
  integer, parameter, public :: CUFFT_R2C = Z'2a' ! Real to Complex (interleaved)
  integer, parameter, public :: CUFFT_C2R = Z'2c' ! Complex (interleaved) to Real
  integer, parameter, public :: CUFFT_C2C = Z'29' ! Complex to Complex, interleaved
  integer, parameter, public :: CUFFT_D2Z = Z'6a' ! Double to Double-Complex
  integer, parameter, public :: CUFFT_Z2D = Z'6c' ! Double-Complex to Double
  integer, parameter, public :: CUFFT_Z2Z = Z'69' ! Double-Complex to Double-Complex

  ! -----------------
  ! cufftPlan[1|2|3]d
  ! -----------------

  interface cufftPlan1d
     subroutine cufftPlan1d(plan, nx, type, batch) bind(C,name='cufftPlan1d')
       use iso_c_binding
       integer(c_int):: plan
       integer(c_int),value:: nx, batch,type
     end subroutine cufftPlan1d
  end interface

  interface cufftPlan2d
     subroutine cufftPlan2d(plan, nx, ny, type) bind(C,name='cufftPlan2d')
       use iso_c_binding
       integer (c_int) :: plan
       integer (c_int), value :: nx, ny, type
     end subroutine cufftPlan2d
  end interface

  interface cufftPlan3d
     subroutine cufftPlan3d(plan, nx, ny, nz, type) bind(C,name='cufftPlan3d')
       use iso_c_binding
       integer (c_int) :: plan
       integer (c_int), value :: nx, ny, nz, type
     end subroutine cufftPlan3d
  end interface


!interface cufftPlanMany
!     subroutine cufftPlanMany(plan, rank, nx,&
!                                  in, is, id,&
!                                  on, os, od,&
!                              type, batch) &
!      bind(C,name='cufftPlanMany')
!      use iso_c_binding
!      integer(c_int):: plan
!      integer(c_int),value:: rank
!      integer(c_int):: nx, in, on
!      integer(c_int),value:: is, id
!      integer(c_int),value:: os, od, type, batch
!     end subroutine cufftPlanMany
! end interface 

interface cufftPlanMany
subroutine cufftPlanMany(plan, rank, n, &
                         inembed, istride, idist, &
                         onembed, ostride, odist, &
                         type, batch) &
 bind(C,name='cufftPlanMany')
use iso_c_binding
integer(c_int):: plan
integer(c_int),value:: rank, type, batch
integer(c_int),value:: istride, idist, ostride, odist
integer(c_int):: n(*)
integer(c_int):: inembed(*), onembed(*)
end subroutine cufftPlanMany
end interface 




  ! ------------
  ! cufftDestroy
  ! ------------

  interface cufftDestroy
     subroutine cufftDestroy(plan) bind(C,name='cufftDestroy')
       use iso_c_binding
       integer(c_int),value:: plan
     end subroutine cufftDestroy
  end interface

  ! A generic interface such as "cufftExec" which contained all the R2C, C2C, ...
  ! subroutine variants would work fine as long as transforms were not done in place.
  ! The types and kinds associated with idata and odata would determine the 
  ! correct routine.
  !
  ! However, this appoach would break down when performing in-place transforms.
  ! For example, suppose one wanted to perform a real to complex transform in-place.
  ! The soubroutine call:
  !    call cufftExec(plan, data, data) 
  ! would result in a compile time error if "data" were of type real or
  ! the cufftExecC2C would be executed if "data" were of type complex 

  ! So, we define a separate interface for each R2C, C2C, ... variant and 
  ! use the "!dir$ ignore_tkr" directive to ignore compile-time checks of 
  ! data type, kind, and ranks.

  ! ----------------------------------
  ! cufftExec[C2C|R2C|C2R|Z2Z|Z2D|D2Z]
  ! ----------------------------------

  interface cufftExecC2C
     subroutine cufftExecC2C(plan, idata, odata, direction) &
          & bind(C,name='cufftExecC2C')
       use iso_c_binding
       use precision_m
       integer(c_int), value:: plan
       !pgi$ ignore_tkr idata, odata
       complex(singlePrecision), device :: idata(*), odata(*)
       integer(c_int), value:: direction
     end subroutine cufftExecC2C
  end interface

  interface cufftExecR2C
     subroutine cufftExecR2C(plan, idata, odata) bind(C,name='cufftExecR2C')
       use iso_c_binding
       use precision_m
       integer (c_int), value :: plan
       !pgi$ ignore_tkr idata, odata
       real (singlePrecision), device :: idata(*)
       complex (singlePrecision), device :: odata(*)
     end subroutine cufftExecR2C
  end interface

  interface cufftExecC2R
     subroutine cufftExecC2R(plan, idata, odata) bind(C,name='cufftExecC2R')
       use iso_c_binding
       use precision_m
       integer (c_int), value :: plan
       !pgi$ ignore_tkr idata, odata
       complex (singlePrecision), device :: idata(*)
       real (singlePrecision), device :: odata(*)
     end subroutine cufftExecC2R
  end interface
        
  interface cufftExecZ2Z
     subroutine cufftExecZ2Z(plan, idata, odata, direction) &
          & bind(C,name='cufftExecZ2Z')
       use iso_c_binding
       use precision_m
       integer(c_int),value:: plan
       !pgi$ ignore_tkr idata, odata
       complex(doublePrecision), device:: idata(*), odata(*)
       integer(c_int),value:: direction
     end subroutine cufftExecZ2Z
  end interface

  interface cufftExecD2Z
     subroutine cufftExecD2Z(plan, idata, odata) bind(C,name='cufftExecD2Z')
       use iso_c_binding
       use precision_m
       integer (c_int), value :: plan
       !pgi$ ignore_tkr idata, odata
       real (doublePrecision), device :: idata(*)
       complex (doublePrecision), device :: odata(*)
     end subroutine cufftExecD2Z
  end interface

  interface cufftExecZ2D
     subroutine cufftExecZ2D(plan, idata, odata) bind(C,name='cufftExecZ2D')
       use iso_c_binding
       use precision_m
       integer (c_int), value :: plan
       !pgi$ ignore_tkr idata, odata
       complex (doublePrecision), device :: idata(*)
       real (doublePrecision), device :: odata(*)
     end subroutine cufftExecZ2D
  end interface

  interface cufftExec
     subroutine cufftExec(plan, transform, idata, odata, direction)
       integer :: plan, transform
       !pgi$ ignore_tkr idata, odata
       real, device :: idata(*), odata(*)
       integer, optional :: direction
     end subroutine cufftExec
  end interface

end module cufft_m



! A general subroutine that uses the transform type passed in "transform" to 
! determine which routine to call.  The interfaces for the routines specified
! above use the "!dir$ ignore_tkr" directive to remove checking of
! type, kind, and rank to facilitate certain in-place transforms (eg. cufftExecR2C)
! and multidimensional arrays

subroutine cufftExec(plan, transform, idata, odata, direction)
  use cufft_m
  implicit none

  integer :: plan, transform
  real, device :: idata(*), odata(*)
  integer, optional :: direction
   

  select case (transform)
  case (CUFFT_R2C) 
     call cufftExecR2C(plan, idata, odata)
  case (CUFFT_C2R) 
     call cufftExecC2R(plan, idata, odata)
  case (CUFFT_C2C)
     if (.not. present(direction)) then
        write(*,*) 'Error in cufft: C2C transform called without specifying direction'
        stop
     end if
     call cufftExecC2C(plan, idata, odata, direction)
  case (CUFFT_D2Z) 
     call cufftExecD2Z(plan, idata, odata)
  case (CUFFT_Z2D) 
     call cufftExecZ2D(plan, idata, odata)
  case (CUFFT_Z2Z)
     if (.not. present(direction)) then
        write(*,*) 'Error in cufft: Z2Z transform called without specifying direction'
        stop
     end if
     call cufftExecZ2Z(plan, idata, odata, direction)  
  case default
     write(*,*) 'Invalid transform type passed to cufftExec'
  end select
end subroutine cufftExec

