************************************************************************
*
* SUBROUTINE RK8
*
* MCTDH interface to the call of the Runge-Kutta (order 8) integrator
* from the ODESIM project.
*
************************************************************************
*
* Input:
*     psi         : initial values of the dynamic variables
*     dtpsi       : initial time derivative of the dyn. vars
*     psidim      : number of dyn. vars
*     intperiod   : length of integration interval
*     abstime     : initial time
*     initstep    : initial stepsize
*     accuracy    : parameter controlling the accuracy
*     auxpsi      : storage for intermediate results, workspace...
*                   size: 12*psidim; can be uninitialized
*     func        : routine for evaluating the RHS
*     mc,mr,mi,ml : the MCTDH-arrays
*     lpost       : whether to call a post-processor
*     postfn      : post-processing routine
*
* Output:
*     psi         : final values of the dyn. vars
*     dtpsi       : undetermined (i.e. it gets overwritten!)
*     initstep    : optimal stepsize at end of integration
*     error       : errorflag; 0 means OK (currently it's always 0!)
*
* FO, 08/2002
*
************************************************************************

      subroutine rk8(psi,dtpsi,psidim,intperiod,abstime,initstep,
     $               accuracy,error,auxpsi,func,
     $               lpost,postfn,rksteps)


      implicit none

#include "global.inc"
      integer rk5steps,rk5accpt,rk5rejct,rk5fncns,
     $        rk8steps,rk8accpt,rk8rejct,rk8fncns



      integer    psidim,error,rksteps
      real*8     intperiod,abstime,initstep,accuracy,wt(psidim)
      complex*16 psi(psidim),dtpsi(psidim),auxpsi(psidim*12)
      external   func,postfn
      logical    lpost

      real*8     tend,hmax,h,hopt,safe,beta,facdec,facinc,
     $           rerr,time
      integer    intgfl,ierr
      logical    initialized

      parameter (
     $     safe   = 0.9d0,
     $     beta   = 0.04d0,
     $     facdec = 0.333333d0,
     $     facinc = 6.0d0 )

CCCCCC   Axel threw away the wheel of the bicycle
!      save initialized
!      data initialized /.false./

*     initialize statistics on first call
         initialized = .FALSE.
      if (.not. initialized) then
         rk8steps = 0
         rk8accpt = 0
         rk8rejct = 0
         rk8fncns = 0
         initialized = .true.
      endif

*     set parameters for call to core integrator
      h       = initstep
      time    = abstime
      tend    = abstime+intperiod
      hmax    = intperiod
      intgfl  = 0


*     call core integrator
      call OD85CR(
     $     func, psidim, psi, time, tend, hmax, h, hopt, accuracy,
     $     macheps, safe, beta, facdec, facinc,
     $     auxpsi(1), auxpsi(psidim+1), auxpsi(2*psidim+1),
     $     auxpsi(3*psidim+1), dtpsi, auxpsi(4*psidim+1),
     $     auxpsi(5*psidim+1), auxpsi(6*psidim+1), auxpsi(7*psidim+1),
     $     auxpsi(8*psidim+1), auxpsi(9*psidim+1), auxpsi(10*psidim+1),
     $     auxpsi(11*psidim+1), wt , rk8fncns, rk8steps,
     $     rk8accpt, rk8rejct, intgfl, ierr, rerr,
     $     lpost, postfn)

*     check error condition
      if (intgfl .ne. 1) then
         write (6,*) 'RUNGEKUTTA order 8 Problem: Errorflag:', intgfl  
      else
         error = 0
      endif


*     return the optimal stepsize and number of steps
      initstep = hopt
      rksteps=rk8accpt+rk8rejct
      return

 100  format ('Error number ',i2,', parameter ',f16.10)

      end
