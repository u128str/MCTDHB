C#######################################################################
C
C#######################################################################


C-----------------------------------------------------------------------
C PARAMETERS FOR ARRAY DIMENSIONS
C
C  sbaspar:    maximum number of parameters used to define bases
C-----------------------------------------------------------------------
      integer       sbaspar
      parameter     (sbaspar=6)

C-----------------------------------------------------------------------
C array sizes
C
C spf1ddim: length of spf1d array for 1-D "spf" functions
C-----------------------------------------------------------------------
      integer spf1ddim

      common /wf1ddim/ spf1ddim

C-----------------------------------------------------------------------
C pointers for arrays for density operators of type I
C
C dzspf1d(f,s,t) : pointers for spf1d array for 1-D "spf" functions
C ddimf(f,s,t)   : no. of 1d spfs to be generated
C dzpsi1d(s,t)   : pointer to state wavefunction of 1d spfs
C-----------------------------------------------------------------------
      integer dzspf1d(maxdim,maxsta,maxsta),ddimf(maxdim,maxsta,maxsta),
     +        dzpsi1d(maxsta,maxsta)

      common /do1dpoint/ dzspf1d,ddimf,dzpsi1d

C-----------------------------------------------------------------------
C pointers for arrays
C
C zspf1d(f,s) : pointers for spf1d array for 1-D "spf" functions
C dimf(f,s)   : no. of 1d spfs to be generated
C-----------------------------------------------------------------------
      integer zspf1d(maxdim,maxsta),dimf(maxdim,maxsta)

      common /wf1dpoint/ zspf1d,dimf

C-----------------------------------------------------------------------
C FLAGS TO CONTROL INWF OPTIONS
C
C  inwfop : no. of operator used to modify an initial wavepacket.
C  lad    : adiabatic corrections
C  ldia   : diabatic corrections
C  lhh2   : special (a)diabatic corrections for H + H2, H + D2
C  lbkmp2 : special (a)diabatic corrections with BKMP2-PES
C  ledstr : calculate adiabatically corrected energy distribution
C  lcallwf : call the user supplied subroutine 'initWF' to build the
C            initial wavefunction.
C  lbuildwf: the initial wavefunction is build by genphi and gencoeff.
C  lorthopsi: wavefunction is Schmidt-orthogonalised at start of restart.
C  lveigen  : the eigenvalues and vectors of a 1d-operator are
C             written to the file veigen_modelabel.
C  lnodirect: Switches off the direct calculation of operate (inwfop).
C             Only the iterative calculation is performed.  
C-----------------------------------------------------------------------
      logical ladiab,lad,ldia,ledstr,lhh2,lcallwf,lbuildwf,lorthopsi,
     +        lveigen,lbkmp2,lnodirect

      common /linwfopt/ ladiab,lad,ldia,ledstr,lhh2,lcallwf,lbuildwf,
     +        lorthopsi,lveigen,lbkmp2,lnodirect

C-----------------------------------------------------------------------
C FLAGS TO CONTROL THE OPERATION OF AN OPERATOR ON THE INWF
C  inwfop : the operator(s) used to generate an initial wavepacket.
C  inwfopnit : Number of iterations (default=10).
C  inwfoptol : Error tolerance, default = 1.0d-8.
C-----------------------------------------------------------------------
      integer inwfop(16),inwfopnit
      real*8  inwfoptol

      common /iinwfop/ inwfop,inwfopnit
      common /rinwfop/ inwfoptol

C-----------------------------------------------------------------------
C addim    : dim(2,1)*dim(3,1)
C azahl  : number of a-coefficients to be read
C iza    : input file line for start of a-coefficient list
C-----------------------------------------------------------------------
      integer addim,azahl,iza,izcll,czahl

      common /inwfdim/ addim,azahl,iza,izcll,czahl

C-----------------------------------------------------------------------
C PARAMETERS NEEDED BY VEIGEN AND ADIABATIC ROUTINES
C-----------------------------------------------------------------------
      integer maxedim
      parameter (maxedim=100)
      real*8 e0eigen,eeigen(maxedim)

      common /reeigen/ e0eigen,eeigen

C-----------------------------------------------------------------------
C PARAMETERS DEFINING SINGLE-PARTICLE FUNCTION BASIS
C
C-----------------------------------------------------------------------
      integer sbasis(maxdim),isbaspar(sbaspar,maxdim,maxsta*maxsta)
      real*8  rsbaspar(sbaspar,maxdim,maxsta*maxsta)

      common /sbasi/ sbasis,isbaspar
      common /sbasr/ rsbaspar

C-----------------------------------------------------------------------
C common block correction routines. (needed in adiab routines and LSTH)
C-----------------------------------------------------------------------
      real*8  tfac
      integer ftr,  statetr
      common /adiabr/ tfac
      common /adiabi/ ftr,  statetr


