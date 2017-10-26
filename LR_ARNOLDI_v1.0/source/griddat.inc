C***********************************************************************
C                                                                     
C              MCTDH include file griddat.inc 
C                                                                     
C  MCTDH include file defining array sizes,  parameters etc.  
C  which are relevant for the primitive basis (grid).  
C                                                                     
C                                                                     
C***********************************************************************

C-----------------------------------------------------------------------
C VARIOUS INTEGERS DEFINING SYSTEM
C
C ndof     : total number of degrees of freedom
C feb      : number of degree of freedom described by an electronic basis
C fkoe     : dof containing coefficient information (ndof+1)
C ftime    : dof containing time information (ndof+2)
C gdim(f)  : number of grid points for dof f
C maxgdim  : Max[gdim(f)]
C maxgdim2 : Max[gdim(f)**2]
C-----------------------------------------------------------------------
      integer ndof,feb,fkoe,gdim(maxdim),maxgdim,maxgdim2,ftime

      common /pgriddat/ ndof,feb,fkoe,gdim,maxgdim,maxgdim2,ftime

C-----------------------------------------------------------------------
C HOW ELECTRONIC BASIS IS HANDLED
C
C  leb       : set to true if operator contains an electronic basis
C  lmult     : multiple state calculation
C  lmulpack  : two or more wavepackets are propagated simultaneously
C-----------------------------------------------------------------------
      logical leb,lmult,lmulpack

      common /elbopts/ leb,lmult,lmulpack

C-----------------------------------------------------------------------
C VARIOUS INTEGERS DEFINING SIZE OF DVR
C
C ortdim   : Length of array containing grid positions
C            ortdim = Sum[subdim(m)]
C dvrdim   : Length of array containing DVR representations
C            dvrdim = Sum[subdim(m)*subdim(m)]
C            where sum is over modes using DVR
C fftdim   : Length of array containing FFT representations
C            fftdim = Sum[subdim(m)]
C            where sum is over modes using FFT
C expdim   : Length of complex transformation array of EXP-DVRs
C            expdim = Sum[gdim(f)*kdim]
C            where sum is over modes using EXP-DVR
C sphdim: length of arrays related to spherical harmonics primitive
C         basis: sphdim = Sum[subdim(m)]
C-----------------------------------------------------------------------
      integer dvrord
      parameter (dvrord=2)

      integer ortdim,dvrdim,fftdim,expdim,sphdim

      common /dvrdims/ ortdim,dvrdim,fftdim,expdim,sphdim

C-----------------------------------------------------------------------
C PARAMETERS FOR ARRAY DIMENSIONS
C
C  mbaspar:    maximum number of parameters used to define bases
C-----------------------------------------------------------------------
      integer       mbaspar
      parameter     (mbaspar=6)

C-----------------------------------------------------------------------
C PARAMETERS DEFINING PRIMITIVE BASIS
C
C basis(f) : basis type for dof f
C rpbaspar(n,f): real parameters needed to define primitive basis
C ipbaspar(n,f): integer parameters needed to define primitive basis
C-----------------------------------------------------------------------
      integer basis(maxdim),ipbaspar(mbaspar,maxdim)
      real*8  rpbaspar(mbaspar,maxdim)
      real*8  xend(2,maxdim)
      logical ldvr(maxdim)
      character*(c1) modelabel(maxdim)

      common /ipbas/ basis,ipbaspar
      common /lpbas/ ldvr
      common /ppar/ rpbaspar,xend
      common /mdlbl/ modelabel

C-----------------------------------------------------------------------
C POINTERS FOR ARRAYS
C
C zort(f)    : pointer for ort matrix for dof f
C zdvr(f)    : pointer for dvr matrices (dif1mat, dif2mat, trafo) for dof f
C zfft(f)    : pointers for FFT grid matrices for dof f 
C zexp(f)    : pointers for transformation matrixes for EXP-DVR for 
C              dof f
C zsph(f)    : pointer for jsph, msph, kinsph matrices for dof f 
C-----------------------------------------------------------------------
      integer zort(maxdim),zdvr(maxdim),zfft(maxdim),zexp(maxdim),
     +        zsph(maxdim)

      common /dvrpoint/ zort,zdvr,zfft,zexp,zsph

C-----------------------------------------------------------------------
C PARAMETERS DEFINING "COMBINED" GRIDS
C
C lconm    : calculation contains combined spf_s
C nmode    : number of spf degrees of freedom
C meb      : mode including an electronic basis (feb)
C mkoe     : "mode" including Hamiltonian constants
C vgdim(f) : number of grid points for dofs preceeding f
C ngdim(f) : number of grid points for dofs following f
C subdim(m): number of grid points for mode m
C supgdim(m): number of grid points for mode m, subdim(m) = supgdim(m)^2, 
C             for densities of type I only
C npacket  : number of wavepackets to be propagated simultaneously
C nstate   : number of separate states in calculation
C sgdim    : Max[supgdim(m)] (for densities of type I only)
C sdim     : Max[subdim(m)]
C sdim2    : Max[subdim(m)**2] (for densities of type I only)
C-----------------------------------------------------------------------
      logical lconm
      integer nmode,meb,nstate,npacket,sdim,sdim2,sgdim,mkoe
      integer subdim(maxdim),supgdim(maxdim),vgdim(maxdim),ngdim(maxdim)

      common /igriddat/ nmode,meb,nstate,npacket,subdim,supgdim,vgdim,
     +                  ngdim,sdim,sdim2,sgdim,mkoe
      common /lgriddat/ lconm
 
      integer nspfdof(maxdim),spfdof(maxdim,maxdim),dofspf(maxdim)
      common /spfdofs/ nspfdof,spfdof,dofspf

C-----------------------------------------------------------------------
C Special parameters needed in various parts of the program
C
C mass(f)  : (reduced) mass for dof f
C jtot     : Total angular momentum 
C jbf      : (fixed) angular m or K quantum number. (Leg, coupled states).
C csmass   : 
C-----------------------------------------------------------------------
      real*8 mass(maxdim),jtot,jbf,csmass

      common /specpar/ mass,jtot,jbf,csmass
    
