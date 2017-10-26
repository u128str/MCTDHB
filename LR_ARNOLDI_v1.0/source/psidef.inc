C***********************************************************************
C                                                                     
C              MCTDH include file   MAIN 
C                                                                     
C  MCTDH include file defining array sizes, common blocks etc. common 
C  to several parts of the program.  
C  Declarations independent of modes and states should be placed
C  in include file "global".     
C                                                                     
C  11/96 AJ                                                         
C                                                                     
C***********************************************************************

C-----------------------------------------------------------------------
C VARIOUS INTEGERS DEFINING SIZE OF WAVEFUNCTION
C
C dgldim   : Length of array containing wavefunction
C            dgldim = block + Sum[dim(m)*subdim(m)]
C dmatdim  : Length of array containing the density matrices
C            dmatdim = Sum[dim(m)*dim(m)]
C d3matdim : Length of array containing dicht3 matrix
C            d3matdim = Sum[2*dim(m)]
C adim     : Cumulative length of all A vectors.
C            adim = Sum[block(s)]
C maxspf   : Largest number of single-particle functions
C            maxspf = Max_(m,s) dim(m,s)
C maxjindx : maximum no. of J indices needed in slected CI calculation.
C-----------------------------------------------------------------------
      integer dgldim,dmatdim,d3matdim,adim,maxspf,maxjindx

      common /psidims/ dgldim,dmatdim,d3matdim,adim,maxspf,maxjindx

C-----------------------------------------------------------------------
C ARRAY DIMENSIONS DEPENDENT ON NMODE
C
C PhiDim(M): Length of phi-vector for mode m
C            PhiDim(M) = Sum[Dim(S,M)*SubDim(M)]
C maxphidim: Size of largest phi vector
C            maxphidim = max[phidim(m)]
C totphidim: overall size of single-particle functions
C            totphidim = Sum_m [phidim(m)] = dgldim-adim
C-----------------------------------------------------------------------
      integer phidim(maxdim),maxphidim,totphidim
      
      common /marraydim/ phidim,maxphidim,totphidim

C-----------------------------------------------------------------------
C INFORMATION ABOUT WAVEFUNCTION TYPE
C
C psitype: 0 = MCTDH (default)
C          1 = numerically exact (full grid representation)
C-----------------------------------------------------------------------
      integer psitype

      common /ipsitype/ psitype

C-----------------------------------------------------------------------
C parameters needed for S-MCTDH
C
C citype   : type of selection used
C            1 = static
C cicut    : cutoff for selection
C ciweight(m) : weight for selection
C-----------------------------------------------------------------------
      integer citype
      real*8  cicut,ciweight(maxdim)

      common /icidef/ citype
      common /rcidef/ cicut,ciweight

C-----------------------------------------------------------------------
C POINTERS FOR ARRAYS
C
C dim(m,s)   : number of single particle functions for mode m
C vdim(m,s)  : cumulative number of preceeding degrees of freedom
C ndim(m,s)  : cumulative number of following degrees of freedom
C block(s)   : number of MCTDH A coefficients for state s
C zetf(m,s)  : pointer to start of single particle functions for mode
C              m in psi
C dmat(m,s)  : pointer to start of dichte matrices for mode m
C d3mat(m,s) : pointer to start of dichte3 matrices for mode m
C zpsi(s)    : pointer to wavefunction for state s
C jvdim(m)   : vdim for full jspace
C jdim(m)    : dim for full jspace
C jndim(m)   : ndim for full jspace
C jblock     : total size of full jspace
C maxblock   : Largest size of "block(s)". Also the largest vector 
C              required when operating on psi, e.g. in matrizen2 or 3
C              (see "maxmatwk" zeigprop for the exact size of "maxblock")
C-----------------------------------------------------------------------
      integer dim(maxdim,maxsta),vdim(maxdim,maxsta),ndim(maxdim,maxsta)
      integer zetf(maxdim,maxsta),d3mat(maxdim,maxsta)
      integer dmat(maxdim,maxsta)
      integer zpsi(maxsta),block(maxsta),maxblock
      integer jdim(maxdim),jvdim(maxdim),jndim(maxdim),jblock

      common /psipoint/ dim,vdim,ndim,zetf,dmat,d3mat,zpsi,block,
     +       jdim,jvdim,jndim,jblock,maxblock

C-----------------------------------------------------------------------
C POINTERS FOR ARRAYS USED IN OPWF
C-----------------------------------------------------------------------
      integer d3matsdim,dmatsdim
      integer dims(maxdim,maxsta),zetfs(maxdim,maxsta)
      integer d3mats(maxdim,maxsta),dmats(maxdim,maxsta)
      integer workcdims,workrdims

      common /oppsipoint/ d3matsdim,dmatsdim,dims,zetfs,d3mats,dmats,
     +        workcdims,workrdims

C-----------------------------------------------------------------------
C POINTERS FOR ARRAYS FOR DENSITY OPERATORS
C
C ddim(m,s,t)   : number of single particle functions for mode m
C vddim(m,s,t)  : cumulative number of preceeding degrees of freedom
C nddim(m,s,t)  : cumulative number of following degrees of freedom
C dblock(s,t)   : number of MCTDH A coefficients for state s 
C dzetf(m,s,t)  : pointer to start of single particle functions for mode
C                 m in psi
C ddmat(m,s,t)  : pointer to start of dichte matrices for mode m
C dd3mat(m,s,t) : pointer to start of dichte3 matrices for mode m
C dzpsi(s,t)    : pointer to wavefunction for state s
C-----------------------------------------------------------------------
      integer ddim(maxdim,maxsta,maxsta),vddim(maxdim,maxsta,maxsta)
      integer nddim(maxdim,maxsta,maxsta),dzetf(maxdim,maxsta,maxsta)
      integer dd3mat(maxdim,maxsta,maxsta),ddmat(maxdim,maxsta,maxsta)
      integer dzpsi(maxsta,maxsta),dblock(maxsta,maxsta)

      common /dpsipoint/ ddim,vddim,nddim,
     +        dzetf,ddmat,dd3mat,dzpsi,dblock

C-----------------------------------------------------------------------
C TEMPERATURE, INITIAL STATE AND GROUNDSTATE
C-----------------------------------------------------------------------
      real*8     temp, intrace
      complex*16 inctr
      integer    dlstate, drstate

      common /temper/ temp, intrace
      common /dstates/ dlstate, drstate
      common /incorrtrace/ inctr



