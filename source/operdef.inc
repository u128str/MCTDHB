C***********************************************************************
C                                                                     
C              MCTDH include file  operdef.inc 
C                                                                     
C***********************************************************************

      integer maxhtm,maxdef,maxfac,maxmuld
      parameter (maxhtm=1200)
      parameter (maxfac=1200)
      parameter (maxdef=100)
      parameter (maxmuld=10)

C-----------------------------------------------------------------------
C parameters defining operator modes (may be different from system)
C-----------------------------------------------------------------------
       integer nopmode,nopdof

       common /opmodes/ nopmode,nopdof

C-----------------------------------------------------------------------
C POINTERS FOR ARRAYS
C
C-----------------------------------------------------------------------
      integer zhop(maxhtm)

      common /operpoint/ zhop

C-----------------------------------------------------------------------
C the hterm and associated arrays
C-----------------------------------------------------------------------
      integer hterm(maxfac,maxhtm),htmdof(maxhtm),htmmode(maxhtm),
     +        htmmuld(maxhtm),htmfac(maxhtm),htmform(maxhtm),
     +        htmsym(maxhtm)

      common /htmdef/ hterm,htmdof,htmmode,htmmuld,htmfac,htmform,
     +        htmsym

C-----------------------------------------------------------------------
C SOME MORE PARAMETERS NEEDED FOR HAMILTONIAN
C
C nopts(n)  : number of parameters storred in pesopts(n) for nth surface
C-----------------------------------------------------------------------
      integer maxhpar,maxhop,maxpes
      parameter (maxhpar=5)
      parameter (maxhop=1000)
      parameter (maxpes=20)

      integer hopfile(maxhop),hopdim,hopipar(maxhpar,maxhop),
     +        hopilab(maxhop),npesopt,nhoppar,zpesopts(maxhpar,maxpes),
     +        nopts(maxpes)
      real*8 hoprpar(maxhpar,maxhop)
      real*8 hoppar(maxhpar,maxhop)
      character*(c2) hoplab(maxhop)
      character*(c5) pesopts(maxpes)

      common /ihoplab/ hopfile,nhoppar,hopdim,hopilab,hopipar,
     +        npesopt,zpesopts,nopts
      common /ahoplab/ hoplab,pesopts
      common /rhoplab/ hoprpar,hoppar

C-----------------------------------------------------------------------
C parameters used to cut high/low values of correlated PES surface
C-----------------------------------------------------------------------
      logical lcutpes
      real*8  pesmin,pesmax

      common /lcutpescb/ lcutpes
      common /rcutpes/ pesmin,pesmax

C-----------------------------------------------------------------------
C properties of the PES surface
C
C iqconst(f) : if 1, and dof is not in specified system, coordinate
C              is set to value in rqconst(f)
C              if 0, and dof is not in specified system, coordinate
C              is set to 0.0d0
C              if dof is in specified system, this flag is ignored.
C-----------------------------------------------------------------------
      integer iqconst(maxdim)
      real*8 vmin,vmax,rqconst(maxdim)

      common /ipesprop/ iqconst
      common /rpesprop/ vmin,vmax,rqconst

C-----------------------------------------------------------------------
C arrays storing information for multi-D operators
C
C mulddim: no. of multi-D sets of coordinates
C
C-----------------------------------------------------------------------
      integer mulddim,mulddof(maxdim,maxmuld),muldsize(maxmuld),
     +        nmulddof(maxmuld) 

      common /mulddef/ mulddim,mulddof,muldsize,nmulddof

C-----------------------------------------------------------------------
C PARAMETERS NEEDED FOR HAMILTONIAN
C
C hopsdim   : length of array containing hops
C  lcap   : calculation employing CAPs
C  lcdvr  : any multidimensional Hamiltonian terms use the CDVR method
C  ltdvr  : any multidimensional Hamiltonian terms use the TDVR method
C  lanalpes : an analytic PES will be generated (e.g. for CDVR)
C-----------------------------------------------------------------------
      integer      hopsdim,ktot
      logical      lcap,lcdvr,ltdvr,lanalpes

      common /iopinitg/ hopsdim,ktot
      common /lopinitg/ lcap,lcdvr,ltdvr,lanalpes

C-----------------------------------------------------------------------
C SOME MORE PARAMETERS NEEDED FOR HAMILTONIAN
C
C kzahl     : number of correlated Hamiltonian expansion coefficients
C koeff(k)  : expansion coefficients of correlated Hamiltonian
C ki(k),kf(k): defines states for Hamiltonian term <kf|H|ki>
C complete(m) : if true, spf basis is "complete" and (1-P)H phi =0
C diag(m,k)   : array of unitary spf hamiltonian matrices
C hsym(m,k)   : array of operator symmetries
C hlsym(m,k)  : array of hloch symmetries
C diagh(m)    : auxiliary diag array
C kinop(f,k)  : logical set to true if a correlated Hamiltonian term is
C               a kinetic (i.e. non-diagonal) operator
C laddmode(m,nh)  : logical set to true if the mode m is to be included
C                   into the operator nh. (See keyword 'addmode').
C khzahl(m,s)  : number of uncorrelated coefficients for dimension m
C                and state s
C hamlab(nham) : string with name of nham-th Hamiltonian.
C hamilton(f,k) : integer with no. of kth Hamiltonian term for dof f.
C zhun(m,s) : pointer to start of uncorrelated terms in hamilton array
C ldissop: indicates if dissipative operators are involved
C ldissin: indicates if input wishes dissipative treatment
C ldproj:  indicates if the equations of motion of the parameters
C          shall be modified
C-----------------------------------------------------------------------
      integer maxham
      parameter (maxham=10)

      integer kzahl(maxham),nmulpot(maxham),khtot(maxham)
      integer ki(maxkoe),kf(maxkoe)
      complex*16 koeff(maxkoe)
      logical complete(maxdim),diag(maxdim,maxkoe),diagh(maxdim),
     +        kinop(maxdim,maxkoe),laddmode(maxdim,maxham),
     +        ldissop,ldissin,ldproj

      integer khzahl(maxdim,maxsta,maxham),
     +        zhun(maxdim,maxsta,maxham),hamilton(maxdim,maxkoe),
     +        htmtype(maxhtm),htmdim,zham(maxham),hamdim,
     +	      hamnum(maxkoe),hsym(maxdim,maxkoe),hlsym(maxdim,maxkoe)
      character*(c2) hamlab(maxham)

      common /iop/ kzahl,khzahl,hamilton,htmtype,zhun,khtot,
     +             nmulpot,ki,kf,htmdim,zham,hamdim,hamnum,
     +             hsym,hlsym
      common /cop/ koeff
      common /lop/ complete,diag,diagh,kinop,laddmode,ldissop,
     +             ldissin,ldproj
      common /aop/ hamlab

C-----------------------------------------------------------------------
C flag to turn off use of diag flag and extraction of huncorr terms.
C-----------------------------------------------------------------------
       logical nodiag(maxham)

       common /hamflags/ nodiag

C-----------------------------------------------------------------------
C logical variables for FT in potfit
C-----------------------------------------------------------------------
       logical lft,lftsym
 
       common/plftt/ lft,lftsym

C-----------------------------------------------------------------------
C
C Flags for operator action
C
C kfirst, kact, ktwin
C
C ktwin points to the twin operator in a dissipative part. 
C kact: 1: for wave-functions
C       2: pure-left
C       3: pure-right
C       4: left-right	
C
C-----------------------------------------------------------------------
      integer   kfirst(maxdim,maxkoe),kact(maxkoe),ktwin(maxkoe)

      common /iop1/ kfirst,kact,ktwin

C-----------------------------------------------------------------------
C Variables needed for autocap option

C countpsi : number of output of psi in psi file
C autocap(k) : if > 0, k the term includes an autocap on dof f=autocap(k).
C              The cap is then activated when it is populated.
C indicap(f) : number of grid point where the CAP starts
C gpopcap(f) : grid population at the point where the CAP starts
C-----------------------------------------------------------------------
      integer countpsi,indicap(maxdim),autocap(maxkoe)
      real*8  gpopcap(maxdim)

      common /iautcap/ countpsi,indicap,autocap
      common /rautcap/ gpopcap

C-----------------------------------------------------------------------
C Information about natpot terms:
C
C modc    - the contracted mode of the natpots
C nterm   - the number of the first natpot term in the correlated terms
C           of hamiltonian
C natcdof - the dof for which natpot terms are pre-multiplied (in "fast"
C           algorithm)
C-----------------------------------------------------------------------
      integer vpdim(maxdim),pdim(maxdim),npdim(maxdim),modc,nterm,
     +        natcdof

      common /natdata/ vpdim,pdim,npdim,modc,nterm,natcdof

C-----------------------------------------------------------------------
C Variables for fourier transformation of some dofs in potfit program
C
C pdofft  : dof, which should be fourier transformed
C npdofft : number of grid points used in fourier transformation
C----------------------------------------------------------------------- 
      integer pdofft,npdofft

      common /pft/ pdofft,npdofft

C-----------------------------------------------------------------------
C SOME MORE PARAMETERS NEEDED FOR EXPECTATION VALUES
C
C nexpect     : number of operators to derive the expectation-value
C expectlab() : labels for the expectation operators
C-----------------------------------------------------------------------
      integer nexpect
      complex*16 expvalue(maxham)
      character*(c2) expectlab(maxham)

      common /lexpectation/ expectlab
      common /iexpectation/ nexpect
      common /cexpectation/ expvalue



