C***********************************************************************
C
C                         MCTDH include file GLOBAL 
C                                                                     
C  MCTDH include file defining various parameters and variables       
C  global to the program.     
C                                                                    
C   V6.0 MB
C                                                                    
C***********************************************************************

C-----------------------------------------------------------------------
C PARAMETERS FOR CHARACTER STRING LENGTHS
C
C c1, c2, c3, c4, c5: various lengths
C DO NOT CHANGE c3. Must match keyword length in rdinpf library routine
C-----------------------------------------------------------------------
      integer c1, c2, c3, c4, c5

      parameter (c1=16)
      parameter (c2=32)
      parameter (c3=60)
      parameter (c4=80)
      parameter (c5=160)

C-----------------------------------------------------------------------
C progver: code version number, defines output formats etc.
C tiny is added to all these numbers to ensure that version number is
C .ge. to that desired!
C Tiny is also used as a tolerance when real or complex numbers are
C compared.
C filever(i): file version number for file on channel i.
C             Defines output formats etc.
C  NB: filever must be dimensioned larger than largest channel no. used.
C-----------------------------------------------------------------------
      real*8 tiny
      parameter (tiny=1.0d-8)
      real*8 progver,filever(100)
      common /version/ progver,filever

C-----------------------------------------------------------------------
C COMMON FOR INPUT / OUTPUT CHANNELS
C initial channels are specified in main program, then transferred by
C this common block to allow modifications to be made in the program
C e.g. the output file is on channel 6 if printed to the screen, or
C channel 13 if written to file
C
C iin:   .inp file 
C ilog:  log file 
C iolog: op.log file
C irst:  restart file
C itim:  timing file
C ioper: operator file
C idvr:  dvr file
C iout:  output file
C ichk:  check file
C irho:  gridpop file
C ipden: pdensity file (one particle density)
C ictr:  ctrace file
C icroco: cross file (cross-correlation function)
C iaut:  auto file
C iaute: autoe file 
C iaut1: auto1 file (contains c_1 = <psi(0)| H |psi(t)>)
C iaut2: auto2 file (contains c_2 = <psi(0)|H^2|psi(t)>)
C istp:  intstep file
C istop: stop file
C iupd:  update file
C ispeed: speed file
C irlx:  rlx_info file
C idat:  daten file
C iop:   operator text file
C ipar:  Hamiltonian parameter file
C itmp:  scratch file used in building of operator 
C iadwkb: diabatically corrected initial energy distribution
C inpot   : 'natpot' file
C ippsi   : 'ppsi' file
C iit     : 'iteration' file
C ipw     : 'prodwei' file
C igit    : 'pgitter' file
C icr1:   seems to be unused
C iwtt:  <Psi(t)|W_gamma|Psi(t)> data
C iflx:  reactive flux data
C iexpect: file of expectation values
C iexim: implicit psi file data
C irp:   reaction probabilities
C ipsi1: second psi file
C iaus:  general output file used in analyse programs
C ieigval: file of eigenvalues and intensities
C ilancz:  file of Lanczos vectors
C ieigvec: file of eigenvectors of tridiagonal Lanczos matrix
C ivpot   : 'vpot' file
C iveig   : veigen file
C iorb    : orben file, orbital energies

C-----------------------------------------------------------------------
      integer iin,ilog,irst,itim,ioper,idvr
      integer iout,ipsi,ichk,irho,iaut,iaute,istp,istop,
     +        iupd,ispeed,irlx,iaut1,iaut2,ictr,ipden
      integer idat,iop,ipar,itmp,iadwkb,iolog
      integer inpot,iit,ipw,ippsi,igit,iaus
      integer icr1,iwtt,iflx,iexpect,iexim,irp,ipsi1,ieigval,ilancz
      integer ieigvec,ivpot,iveig,iorb,icroco

      common /channels/ iin,ilog,iolog,irst,itim,ioper,idvr,
     +                  ipsi,ichk,iout,irho,iaut,iaute,istp,
     +                  istop,iupd,ispeed,irlx,iaut1,iaut2,ictr,
     +                  ipden,idat,iop,ipar,itmp,iadwkb,
     +                  inpot,iit,ipw,ippsi,igit,iaus,
     +                  icr1,iwtt,iflx,iexpect,iexim,irp,ipsi1,
     +                  ieigval,ilancz,ieigvec,ivpot,iveig,iorb,icroco


C-----------------------------------------------------------------------
C CHARACTER STRINGS TO PASS MESSAGES TO ERROR ROUTINE
C-----------------------------------------------------------------------
      character*(c5) routine, message, inptit

      common /errorhandling/ routine, message, inptit


C-----------------------------------------------------------------------
C chkerr: array used for checking system definition
C
C chkerr(1)   : Error reading file
C chkerr(2)   : Dofs differ from system
C chkerr(3)   : DVR bases differ from system
C chkerr(4)   : Mode combinations differ from system
C chkerr(5)   : Electronic basis type differs from system
C chkerr(6)   : Spf basis differs from system 
C chkerr(7)   : Dofs not subset of system dofs
C chkerr(8)   : ? not for this system
C chkerr(9)   : Psi type differs from system
C chkerr(10)  : Density matrix type differs from system
C 
C-----------------------------------------------------------------------
      integer nerr
      parameter (nerr=15)

      integer chkerr(nerr)

      common /rwerror/ chkerr


C-----------------------------------------------------------------------
C DIRECTORY DEFINITION FOR INPUT/OUTPUT FILES
C
C name      :  name of directory in which the output files are stored
C inputfile : path of the inputfile (~.inp)
C options   : string containing the list of options
C laenge    : length of string name
C inplaenge : length of string inputfile
C optlaenge : length of string options
C rname     : directory from which to take restart file
C rlaenge   : length of string rname
C oname     : directory from which to take operator file
C olaenge   : length of string oname
C dname     : directory from which to take dvr file
C dlaenge   : length of string dname
C vname     : name of a vpot file to read
C vlength   : length of string vname
C operfile  : path and name of oper file
C opflaenge : length of string operfile
C crossname : directory from which to take restart file used for
C             calculating the cross-correlation-function
C crosslng  : length of string crossname
C-----------------------------------------------------------------------
      character*(c5)  name,inputfile,rname,dname,oname,vname,operfile
      character*(c5)  options,runtitle,crossname
      integer         laenge,inplaenge,optlaenge,rlaenge,dlaenge,
     +                olaenge,vlength,opflaenge,crosslng

      common /names/ name,inputfile,options,runtitle,rname,dname,oname,
     +               vname,operfile,crossname
      common /lengths/ laenge,inplaenge,optlaenge,rlaenge,dlaenge,
     +                 olaenge,vlength,opflaenge,crosslng

C-----------------------------------------------------------------------
C PARAMETERS FOR CALCULATIONS
C
C pi : pi is defined by 4d0*atan(1d0) in subroutine default.
C fs : time conversion from a.u. to fs (x[a.u.] = fs*x[fs]).
C ev : energy conversion from a.u. to ev (x[a.u.] = x[eV]/ev).
C   Note: If a.u. are output, fs=ev=1
C one: the variable 1
C oner: the variable 1.0d0
C True,False: the logicals .true. and .false., respectively
C-----------------------------------------------------------------------
      integer   one,zeroi,onei,twoi,threei,fouri,fivei,sixi,seveni
     +          ,eighti,ninei,teni
      logical   true, false

      real*8    oner
      parameter (oner = 1.0d0)
      parameter (one=1,zeroi=0,onei=1,twoi=2,threei=3,fouri=4
     +  ,fivei=5,sixi=6,seveni=7,eighti=8,ninei=9,teni=10)
      parameter (true = .true., false = .false.)

      real*8    pi, fs, ev
      common /constants/ pi, fs, ev

C-----------------------------------------------------------------------
C PARAMETERS read from input 
C
C npar      : number of parameters read
C apar(n)   : string with name of nth parameter
C rpar(n)   : (real) value of nth parameter
C-----------------------------------------------------------------------
      integer maxpar
      parameter (maxpar=350)

      integer npar
      real*8 rpar(maxpar)
      character*(c1) apar(maxpar)

      common /iparam/ npar
      common /rparam/ rpar
      common /aparam/ apar

C-----------------------------------------------------------------------
C  POINTERS for the start positions of
C  the arrays stored in the memory arrays mc, mr, ms, mi and ml.
C
C see zeiger routines for definitions
C
C-----------------------------------------------------------------------

C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C !! If you add any further constant, change the mccount,..,mfcount   !!
C !! values in subroutine checkmemprop !!                             !!
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer mcdim,mchpsi,mchunpsi,mcintpsi,mchin,mcrueck,mcexphin,
     +        mcexprueck,
     +        mchteil,mchloch,mcdicht1,mcdicht2,mcdicht4,mcdicht4o,
     +        mchnatur,mcwmat,mcipbt,mcipft,
     +        mcworkc,mccdvr2,mccdvr4,mcmfsum,mcoldpsi,
     +        mcpredpsi,mcautopsi,mcautodtpsi,mchc,mcspf1d,mcphi,
     +        mcpsiref,mcccdover
      integer mrdim,mrdicht3,mrort,mrdvrmat,mrtrafo,mrfftp,
     +        mrworkr,mrgpop1,mrgpop2,mrcdvr3,mrhunew,mrhunev,
     +        mrdicht3o,mrkinsph,mrhops,mrhr,mreigval,mreigerr,
     +        mreigint
      integer midim,miworki,mijsph,mimsph,mifftfak,mihi,mijindx,
     +        miccjindx
      integer mldim,mlworkl,mlhl
      integer msdim,msworks
      integer mfdim,mfworkf

      common /mcpoint/ mcdim,mchpsi,mchunpsi,mcintpsi,mchin,mcrueck,
     +        mchteil,mchloch,mcdicht1,mcdicht2,mcdicht4,mchnatur,
     +        mcdicht4o,mcwmat,mcipbt,mcipft,mcworkc,
     +        mccdvr2,mccdvr4,mcmfsum,mcoldpsi,mcpredpsi,mcautopsi,
     +        mcautodtpsi,mchc,mcspf1d,mcphi,mcexphin,mcexprueck,
     +        mcpsiref,mcccdover 

      common /mrpoint/ mrdim,mrdicht3,mrworkr,mrort,mrdvrmat,
     +        mrtrafo,mrfftp,mrdicht3o,mrgpop1,mrgpop2,mrcdvr3,
     +        mrhunew,mrhunev,mrkinsph,mrhops,mrhr,mreigval,mreigerr,
     +        mreigint

      common /mipoint/ midim,miworki,mijsph,mimsph,mifftfak,mihi,
     +        mijindx,miccjindx

      common /mlpoint/ mldim,mlworkl,mlhl

      common /mspoint/ msdim,msworks

      common /mfpoint/ mfdim,mfworkf

C-----------------------------------------------------------------------
C WORK ARRAY SIZES
C
C workcdim : Length of the complex*16 scratch array
C workrdim : Size of real*8 scratch array 
C workldim : Size of logical scratch array
C workidim : Size of integer scratch array
C worksdim : Size of complex*8 scratch array 
C workfdim : Size of real*4 scratch array 
C-----------------------------------------------------------------------
      integer workcdim,workrdim,workldim,workidim,worksdim,workfdim

      common /arraydim/ workcdim,workrdim,workldim,workidim,
     +        worksdim,workfdim

C-----------------------------------------------------------------------
C FLAG FOR WAVEFUNCTION / DENSITY MATRIX MODE
C
C dentype: 0 wavefunction
C          1 density matrix type I
C          2 density matrix type II
C-----------------------------------------------------------------------
      integer dentype

      common /method/ dentype

C-----------------------------------------------------------------------
C LOGICALS TO CONTROL INPUT / OUTPUT OPTIONS
C
C If set to .true.,
C     neu      : calculation is a new run
C     cont     : calculation is a continuation run
C     restart  : calculation is a restart run
C     test     : runs up to propagation to test memory etc.
C     test0    : short test-run (via option -s)
C     lovwr    : files in directory are overwritten if neu 
C     logisopen: the log file is open
C     lwarn    : is set to 'true' in subroutine 'warnmsg'
C     liopt    : is set to 'true' when input is read from name/input. 
C     ltime    : controls timing output
C     lnormstop: if true, the keyword 'normstop' in run-section is used
C     lfast    : if true, the fast algorithm for H(npot)*A is used
C-----------------------------------------------------------------------
      logical neu,cont,restart,test,lovwr,lwarn,lnpot,ltime
      logical test0,logisopen,liopt,lfast
      character*8 sstop,scpu
      real*8  nrmstop
      logical lnormstop, lstop

      common /blstop/ lnormstop, lstop
      common /stopnorm/ nrmstop
      common /stopstring/ sstop,scpu
      common /strtopt/ neu,cont,restart,test,lovwr,lnpot,test0,lfast
      common /ioopt/   logisopen,ltime,lwarn,liopt

C-----------------------------------------------------------------------
C 
C FLAGS TO SPECIFY THE RUN-TYPE AND RW-FILE GENERATION
C
C  If set to .true., 
C     lrundvr  : A DVR file will be generated
C     lrunoper : An OPER file will be generated
C     lruninwf : A RESTART file will be generated
C     lrunpes  : A PES file will be generated
C     lrunprop : A propagation will be made
C     lrundiag : if true, the Hamiltonian is diagonalised using Lanczos
C     lrddvr   : A DVR file will be read
C     lrdoper  : An OPER file will be read
C     lrdinwf  : A RESTART file will be read
C     ldeldvr  : DVR file is deleted after use
C     ldelop   : OPER file is deleted after use
C     lexact   : numerically exact calculation
C-----------------------------------------------------------------------
      logical lrundvr,lrunoper,lruninwf,lrunprop,lrddvr,lrdoper,
     +        lrdinwf,lrunpes,lrundiag,lexact,ldelop,ldeldvr

      common /runopts/ lrundvr,lrunoper,lruninwf,lrunprop,lrddvr,
     +                 lrdoper,lrdinwf,lexact,lrundiag,lrunpes,
     +                 ldelop,ldeldvr

C-----------------------------------------------------------------------
C ACCURACY VARIABLES
C
C macheps: machine precision, i. e. smallest positive number for which
C          1.0+macheps > 1.0  (macheps=2.2d-16)
C-----------------------------------------------------------------------
       real*8 macheps
       common /macheps/ macheps

C-----------------------------------------------------------------------
C VARIABLES USED TO READ INPUT FILE
C (defined here so that same size arrays are used in different routines)
C-----------------------------------------------------------------------
      integer maxkey
      parameter (maxkey=30)

      character*(c3) keyword(maxkey),keyorig(maxkey)
      integer      lc(maxkey)

      common /akey/ keyword,keyorig
      common /ikey/ lc


C-----------------------------------------------------------------------
C array needed to specify which parts of read-write files are to be read
C-----------------------------------------------------------------------
      integer nrwdata
      parameter (nrwdata=15)

      logical dvrdata(nrwdata)

      common /rwlogicals/ dvrdata


C-----------------------------------------------------------------------
C integers for outputcontrol of pdensity: 
C npd1...npd4 -- if set -- specify the modes for which the pdensity is 
C written.
C-----------------------------------------------------------------------
      integer npd1,npd2,npd3,npd4
      common /npd/ npd1,npd2,npd3,npd4











