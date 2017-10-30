! Subroutine to provide jobs to different process 
! Operates in the assumptio that number of processes (processors) in not
! large than Morb+Morb*(Morb+1)/2;
! in principle there sould be several possibilities:
!                          a) delta-function - only linear part is
!                          needed to be parallelized;
!                          b)  True separable potential 


          SUBROUTINE  MNGR_Orb_Part(Nproc)
          USE SHARED_DIMS
          USE W_INTERPARTICLE
          USE DVR_ALL
          USE Parallel_Orb
 
          IMPLICIT NONE
c==============================================================================
          INTEGER :: I,J,K,icntr,Nproc,M,nn,FromN,TillN,PRL_TYPE
          INTEGER :: NJOBS,NJOBS1,NJOBS2,n_per_proc,nrest,nfnl
          INTEGER,  DIMENSION(10000) ::  nwork
c==============================================================================
!       I USE "Parallel" arrays: Proc_Job(Nproc) 
!                                Proc_Iorb_Starts -> summation over Iorb in Func_Str 
!                                Proc_Iorb_Finish
!                                Proc_Jdim_Starts -> summation over Dimensionalities in Func_St
!                                Proc_Jdim_Finish
!                                Proc_From
!                                Proc_Till -> part which is computed in   Get_Op_PSI(Vec,WOp_X,NDX,1,From,Till)
c==============================================================================
c============================== If (.FALSE.) W(x-x')=delta(x-x')*xlambda0 is used
c============================== or if (.TRUE.) - true interparticle potential 
c============================== W(r-r')*xlambda0 is used
c============================== defined by Wxx(I,J)=W(r,r') in DVR grid
c============================== or by vector W2xx(J)=W(ort(1)-ort(J))
c============================== or for separable case: WOP_X*WOP_Y*WOP_Z
            M=Morb
            NJOBS1= DIM_MCTDHB*M
            NJOBS1= M
            NJOBS2= DIM_MCTDHB*(M*(M+1))/2
       WRITE(6,*)"ORB:Parallel Menager for Orbital Part ",Nproc,Wxx_TYPE

c=========================== Delta case ========================================
!           PRL_TYPE=0  ! WSL is computed on Master's node and KIN+V  are MPI      -DISBALANCED 
!           PRL_TYPE=1  ! WSL and KIN+V are BOTH MPI parallelized               -   BALANCED 
         PRL_TYPE=0 
         IF(NDX*NDY*NDZ.ge.1024) PRL_TYPE=1
         IF(Wxx_TYPE.ne.0) PRL_TYPE=1
         PRLTYPE: SELECT CASE (PRL_TYPE)
c=========================== Delta case ========================================
         CASE (0)
       write(6,*)"ORB: Parallel menager for Orbitals is DISBALANCED:" 
       write(6,*)"ORB: DISBALANCED: NonLinear on Master, Kinetics - MPI"
          Proc_Iorb_Starts=0
          Proc_Iorb_Finish=0
          nwork=0
      IF(Morb.eq.1) THEN
      Proc_Iorb_Starts(1)=1
      Proc_Iorb_Finish(1)=1
      ELSE
         NJOBS=NJOBS1+1 
        IF(Nproc.ge.NJOBS) THEN
        icntr=2
        Do I=1,M
        Proc_Iorb_Starts(icntr)=I
        Proc_Iorb_Finish(icntr)=I
        icntr=icntr+1
        ENDDO
        ENDIF

         NJOBS=NJOBS1+1 
         IF(Nproc.lt.NJOBS) THEN
        write(6,*)"ORB: NPROC",nproc
        write(6,*)"ORB: NJOBS",njobs
        write(6,*)"ORB: NORB",morb
         n_per_proc=INT(NJOBS/Nproc)
         nrest=MOD(NJOBS,Nproc)
         nwork=INT(NJOBS/Nproc)
         DO I=1,nrest
         nwork(I+1)=nwork(I+1)+1
         ENDDO
            nfnl=-1
         DO I=1,Nproc
          Proc_Iorb_Starts(I)=nfnl+1
          Proc_Iorb_Finish(I)=nwork(I)+nfnl
          nfnl= Proc_Iorb_Finish(I)
         ENDDO
          Proc_Iorb_Starts(1)=1
       ENDIF
         ENDIF        


         nn=NDX*NDY*NDZ
        Orb_Block=0
         Do I=1,Nproc
         IF(Proc_Iorb_Finish(I).ne.0)
     .   Orb_Block(I)=nn*(Proc_Iorb_Finish(I)-Proc_Iorb_Starts(I)+1)
         EnDDo
         Orb_Displ(1)=0
         Do I=2,Nproc
         IF(Proc_Iorb_Finish(I).ne.0)
     .    Orb_Displ(I)=Orb_Displ(I-1)+Orb_Block(I-1)
         EnDDo
         Do I=1,Nproc
         IF(Proc_Iorb_Finish(I).eq.0) Proc_Iorb_Starts(I)=1
         EnDDo
         
         Do icntr=1,Nproc
       write(6,*)icntr,"Fr",Proc_Iorb_Starts(icntr), " Tl",
     . Proc_Iorb_Finish(icntr),"BLC",Orb_Block(icntr)," DS",
     . Orb_Displ(icntr)
         EnDDo
c============================ WLS - is computed on master node
         Proc_WSL_Starts=0
         Proc_WSL_Finish=-1
            NJOBS= (M*(M+1))/2
        write(6,*)"ORB: WSL: NPROC",nproc
        write(6,*)"ORB: WSL: NJOBS",njobs
        write(6,*)"ORB: WSL: NORB",morb
          Proc_WSL_Starts(1)=1
          Proc_WSL_Finish(1)=njobs
         Do icntr=1,Nproc
       write(6,*)icntr,"ORB: WSL From",Proc_WSL_Starts(icntr), " Tl",
     . Proc_WSL_Finish(icntr),"out of ",M*(M+1)/2
         EnDDo
c======================================== ALL are MPI -Paralellized 
         CASE (1)
       write(6,*)"ORB: Parallel menager for Orbitals is BALANCED:" 
       write(6,*)"ORB: BALANCED: Non-Linear  and Kinetics are MPI"
          Proc_Iorb_Starts= 0
          Proc_Iorb_Finish=-1
c         WSL=xlambda0*Conjg(psiS)*PsiL/weight/weight
         NJOBS=NJOBS1 ! ONLY Kinetic Part is parallelized 2-body is not expensive
        write(6,*)"ORB NJOBS ",NJOBS," on ",Nproc," Proc's"
         Proc_Job=1
         IF(Nproc.ge.NJOBS) THEN
         icntr=1
         Do I=1,M
          Proc_Iorb_Starts(icntr)=I
          Proc_Iorb_Finish(icntr)=I
          icntr=icntr+1
         ENDDO
c        write(6,*)" NJOBS",NJOBS,"?equiv?",icntr," on ",Nproc," Proc's"
         ENDIF

         IF(Nproc.lt.NJOBS) THEN
         n_per_proc=NJOBS/Nproc
c        write(6,*)" N per Proc",n_per_proc
         Do I=1,Nproc
         n_per_proc=NJOBS/Nproc
         FromN=(I-1)*n_per_proc+1
         Proc_Iorb_Starts(I)=FromN
      IF(I.gt.Nproc-(NJOBS-n_per_proc*Nproc)) THEN
       FromN=FromN +(I-1-(Nproc-(NJOBS-n_per_proc*Nproc)))
!      FromN=FromN+(MYID-(NPROC-(MaxTrm2b-Nterms*NPROC))
       Proc_Iorb_Starts(I)=FromN
       n_per_proc=NJOBS/Nproc+1
      ENDIF
      TillN=FromN+n_per_proc-1
      IF(TillN.ge.NJOBS) TillN=NJOBS
      Proc_Iorb_Finish(I)=TillN
c       write(6,*)I,"From",Proc_Iorb_Starts(I), " Till",
c     . Proc_Iorb_Finish(I)
         ENDDO
         ENDIF
c============================= Hand MAde 
         nn=NDX*NDY*NDZ
         Do I=1,Nproc
         Orb_Block(I)=nn*(Proc_Iorb_Finish(I)-Proc_Iorb_Starts(I)+1)
         EnDDo
         Orb_Displ(1)=0
         Do I=2,Nproc
         Orb_Displ(I)=Orb_Displ(I-1)+Orb_Block(I-1)
         EnDDo

         Do icntr=1,Nproc
       write(6,*)"ORB: ",icntr,"Fr",Proc_Iorb_Starts(icntr), " Tl",
     . Proc_Iorb_Finish(icntr),"BLC",Orb_Block(icntr)," DS",
     . Orb_Displ(icntr)
         EnDDo
c         write(6,*)"Orb_Block",Orb_Block
c         write(6,*)"Orb_Displ",Orb_Displ
c          pause
c============================ WLS manager p
         Proc_WSL_Starts=0
         Proc_WSL_Finish=0
            NJOBS= (M*(M+1))/2
       IF(Nproc.le.NJOBS) THEN
        write(6,*)"ORB WSL: NPROC",nproc
        write(6,*)"ORB WSL: NJOBS",njobs
        write(6,*)"ORB WSL: NORB",morb
         n_per_proc=INT(NJOBS/Nproc)
         nrest=MOD(NJOBS,Nproc)
         nwork=INT(NJOBS/Nproc)
         DO I=1,nrest
         nwork(I+1)=nwork(I+1)+1
         ENDDO
            nfnl=0
         DO I=1,Nproc
          Proc_WSL_Starts(I)=nfnl+1
          Proc_WSL_Finish(I)=nwork(I)+nfnl
          nfnl= Proc_WSL_Finish(I)
         ENDDO
       ENDIF
         Do icntr=1,Nproc
       write(6,*)"ORB:",icntr,"WSL From",Proc_WSL_Starts(icntr), " Tl",
     . Proc_WSL_Finish(icntr),"out of ",M*(M+1)/2
         EnDDo
c============================General Cases ======================================
c=================== Case Separable W WOp_X
        CASE (10,11)
       WRITE(6,*)" Is not yet implemented !!!!Wxx_TYPE",Wxx_TYPE
       stop 
c=================== Case Unseparable W(x-x') - only one coloumn is needed
c       CASE (2)
c       WRITE(6,*)" Is not yet implemented !!!!Wxx_TYPE",Wxx_TYPE
c       stop 
c       RETURN
c=================== Case Unseparable W(x,x') - all elements Wxx are needed 
       CASE (4)
       WRITE(6,*)" Is not yet implemented !!!!Wxx_TYPE",Wxx_TYPE
       stop 
       RETURN
c===============================================================================
       CASE DEFAULT
       WRITE(6,*)" In Menager_Orb_Part case Default???? Why?,Wxx_TYPE"
     . ,Wxx_TYPE
       WRITE(6,*)" STOP ERRRR "
       stop 
       END SELECT PRLTYPE
       END  SUBROUTINE MNGR_Orb_Part
