         SUBROUTINE Integrator_ABM &
        (PSI,OPSI,AuxPsi,AbsTime,xIntPeriod,NewRestart,&
         ABMError,abmvar)
         USE SHARED_DIMS
         IMPLICIT NONE
!c======================================================
        external abmvar
        COMPLEX*16, DIMENSION(NDX*NDY*NDZ,Morb) :: Psi,OPsi
        COMPLEX*16, DIMENSION(NDX*NDY*NDZ*Morb,18) ::AuxPsi

        REAL*8 :: xIntPeriod, AbsTime,xInitStep,TolError,ABMError
        INTEGER ::IntOrder,NSteps,NRepeatedSteps,iErrorCode,ierr,PsiDim
        LOGICAL :: NewRestart
        EXTERNAL FUNCSTR,AbsBSError,PolyExtrapol,AbsABMError,dummypost

        REAL*8 :: RelTime,xIntPeriod_ORG
        REAL*8 :: ActualLargeStep, NextLargeStep
       INTEGER :: SmallSteps

!c==================================== RUNGE KUTTA STUFF
        REAL*8 :: hopt ! last optimal step

         xIntPeriod_ORG=xIntPeriod
         SmallSteps=0


         TolError=ABMError
         xInitStep=xIntPeriod
         IntOrder=Time_intgr_order
         NSteps=0
         NRepeatedSteps=10
         iErrorCode=0
         PSiDim=NDX*NDY*NDZ*Morb
         RelTime=AbsTime


      IF (trim(Time_intgr).eq.'OMPABM') then
        if(Time_intgr_order.gt.7) then 
     write(6,*)"ADAMS-BASHFORTH-MOULTON IS IMPLEMENTED TILL 8 ORDER"
     write(6,*)"CHANGE Time_intgr_order<=7"
        stop "Integrator ABM"
        endif
      call  OMPABM(Psi,OPsi,PsiDim,xIntPeriod,RelTime,IntOrder, &
                  xInitStep,TolError,NewRestart,NSteps,NRepeatedSteps,&
                  iErrorCode,AuxPsi,FUNCSTR,AbsABMError)
      IF(iErrorCode.ne.0) write(6,*)"IN OMPABM !!!!!!! I -ERROR",iErrorCode
         Intgr_total_steps=Intgr_total_steps+NSteps
         IF(print_level.ge.3) then
       write(6,'(a30,i2,a22,i6,i16)')"OMP ADAMS-BASHFORTH-MOULTON of" &
       ,Time_intgr_order," order takes SmallSteps", Nsteps, Intgr_total_steps
         ENDIF !PRINT_LEVEL
       endIF
!c================================ END OMPABM ==============================
!c================================ END OMPABM ==============================
!c================================ END OMPABM ==============================

      IF (trim(Time_intgr).eq.'ABM') then
        if(Time_intgr_order.gt.7) then 
     write(6,*)"ADAMS-BASHFORTH-MOULTON IS IMPLEMENTED TILL 8 ORDER"
     write(6,*)"CHANGE Time_intgr_order<=7"
        stop "Integrator ABM"
        endif
      call  ABM(Psi,OPsi,PsiDim,xIntPeriod,RelTime,IntOrder, &
                  xInitStep,TolError,NewRestart,NSteps,NRepeatedSteps,&
                  iErrorCode,AuxPsi,FUNCSTR,AbsABMError)
      IF(iErrorCode.ne.0) write(6,*)"IN ABM !!!!!!! I -ERROR",iErrorCode
         Intgr_total_steps=Intgr_total_steps+NSteps
         IF(print_level.ge.3) then
       write(6,'(a30,i3,a22,i6,i16)')"ADAMS-BASHFORTH-MOULTON of" &
      ,Time_intgr_order," order takes SmallSteps", Nsteps,Intgr_total_steps
         ENDIF !PRINT_LEVEL
       endIF
!c================================ ABM END ==============================
!c================================ ABM END ==============================
!c================================ ABM END ==============================


!c================================ RUNGE KUTTA==============================
!c================================ RUNGE KUTTA==============================
!c================================ RUNGE KUTTA==============================
       IF (trim(Time_intgr).eq.'RK') then
        if((Time_intgr_order.ne.5).and.&
           (Time_intgr_order.ne.8)) then 
           write(6,*)"RUNGE-KUTTA works with order 5 or 8"
           write(6,*)"CHANGE Time_intgr_order=8 or =5"
        stop "Integrator RK5/8"
        endif
 
        if(Time_intgr_order.eq.5) then
      call rk5(Psi,OPsi,PsiDim,xIntPeriod,RelTime,xInitStep,&
                    TolError,iErrorCode,AuxPsi(1:PsiDim,1:7),FUNCSTR,&
                   .false.,dummypost,NSteps)
        elseif (Time_intgr_order.eq.8) then
      call rk8(Psi,OPsi,PsiDim,xIntPeriod,RelTime,xInitStep,&
                    TolError,iErrorCode,AuxPsi(1:PsiDim,1:12),FUNCSTR,&
                    .false.,dummypost,NSteps)
        endif 
        Reltime=Reltime+xIntPeriod 
        CALL FUNCSTR(Reltime,Psi,Opsi) 
        IF(iErrorCode.ne.0) write(6,*)"IN RK !!!!!!! I -ERROR",iErrorCode
        Intgr_total_steps=Intgr_total_steps+NSteps
         IF(print_level.ge.3) then
        write(6,'(a30,i3,a22,i6,i16)')" Runge-Kutta of" &
       ,Time_intgr_order," order takes SmallSteps", NSteps,Intgr_total_steps
       xIntPeriod=xIntPeriod_ORG
         ENDIF !PRINT_LEVEL
       endIF
!c================================ END RUNGE KUTTA==============================
!c================================ END RUNGE KUTTA==============================
!c================================ END RUNGE KUTTA==============================

      IF (trim(Time_intgr).eq.'BS') then
        if(Time_intgr_order.gt.16) then 
     write(6,*)"BULIRSCH-STOER IS IMPLEMENTED TILL 16 ORDER"
     write(6,*)"CHANGE Time_intgr_order<=16"
        stop "Integrator BS"
        endif
89      continue 
            SmallSteps=0
          call  BSStep(Psi,OPsi,PsiDim,xIntPeriod,RelTime,IntOrder, &
                  xInitStep,TolError,& 
           ActualLargeStep, NextLargeStep,SmallSteps, &
                  iErrorCode,AuxPsi,FUNCSTR,AbsBSError,PolyExtrapol)
         Intgr_total_steps=Intgr_total_steps+SmallSteps
        IF((xIntPeriod-ActualLargeStep).ge.TolError) then
       RelTime=RelTime+ActualLargeStep
       CALL FUNCSTR(Reltime,Psi,Opsi) 
         xIntPeriod=xIntPeriod-ActualLargeStep
         xInitStep=min(NextLargeStep,xIntPeriod)
      IF(iErrorCode.ne.0) write(6,*)"IN BS !!!!!!! I -ERROR",iErrorCode
         goto 89
        ENDIF
         IF(print_level.ge.3) then
       write(6,'(a30,i3,a22,i6,i16)')" BULIRSCH-STOER of" &
       ,Time_intgr_order," order takes SmallSteps", SmallSteps,Intgr_total_steps
         xIntPeriod=xIntPeriod_ORG
          ENDIF !PRINT_LEVEL
        endIF
!c================================ BS END ==============================
!c================================ BS END ==============================
!c================================ ZVODE ==============================
!c================================ ZVODE ==============================
      IF (trim(Time_intgr).eq.'STIFF') then
891      continue 
      call  ZVODE_wrapper(FUNCSTR,PsiDim,Psi,RelTime,xIntPeriod,&
                TolError,iErrorCode,AuxPsi,Time_intgr_order)
        write(6,*) 'IN INTGR', Reltime, Abstime, xIntperiod
        IF(xIntperiod-(Reltime-Abstime).gt.TolError) then
           xIntperiod=(AbsTime+xIntPeriod_ORG)-Reltime
        
           IF(iErrorCode.ne.2) write(6,*)"IN ZVODE!!!!!!! I -ERROR",iErrorCode
           Intgr_total_steps=Intgr_total_steps+1
           goto 891
         ENDIF
         IF(print_level.ge.3) then
           write(6,'(a20,i3,i16)')"ZVODE took Steps:",1,Intgr_total_steps
         ENDIF !PRINT_LEVEL
           IF(iErrorCode.eq.2)iErrorCode=0
         xIntPeriod=xIntPeriod_ORG
       endIF
!c================================ ZVODE END ==============================
!c================================ ZVODE END ==============================
!


      IF(iErrorCode.ne.0) write(6,*)"I",iErrorCode,NSteps,NRepeatedSteps

         END SUBROUTINE Integrator_ABM
