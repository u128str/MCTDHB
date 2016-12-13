C **********************************************************************
C *                                                                    *
C *                   ADAMS-BASHFORTH-MOULTON (abmlib.f)               *
C *                                                                    *
C * Library module containing an Adams-Bashforth-Moulton predictor-    *
C * corrector integrator.                                              *
C *                                                                    *
C * Contains:                                                          *
C *   ABM:          The ABM integration routine.                       *
C *   AbsABMError:  Computes the absolute error of the current ABM     *
C *                 integration step.                                  *
C *   RelABMError:  Computes the relative error of the current ABM     *
C *                 integration step.                                  *
C *   WriteABMStep: In the current form WriteABMStep is just a dummy   *
C *                 routine doing absolutely nothing; it is included   *
C *                 here for formal reasons, but can (if desired) be   *
C *                 extended easily such that it writes the size and   *
C *                 error of the current integration step to a file.   *
C *   ABMErrorMsg:  Returns for a given error number a corresponding   *
C *                 error message.                                     *
C *                                                                    *
C **********************************************************************


C **********************************************************************
C *                                                                    *
C *                           SUBROUTINE ABM                           *
C *                                                                    *
C * Integrates a system of complex first order differential equations  *
C * employing the Adams-Bashforth-Moulton predictor-corrector method.  *
C * The ABM routine runs with variable step sizes and (except for the  *
C * beginning) a fixed integration order. The ODE is of the form       *
C * dPsi/dt = Func(AbsTime,Psi) =: DtPsi. All computations are         *
C * performed with double precision. The routine allows a previous     *
C * integration to be resumed (e. g. after the integration had been    *
C * stopped to write Psi to a file), as long as the AuxPsi array       *
C * hasn't been overwritten and the integration order and error        *
C * tolerance hasn't been changed. To do so simply set the restart     *
C * flag "RestartABM" to ".false.".                                    *
C *                                                                    *
C * Input parameters:                                                  *
C *   Psi:           The (complex) initial-value vector.               *
C *   DtPsi:         Time derivative of the initial-value vector (must *
C *                  be passed to the integrator regardless of the     *
C *                  flag "RestartABM").                               *
C *   PsiDim         Length of Psi and DtPsi vectors.                  *
C *   IntPeriod:     Lenght of time interval to be integrated.         *
C *   AbsTime:       Absolute time, i. e. Psi_initial=Psi(AbsTime).    *
C *   IntOrder:      Desired integration order.                        *
C *   InitStep:      Size of first integration step.                   *
C *   TolError:      Maximum error that is tolerated.                  *
C *   RestartABM:    Restart flag; if false, a previous integration is *
C *                  continued, otherwise a new integration is started *
C *                  (see comment above for details).                  *
C *                                                                    *
C * Output parameters:                                                 *
C *   Psi:           Solution of the ODE at time AbsTime+IntPeriod.    *
C *   DtPsi:         DtPsi is undetermined.                            *
C *   Steps:         Counter for the overall number of integration     *
C *                  steps. (Note that Steps isn't set to zero at the  *
C *                  beginning.)                                       *
C *   RepeatedSteps: Counter for the number of failed and thus         *
C *                  repeated integration steps. (Note that            *
C *                  RepeatedSteps isn't set to zero at the            *
C *                  beginning.)                                       *
C *   ErrorCode:     Error code having the following meaning:          *
C *                  0: everything was o. k.,                          *
C *                  1: illegal integration order,                     *
C *                  2: stepsize underflow.                            *
C *                                                                    *
C * Other parameters:                                                  *
C *   AuxPsi:        Auxiliary array of (minimum) size                 *
C *                  PsiDim*(IntOrder+1).                              *
C *                                                                    *
C * External routines:                                                 *
C *   Func:          Computes the time derivative DtPsi of Psi at time *
C *                  AbsTime. Called as Func(AbsTime,Psi,DtPsi).       *
C *   CalcError:     Determines the error of the current integration   *
C *                  step. Called as CalcError(Predicted_Psi,          *
C *                  Corrected_Psi,PsiDim,Error,).                     *
C *   WriteStep:     Writes the stepsize and error to a file.          *
C *                  Called as WriteStep(Step,Order,Stepsize,Error).   *
C *                                                                    *
C * V6.0 MB                                                            *
C *                                                                    *
C **********************************************************************
!      module serabm
!      contains
      Subroutine OMPABM (Psi,DtPsi,PsiDim,IntPeriod,AbsTime,IntOrder,
     +                InitStep,TolError,RestartABM,Steps,RepeatedSteps,
     +                ErrorCode,AuxPsi,Func,CalcError)
        use SHARED_DIMS

      Implicit None

      Real*8    One6th,One30th,One210th,RelativeMinStep
      Integer   MaxOrder
      Parameter (One6th = 1.0D0/6.0D0,One30th = 1.0D0/30.0D0,
     +           One210th = 1.0D0/210.0D0,RelativeMinStep = 1.0D-12,
     +           MaxOrder = 8)

      Logical    RestartABM
      Integer    PsiDim,IntOrder,Steps,RepeatedSteps,ErrorCode
      Real*8     IntPeriod,AbsTime,InitStep,TolError
      Complex*16 Psi(PsiDim),DtPsi(PsiDim),AuxPsi(PsiDim,IntOrder+1)
      External   Func,CalcError
      
      REAL*8     realerror,singleerror,boundary

      Logical    StepIsRepeated,CurOrdEqIntOrd,UseOldH1,conv
      Integer    CurrentOrder,D,P,I,K
      Real*8     Time,NextTime,IntError,MinError,Error,D2,D3,D4,D5,D6,
     +           D7,H1,H1Sqr,MinusH1Cube,MinStepSize,Distance(MaxOrder),
     +           H(MaxOrder-1),PredCoef(MaxOrder),CorrCoef(MaxOrder),
     +           SumOfD(2:MaxOrder-2,2:MaxOrder-1),
     +           ProdOfD(2:MaxOrder-2,2:MaxOrder-1)
      Complex*16 InterimAuxPsi(MaxOrder)
      
      complex*16 temp

      Save       IntError,Error,MinError,H,Distance,CurrentOrder,
     +           CurOrdEqIntOrd,StepIsRepeated,UseOldH1

CSTR 
c         write(6,*) "Integrator time ",AbsTime
CSTR 

         conv = .False.
C --- CHECK INTEGRATION ORDER ---

      If (IntOrder .GT. MaxOrder) Then
         ErrorCode = 1
         Return
      EndIf

C --- INITIALISE VARIABLES ---

      ErrorCode = 0
      Time = AbsTime
      MinStepSize = RelativeMinStep*IntPeriod
      
C --- CONTINUE INTEGRATION ---

      If (.Not. RestartABM) Then
         Goto 200
      EndIf

C --- INITIALISE VARIABLES ---

C "UseOldH1" is a flag set to true if the step size is artificially
C shortened in order to fit it to the remaining integration period
C before leaving the routine. This allows the next stepsize to be
C corrected back to the previous value.

      IntError = 0.40D0*TolError
      MinError = 0.01D0*TolError
      Do P = 1,IntOrder-1
         H(P) = InitStep
      EndDo
      CurrentOrder = 2
      CurOrdEqIntOrd = .False.
      StepIsRepeated = .False.
      UseOldH1 = .False.

C --- INITIALISE AUXPSI ---

!$OMP PARALLEL DO 
!$OMP& FIRSTPRIVATE(Psidim)
      Do D = 1,PsiDim
         AuxPsi(D,1) = Psi(D)
         AuxPsi(D,2) = DtPsi(D)
      EndDo
!$OMP END PARALLEL DO

C --- INTEGRATION LOOP ---

 100  Continue

C --- CHECK WETHER STEPSIZE IS TOO SMALL ---

      If ((H(1) .LT. MinStepSize) .And.
     +   (H(1) .LT. AbsTime+IntPeriod-Time)) Then
         ErrorCode = 2
c STR
        write(6,*)InitStep," ABM STEp",H(1),MinStepSize
        write(6,*)InitStep," ABM STEp",AbsTime+IntPeriod-Time
c STR
         Return
      EndIf

C --- CALCULATE PREDICTOR AND CORRECTOR COEFFICIENTS ---

      Do P = 1,CurrentOrder
         If (P .Eq. 1) Then
            H1          = H(1)
            Distance(1) = H1
            PredCoef(1) = H1
            CorrCoef(1) = H1
         ElseIf (P .Eq. 2) Then
            H1Sqr       = H1*H1
            Distance(2) = Distance(1)+H(2)
            PredCoef(2) = 0.5*H1Sqr
            CorrCoef(2) = -PredCoef(2)
         ElseIf (P .Eq. 3) Then
            MinusH1Cube = -H1Sqr*H1
            D2          = H(2)
            Distance(3) = Distance(2)+H(3)
            PredCoef(3) = One6th*H1Sqr*(3.0*D2+2.0*H1)
            CorrCoef(3) = One6th*MinusH1Cube
         ElseIf (P .Eq. 4) Then
            D3          = D2+H(3)
            Distance(4) = Distance(3)+H(4)
            SumOfD(2,3) = D2+D3
            ProdOfD(2,3) = D2*D3
            PredCoef(4) = One6th*H1Sqr*(3.0*ProdOfD(2,3)
     +                    +(2.0*SumOfD(2,3)+1.5*H1)*H1)
            CorrCoef(4) = One6th*MinusH1Cube*(D2+0.5*H1)
         ElseIf (P .Eq. 5) Then
            D4          = D3+H(4)
            Distance(5) = Distance(4)+H(5)
            SumOfD(3,3) = D3
            ProdOfD(3,3) = D3
            Do I = 2,3
               SumOfD(I,4)  = SumOfD(I,3)+D4
               ProdOfD(I,4) = ProdOfD(I,3)*D4
            EndDo
            PredCoef(5) = One30th*H1Sqr*(15.0*ProdOfD(2,4)
     +                    +(10.0*(D2*SumOfD(3,4)+ProdOfD(3,4))
     +                    +(7.5*SumOfD(2,4)+6.0*H1)*H1)*H1)
            CorrCoef(5) = One30th*MinusH1Cube*(5.0*ProdOfD(2,3)
     +                    +(2.5*SumOfD(2,3)+1.5*H1)*H1)
         ElseIf (P .Eq. 6) Then
            D5          = D4+H(5)
            Distance(6) = Distance(5)+H(6)
            SumOfD(4,4) = D4
            ProdOfD(4,4) = D4
            Do I = 2,4
               SumOfD(I,5)  = SumOfD(I,4)+D5
               ProdOfD(I,5) = ProdOfD(I,4)*D5
            EndDo
            PredCoef(6) = One30th*H1Sqr*(15.0*ProdOfD(2,5)
     +                    +(10.0*(D2*(D3*SumOfD(4,5)+ProdOfD(4,5))
     +                    +ProdOfD(3,5))
     +                    +(7.5*(D2*SumOfD(3,5)+D3*SumOfD(4,5)
     +                    +ProdOfD(4,5))
     +                    +(6.0*SumOfD(2,5)+5.0*H1)*H1)*H1)*H1)
            CorrCoef(6) = One30th*MinusH1Cube*(5.0*ProdOfD(2,4)
     +                    +(2.5*(D2*SumOfD(3,4)+ProdOfD(3,4))
     +                    +(1.5*SumOfD(2,4)+H1)*H1)*H1)
         ElseIf (P .Eq. 7) Then
            D6          = D5+H(6)
            Distance(7) = Distance(6)+H(7)
            SumOfD(5,5) = D5
            ProdOfD(5,5) = D5
            Do I = 2,5
               SumOfD(I,6)  = SumOfD(I,5)+D6
               ProdOfD(I,6) = ProdOfD(I,5)*D6
            EndDo
            PredCoef(7) = One210th*H1Sqr*(105.0*ProdOfD(2,6)
     +                    +(70.0*(D2*(D3*(D4*SumOfD(5,6)+ProdOfD(5,6))
     +                    +ProdOfD(4,6))+ProdOfD(3,6))
     +                    +(52.5*(D2*(D3*SumOfD(4,6)+D4*SumOfD(5,6)
     +                    +ProdOfD(5,6))+D3*(D4*SumOfD(5,6)
     +                    +ProdOfD(5,6))+ProdOfD(4,6))
     +                    +(42.0*(D2*SumOfD(3,6)+D3*SumOfD(4,6)
     +                    +D4*SumOfD(5,6)+ProdOfD(5,6))
     +                    +(35.0*SumOfD(2,6)+30.0*H1)*H1)*H1)*H1)*H1)
         CorrCoef(7) = One210th*MinusH1Cube*(35.0*ProdOfD(2,5)
     +                    +(17.5*(D2*(D3*SumOfD(4,5)+ProdOfD(4,5))
     +                    +ProdOfD(3,5))
     +                    +(10.5*(D2*SumOfD(3,5)+D3*SumOfD(4,5)
     +                    +ProdOfD(4,5))+(7.0*SumOfD(2,5)
     +                    +5.0*H1)*H1)*H1)*H1)
         ElseIf (P .Eq. 8) Then
            D7          = D6+H(7)
            SumOfD(6,6) = D6
            ProdOfD(6,6) = D6
            Do I = 2,6
               SumOfD(I,7)  = SumOfD(I,6)+D7
               ProdOfD(I,7) = ProdOfD(I,6)*D7
            EndDo
            PredCoef(8) = One210th*H1Sqr*(105.0*ProdOfD(2,7)
     +                    +(60.0*(D2*(D3*(D4*(D5*SumOfD(6,7)
     +                    +ProdOfD(6,7))+ProdOfD(5,7))+ProdOfD(4,7))
     +                    +ProdOfD(3,7))
     +                    +(52.5*(D2*(D3*(D4*SumOfD(5,7)+D5*SumOfD(6,7))
     +                    +ProdOfD(6,7))+D4*(D5*SumOfD(6,7)
     +                    +ProdOfD(6,7))+ProdOfD(5,7))
     +                    +D3*(D4*(D5*SumOfD(6,7)+ProdOfD(6,7)
     +                    +ProdOfD(5,7))+ProdOfD(4,7))
     +                    +(42.0*(D2*(D3*SumOfD(4,7)+D4*SumOfD(5,7)
     +                    +D5*SumOfD(6,7)+ProdOfD(6,7))
     +                    +D3*(D4*SumOfD(5,7)+D5*SumOfD(6,7)
     +                    +ProdOfD(6,7))+D4*(D5*SumOfD(6,7)
     +                    +ProdOfD(6,7))+ProdOfD(5,7))
     +                    +(35.0*(D2*SumOfD(3,7)+D3*SumOfD(4,7)
     +                    +D4*SumOfD(5,7)+D5*SumOfD(6,7)+ProdOfD(6,7))
     +                    +(30.0*SumOfD(2,7)
     +                    +26.25*H1)*H1)*H1)*H1)*H1)*H1)
            CorrCoef(8) = One210th*MinusH1Cube*(35.0*ProdOfD(2,6)
     +                    +(17.5*(D2*(D3*(D4*SumOfD(5,6)+ProdOfD(5,6))
     +                    +ProdOfD(4,6))+ProdOfD(3,6))
     +                    +(10.5*(D2*(D3*SumOfD(4,6)+D4*SumOfD(5,6)
     +                    +ProdOfD(5,6))+D3*(D4*SumOfD(5,6)
     +                    +ProdOfD(5,6))+ProdOfD(4,6))
     +                    +(7.0*(D2*SumOfD(3,6)+D3*SumOfD(4,6)
     +                    +D4*SumOfD(5,6)+ProdOfD(5,6))+(5.0*SumOfD(2,6)
     +                    +3.75*H1)*H1)*H1)*H1)*H1)
         Else
            ErrorCode = 1
            Return
         EndIf
      EndDo

C --- RESTORE ORIGINAL PSI WHEN STEP IS REPEATED ---

      If (StepIsRepeated) Then
c STR 2010
!$OMP PARALLEL DO 
!$OMP& FIRSTPRIVATE(Psidim)
         Do D = 1,PsiDim
            Psi(D) = AuxPsi(D,1)
         EndDo
!$OMP END PARALLEL DO
         StepIsRepeated = .False.
      EndIf

C --- PREDICT PSI ---	

C As long as the order p is increased, the predictor must be used as a
C p-step method. Later on, the predictor is run as a (p+1)-step method.
      If (CurOrdEqIntOrd) Then
         K=currentorder
      else
         K=currentorder-1
      endif

!$OMP PARALLEL DO 
!$OMP& PRIVATE(D,P)
!$OMP& SHARED(PSI,PSIDIM,PREDCOEF,AUXPSI,K)
      Do D = 1,PsiDim
          Do P = 1,K
            Psi(D) = Psi(D)+PredCoef(P)*AuxPsi(D,P+1)
         EndDo
      EndDo
!$OMP END PARALLEL DO

      NextTime = Time+H(1)
      Call Func(NextTime,Psi,DtPsi)

C --- SAVE PREDICTED PSI IN DTPSI AND CORRECT PSI ---

!$OMP PARALLEL 
      Boundary = 1.0D0
      Error = 0.0D0
!$OMP DO REDUCTION(-:INTERIMAUXPSI) REDUCTION(MAX:ERROR)
!$OMP& FIRSTPRIVATE(PsiDim,CurrentOrder,Distance)
      Do D = 1,PsiDim
         InterimAuxPsi(1) = DtPsi(D)
         Do P = 2,CurrentOrder
            InterimAuxPsi(P) = (InterimAuxPsi(P-1)-AuxPsi(D,P))
     +                         /Distance(P-1)
         enddo
         DtPsi(D) = Psi(D)
         Psi(D) = AuxPsi(D,1)
         Do P = 1,CurrentOrder
            psi(d) = psi(d)+CorrCoef(P)*InterimAuxPsi(P)
         EndDo
C --- COMPUTE ERROR OF BOTH REAL AND IMAGINARY PART ---
         SingleError = DAbs(Dble(DtPsi(D))-Dble(Psi(D)))
         If (DAbs(Dble(Psi(D))) .GT. Boundary) Then
             SingleError = SingleError/DAbs(Dble(Psi(D)))
         EndIf
         Error = Max(Error,SingleError)
         SingleError = DAbs(DImag(DtPsi(D))-DImag(Psi(D)))
         If (DAbs(DImag(Psi(D))) .GT. Boundary) Then
           SingleError = SingleError/DAbs(DImag(Psi(D)))
         EndIf
         Error = Max(Error,SingleError)
      EndDo
!$OMP END DO
!$OMP END PARALLEL
      If (CurrentOrder .LT. IntOrder) Then
         Error = Error*(2*IntOrder+1.0D0)/(2*CurrentOrder+1.0D0)
      EndIf
      Steps = Steps+1
      If (Error .GT. TolError) Then
         H(1) = 0.8D0*H(1)*(IntError/Error)**(1.0D0/(CurrentOrder+1))
         RepeatedSteps = RepeatedSteps+1
         StepIsRepeated = .True.
         Goto 100
      EndIf

C --- WRITE STEPSIZE AND ERROR ---

Cstr      Call WriteStep(Steps,CurrentOrder,H(1),Error,Time)

C --- RETURN WHEN FINISHED ---

      Time = Time+H(1)
      If (Time .GE. AbsTime+IntPeriod-MinStepSize) Then
         Return
      EndIf
cSTR added to fix sunf90 problem 
      If (conv .EQV. .True.) Then
         Return
      EndIf
cSTR added to fix sunf90 problem 

C --- EVALUATE FUNCTION WITH CORRECTED PSI ---

      Call Func(NextTime,Psi,DtPsi)

 200  Continue

C --- INCREASE INTEGRATION ORDER ---

      If (CurrentOrder .LT. IntOrder) Then
         CurrentOrder = CurrentOrder+1
      Else
         CurOrdEqIntOrd = .True.
      EndIf

C --- COMPUTE NEW DEVIDED DIFFERENCES ---

C As long as the order p is increased, the predictor must be used as a
C p-step method. Later on, the predictor is run as a (p+1)-step method.
      If (CurOrdEqIntOrd) Then
         K=Currentorder+1
      Else
         K=Currentorder
      endif
!$OMP PARALLEL DO REDUCTION(-:INTERIMAUXPSI)
!$OMP& FIRSTPRIVATE(PsiDim,K,Distance)
      Do D = 1,PsiDim
            AuxPsi(D,1) = Psi(D)
            InterimAuxPsi(1) = DtPsi(D)
            Do P = 2,K-1
               InterimAuxPsi(P) = (InterimAuxPsi(P-1)-AuxPsi(D,P))
     +                            /Distance(P-1)
            EndDo
            Do P = 2,K
               AuxPsi(D,P) = InterimAuxPsi(P-1)
            EndDo
      EndDo
!$OMP END PARALLEL DO

C --- CALCULATE NEW STEP SIZE ---

C If the previous step had been shortened (see above), the next stepsize
C may be set (close) to the stepsize before that one.

      Error = Max(Error,MinError)
      Do P = CurrentOrder-1,2,-1
         H(P) = H(P-1)
      EndDo
      H(1) = H(1)*(IntError/Error)**(1.0D0/(CurrentOrder+1))
      If (UseOldH1) Then
         If (H(1) .LT. 0.98D0*H(3)) Then
            H(1) = 0.98D0*H(3)
         EndIf
         UseOldH1 = .False.
      EndIf
      If (Time+H(1) .GT. AbsTime+IntPeriod) Then
         H(1) = AbsTime+IntPeriod-Time
         UseOldH1 = .True.
         conv = .True.
      EndIf

C --- CONTINUE INTEGRATION ---

      Goto 100
      
      End subroutine ompabm
!      end module
