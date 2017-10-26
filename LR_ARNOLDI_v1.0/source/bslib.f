
C **********************************************************************
C *                                                                    *
C *                        BULIRSCH-STOER (bslib.f)                    *
C *                                                                    *
C * Library module containing a Bulirsch-Stoer integrator.             *
C *                                                                    *
C * Contains:                                                          *
C *   BSStep:       The BS integration routine.                        *
C *   ModMidPoint:  Integrates a system of ODEs with modified midpoint *
C *                 method.                                            *
C *   AbsBSError:   Computes the absolute error of the current BS      *
C *                 integration step.                                  *
C *   RelBSError:   Computes the relative error of the current BS      *
C *                 integration step.                                  *
C *   PolyExtrapol: Extrapolates the estimated values for the solution *
C *                 to stepsize zero using a polynomial extrapolation. *
C *   WriteBSStep:  In the current form WriteBSStep is just a dummy    *
C *                 routine doing absolutely nothing; it is included   *
C *                 here for formal reasons, but can (if desired) be   *
C *                 extended easily such that it writes the size and   *
C *                 error of the current integration step to a file.   *
C *   BSErrorMsg:   Returns for a given error number a corresponding   *
C *                 error message.                                     *
C *                                                                    *
C **********************************************************************


C **********************************************************************
C *                                                                    *
C *                            SUBROUTINE BSSTEP                       *
C *                                                                    *
C * Integrates a system of complex first order differential equations  *
C * employing the Bulirsch-Stoer extrapolation method. The routine     *
C * runs with both variable step size and order. BSStep makes only one *
C * single integration step, so it has to be imbedded into a loop that *
C * calls BSStep until the desired time interval is integrated. The    *
C * ODE is of the form dPsi/dt = Func(AbsTime,Psi) =: DtPsi. All       *
C * computations are performed with double precision.                  *
C *                                                                    *
C * Input parameters:                                                  *
C *   Psi:             The (complex) initial-value vector.             *
C *   DtPsi:           Time derivative of the initial-value vector.    *
C *   PsiDim           Length of Psi and DtPsi vectors.                *
C *   IntPeriod:       Lenght of time interval to be integrated.       *
C *   AbsTime:         Absolute time, i. e. Psi_initial=Psi(AbsTime).  *
C *   IntOrder:        Maximum integration order.                      *
C *   LargeStep:       Suggestion for the size of the large step (can  *
C *                    be equal to IntPeriod if no better value is     *
C *                    known).                                         *
C *   TolError:        Maximum error that is tolerated.                *
C *                                                                    *
C * Output parameters:                                                 *
C *   Psi:             Solution of the ODE at time                     *
C *                    AbsTime+ActualLargeStep.                        *
C *   DtPsi:           DtPsi is undetermined.                          *
C *   IntOrder         Optimal column of the extrapolation tableau.    *
C *   ActualLargeStep: Time interval that actually has been integrated *
C *                    (can be lower than LargeStep).                  *
C *   NextLargeStep:   Suggestion for the next large integration step  *
C *                    (can be used as value for LargeStep when BSStep *
C *                    integrates the next large step).                *
C *   SmallSteps:      Counter for the number of small integration     *
C *                    steps. (Note that SmallSteps isn't set to zero  *
C *                    at the beginning.)                              *
C *   ErrorCode:       Error code having the following meaning:        *
C *                    0: everything was o. k.,                        *
C *                    1: illegal integration order,                   *
C *                    2: stepsize underflow.                          *
C *                                                                    *
C * Other parameters:                                                  *
C *   AuxPsi:          Auxiliary array of (minimum) size               *
C *                    PsiDim*(IntOrder+2).                            *
C *                                                                    *
C * External routines:                                                 *
C *   Func:            Computes the time derivative DtPsi of Psi at    *
C *                    time AbsTime. Called as                         *
C *                    Func(AbsTime,Psi,DtPsi).                        *
C *   CalcError:       Determines the error of the current integration *
C *                    step. Called as                                 *
C *                    CalcError(Psi,Delta_Psi,PsiDim,Error).          *
C *   Extrapol:        Extrapolates the estimated values for Psi to    *
C *                    stepsize zero. Called as                        *
C *                    Extrapol(Column,Psi,Delta_Psi,AuxPsi,Tableau,   *
C *                             PsiDim,StepSize**2,IntOrder).          *
C *   WriteStep:       Writes the stepsize and error to a file.        *
C *                    Called as:                                      *
C *                    WriteStep(Iteration,Steps,Stepsize,Error).      *
C *                                                                    *
C * V6.0 MB                                                            *
C * V7.0 GW addition of CData,RData,IData,LData arrays                 *
C *                                                                    *
C **********************************************************************

      Subroutine BSStep (Psi,DtPsi,PsiDim,IntPeriod,AbsTime,IntOrder,
     +                   LargeStep,TolError,ActualLargeStep,
     +                   NextLargeStep,SmallSteps,ErrorCode,AuxPsi,
     +                   Func,CalcError,Extrapol)

      Implicit None

      Real*8    RelativeMinStep,MaxScale,ErrorSafety,ReductionSafety,
     +          MinReduction,MaxReduction,Tiny
      Integer   MaxOrder
      Parameter (RelativeMinStep = 1.0D-10,MaxScale = 0.1D0,
     +          ErrorSafety = 0.25D0,ReductionSafety = 0.7D0,
     +          MinReduction = 0.5D0,MaxReduction = 1.0D-5,
     +          Tiny = 1.0D-25,MaxOrder = 16)

      Integer    PsiDim,IntOrder,SmallSteps,ErrorCode
      Real*8     IntPeriod,AbsTime,LargeStep,TolError,ActualLargeStep,
     +           NextLargeStep
      Complex*16 Psi(PsiDim),DtPsi(PsiDim),AuxPsi(Psidim,IntOrder+2)
      External   Func,CalcError,Extrapol

      Logical CallIsFirstCall,StepIsReduced
      Integer OptimalColumn,OldIntOrder,NumberOfSteps(MaxOrder+1),D,P,
     +        P1,P2
      Real*8  OldAccuracy,IntAccuracy,NewTime,MaxError,MinLargeStep,
     +        SquaredStepSize,Reduction,Work,WorkFactor,MinWork,
     +        ScaleFactor,Error(MaxOrder-1),Expense(MaxOrder+1),
     +        Correction(MaxOrder,MaxOrder)

      Save    OldAccuracy,IntAccuracy,NewTime,Expense,Correction,
     +        OptimalColumn,OldIntOrder,NumberOfSteps,CallIsFirstCall

      Data    OldAccuracy /-1.0D0/,OldIntOrder /-1/,
     +        CallIsFirstCall /.True./


C --- CHECK INTEGRATION ORDER ---
      If ((IntOrder .LT. 2) .Or. (IntOrder .GT. MaxOrder)) Then
         ErrorCode = 1
         Return
      EndIf


C --- INITIALISE VARIABLES ---
      ErrorCode = 0
      MinLargeStep = RelativeMinStep*IntPeriod
      NextLargeStep = -Tiny
      NewTime = -Tiny


C --- REINITIALIZE IF TOLERANCE OR ORDER IS NEW ---
      If ((TolError .NE. OldAccuracy) .Or.
     +   (IntOrder .NE. OldIntOrder)) Then

C    --- INITIALISE VARIABLES ---
         OldAccuracy = TolError
         IntAccuracy = ErrorSafety*TolError
         OldIntOrder = IntOrder
         Do P = 1,IntOrder+1
            NumberOfSteps(P) = 2*P
         EndDo

C    --- COMPUTE EXPENSE COEFFICIENTS ---
         Expense(1) = NumberOfSteps(1)+1
         Do P = 2,IntOrder+1
            Expense(P) = Expense(P-1)+NumberOfSteps(P)
         EndDo

C    --- COMPUTE CORRECTION FACTOR ---
         Do P = 2,IntOrder
            Do P1 = 1,P-1
               Correction(P1,P) = IntAccuracy**((Expense(P1+1)
     +          -Expense(P+1))/((2*P1+1)*(Expense(P+1)
     +          -Expense(1)+1.0D0)))
            EndDo
         EndDo

C    --- DETERMINE OPTIMAL ROW NUMBER FOR CONVERGENCE ---
         Do OptimalColumn = 2,IntOrder-1
            If (Expense(OptimalColumn+1) .GT. Expense(OptimalColumn)
     +         *Correction(OptimalColumn-1,OptimalColumn)) Then
               Goto 10
            EndIf
         EndDo
 10      IntOrder = OptimalColumn
      EndIf

C --- INITIALIZE VARIABLES ---
      ActualLargeStep = LargeStep
      StepIsReduced = .False.
      MinWork = 1.0D35
      Do D = 1,PsiDim
         AuxPsi(D,1) = Psi(D)
      EndDo

C --- RE-ESTABLISH THE ORDER WINDOW ---
      If ((ActualLargeStep .NE. NextLargeStep) .Or.
     +   (AbsTime .NE. NewTime)) Then
         CallIsFirstCall = .True.
         OptimalColumn = IntOrder
      EndIf


C --- MAKE ONE LARGE INTEGRATION STEP ---
 300  Continue


C --- LOOP OVER EACH COLUMN IN THE EXTRAPOLATION TABLEAU ---
      Do P = 1,IntOrder

C    --- INITIALIZE VARIABLES ---
         NewTime = AbsTime+ActualLargeStep
         SquaredStepSize = (ActualLargeStep/NumberOfSteps(P))**2
         SmallSteps = SmallSteps+NumberOfSteps(P)
         
C    --- CHECK IF STEPSIZE IS TOO SMALL ---
         If ((ActualLargeStep .LT. MinLargeStep) .And.
     +      (NewTime .LT. IntPeriod)) Then
            ErrorCode = 2
            Return
         EndIf

C    --- INTEGRATE WITH MODIFIED MIDPOINT METHOD ---

         Call ModMidpoint(AuxPsi(1,1),Psi,DtPsi,AuxPsi(1,3),AuxPsi(1,2),
     +                    PsiDim,ActualLargeStep,AbsTime,
     +                    NumberOfSteps(P),Func)

C    --- EXTRAPOLATE TO ZERO STEP SIZE ---

         Call Extrapol(P,Psi,AuxPsi(1,2),AuxPsi(1,3),AuxPsi(1,4),PsiDim,
     +                 SquaredStepSize,IntOrder)

C    --- COMPUTE ERROR ---

         MaxError = 0.0D0
         If (P .GE. 2) Then
            Call CalcError(Psi,AuxPsi(1,2),PsiDim,MaxError)
            P1 = P-1
            Error(P1) = (MaxError/IntAccuracy)**(1.0D0/(2*P1+1.0D0))
         EndIf

C    --- WRITE STEPSIZE AND ERROR IF DESIRED ---

         Call WriteBSStep(P,NumberOfSteps(P),ActualLargeStep,
     +                  MaxError,AbsTime)

C    --- CHECK IF COLUMN LIES IN ORDER WINDOW ---

         If ((P .GE. 2) .And. ((P .GE. OptimalColumn-1) .Or.
     +      CallIsFirstCall)) Then
 
C       --- EXIT LOOP WHEN CONVERGED ---

            If (MaxError .LT. TolError) Then
               Goto 500
            EndIf

C       --- CHECK FOR POSSIBLE STEPSIZE REDUCTION ---

            If ((P .Eq. IntOrder) .Or. (P .Eq. OptimalColumn+1)) Then
               Reduction = ReductionSafety/Error(P1)
               Goto 400
            ElseIf (P .Eq. OptimalColumn) Then
               If (Correction(OptimalColumn-1,OptimalColumn)
     +            .LT. Error(P1)) Then
                  Reduction = 1.0D0/Error(P1)
                  Goto 400
               EndIf
            ElseIf (OptimalColumn .Eq. IntOrder) Then
               If (Correction(P1,IntOrder-1) .LT. Error(P1)) Then
                  Reduction = ReductionSafety*Correction(P1,IntOrder-1)
     +                        /Error(P1)
                  Goto 400
               EndIf
            ElseIf (Correction(P1,OptimalColumn) .LT. Error(P1)) Then
               Reduction = Correction(P1,OptimalColumn-1)/Error(P1)
               Goto 400
            EndIf
         EndIf
      EndDo

C --- REDUCE STEPSIZE IF STEP WASN'T SUCCESFUL ---

 400  Continue
      Reduction = Min(Reduction,MinReduction)
      Reduction = Max(Reduction,MaxReduction)
      ActualLargeStep = ActualLargeStep*Reduction
      StepIsReduced = .True.

C --- REPEAT INTEGRATION STEP ---

      Goto 300

C --- FINISH INTEGRATION STEP ---

 500  Continue

C --- UPDATE VARIABLES ---

      CallIsFirstCall = .False.

C --- COMPUTE OPTIMAL ROW AND STEPSIZE FOR CONVERGENCE ---

      Do P2 = 1,P1
         WorkFactor = Max(Error(P2),MaxScale)
         Work = WorkFactor*Expense(P2+1)
         If (Work .LT. MinWork) Then
            ScaleFactor = WorkFactor
            MinWork = Work
            OptimalColumn = P2+1
         EndIf
      EndDo
      NextLargeStep = ActualLargeStep/ScaleFactor

C --- CHECK FOR POSSIBLE ORDER INCREASE ---

      If ((OptimalColumn .GE. P) .And. (OptimalColumn .LT. IntOrder)
     +   .And. (.Not. StepIsReduced)) Then
         WorkFactor = Max(ScaleFactor/Correction(OptimalColumn-1,
     +                OptimalColumn),MaxScale)
         If (Expense(OptimalColumn+1)*WorkFactor .LE. MinWork) Then
            NextLargeStep = ActualLargeStep/WorkFactor
            OptimalColumn = OptimalColumn+1
         EndIf
      EndIf

      Return
      End

C **********************************************************************
C *                                                                    *
C *                         SUBROUTINE MODMIDPOINT                     *
C *                                                                    *
C * Integrates a system of ordinary differential equations employing a *
C * modified midpoint method.                                          *
C *                                                                    *
C * Input parameters:                                                  *
C *   SavedPsi:      Initial value vector.                             *
C *   DtPsi:         Time derivative of SavedPsi.                      *
C *   PsiDim:        Length of the psi vectors.                        *
C *   LargeStepSize: Time interval to be integrated.                   *
C *   InitialTime:   Absolute time, i. e. SavedPsi=Psi(InitialTime).   *
C *   NumberOfSteps: Number of integration steps to be made.           *
C *                                                                    *
C * Output parameters:                                                 *
C *   EstimatedPsi:  Solution of the ODE at time                       *
C *                  InitialTime+LargeStepSize.                        *
C *                                                                    *
C * Other parameters:                                                  *
C *   Psi1,Psi2:     Auxiliary arrays each of size PsiDim.             *
C *                                                                    *
C * External routines:                                                 *
C *   Func:          Computes the time derivative DtPsi of Psi and is  *
C *                  called as                                         *
C *                  Func(Time,Psi,DtPsi).                             *
C * V6.0 MB                                                            *
C * V7.0 GW addition of Cdata,Rdata,Idata,Ldata arrays                 *
C *                                                                    *
C **********************************************************************

      Subroutine ModMidpoint (SavedPsi,EstimatedPsi,DtPsi,Psi1,Psi2,
     +                        PsiDim,LargeStepSize,InitialTime,
     +                        NumberOfSteps,Func)

      Implicit None

      Integer    PsiDim,NumberOfSteps
      Real*8     LargeStepSize,InitialTime
      Complex*16 SavedPsi(PsiDim),EstimatedPsi(PsiDim),DtPsi(PsiDim),
     +           Psi1(PsiDim),Psi2(PsiDim)
      External   Func

      Integer    D,P
      Real*8     StepSize,DoubleStepSize,CurrentTime
      Complex*16 Swap

C --- INITIALIZE VARIABLES ---

      StepSize = LargeStepSize/NumberOfSteps
      DoubleStepSize = 2.0D0*StepSize
      CurrentTime = InitialTime+StepSize

C --- STORE PSI AND FIRST ESTIMATE OF NEXT PSI ---

      Do D = 1,PsiDim
         Psi1(D) = SavedPsi(D)
         Psi2(D) = SavedPsi(D)+StepSize*DtPsi(D)
      EndDo
         
C --- EVALUATE FUNCTION WITH ESTIMATED PSI ---

      Call Func(CurrentTime,Psi2,EstimatedPsi)

C --- LOOP OVER NUMBER OF INTEGRATION STEPS ---

      Do P = 2,NumberOfSteps

C    --- PERFORM NEXT INTEGRATION STEP ---

         CurrentTime = CurrentTime+StepSize
         Do D = 1,PsiDim
            Swap = Psi1(D)+DoubleStepSize*EstimatedPsi(D)
            Psi1(D) = Psi2(D)
            Psi2(D) = Swap
         EndDo

C    --- EVALUATE FUNCTION WITH ESTIMATED PSI ---

         Call Func(CurrentTime,Psi2,EstimatedPsi)
      EndDo

C --- TAKE MEAN VALUE FROM LAST AND NEXT ESTIMATE ---

      Do D = 1,PsiDim
         EstimatedPsi(D) = 0.5D0*(Psi1(D)+Psi2(D)
     +                     +StepSize*EstimatedPsi(D))
      EndDo

      Return
      End

C **********************************************************************
C *                                                                    *
C *                         SUBROUTINE ABSBSERROR                      *
C *                                                                    *
C * Estimates the absolute error of the current Bulirsch-Stoer         *
C * iteration.                                                         *
C *                                                                    *
C * Input parameters:                                                  *
C *   Psi:      Vector containing the current solution. (Here Psi is a *
C *             dummy parameter but may be used e. g. if a relative    *
C *             error is desired.)                                     *
C *   PsiError: Vector containing error estimates for each component   *
C *             of the system of ODEs.                                 *
C *   PsiDim:   Length of PsiError.                                    *
C *                                                                    *
C * Output parameters:                                                 *
C *   Error:    Estimated error.                                       *
C *                                                                    *
C * V6.0 MB                                                            *
C *                                                                    *
C **********************************************************************

      Subroutine AbsBSError (Psi,PsiError,PsiDim,Error)

      Implicit None

      Integer    PsiDim
      Real*8     Error
      Complex*16 Psi(PsiDim),PsiError(PsiDim)

      Integer D

C --- ESTIMATE THE ERROR ---

      Error = 0.0D0
      Do D = 1,PsiDim
         Error = Max(Error,Dble(Abs(PsiError(D))))
      EndDo
         
      Return
      End

C **********************************************************************
C *                                                                    *
C *                         SUBROUTINE RELBSERROR                      *
C *                                                                    *
C * Estimates the relative error of the current Bulirsch-Stoer         *
C * iteration.                                                         *
C *                                                                    *
C * Input parameters:                                                  *
C *   Psi:      Vector containing the current solution.                *
C *   PsiError: Vector containing error estimates for each component   *
C *             of the system of ODEs.                                 *
C *   PsiDim:   Length of PsiError.                                    *
C *                                                                    *
C * Output parameters:                                                 *
C *   Error:    Estimated error.                                       *
C *                                                                    *
C *                                                                    *
C * V6.0 MB                                                            *
C * V7.0 MB addition of Data arrays                                    *
C *                                                                    *
C **********************************************************************

      Subroutine RelBSError (Psi,PsiError,PsiDim,Error)

      Implicit None

      Real*8 Tiny
      Parameter (Tiny = 1.0d-30)

      Integer    PsiDim
      Real*8     Error,Maximum
      Complex*16 Psi(PsiDim),PsiError(PsiDim)

      Integer D

C --- ESTIMATE THE ERROR ---

      Maximum = Tiny
      Error = 0.0D0
      Do D = 1,PsiDim
         Maximum = Max(Maximum,Dble(Abs(Psi(D))))
         Error = Max(Error,Dble(Abs(PsiError(D))))
      EndDo
      Error = Error/Maximum
         
      Return
      End

C **********************************************************************
C *                                                                    *
C *                         SUBROUTINE POLYEXTRAPOL                    *
C *                                                                    *
C * Extrapolates the estimated values for psi to stepsize zero using a *
C * polynomial extrapolation.                                          *
C *                                                                    *
C * Input parameters:                                                  *
C *   Column:          Column in the extrapolation tableau to be       *
C *                    filled.                                         *
C *   Psi:             Vector to be filled in extrapolation tableau.   *
C *   PsiDim:          Length of Psi.                                  *
C *   Tableau:         Extrapolation tableau of size                   *
C *                    PsiDim*(IntOrder-1).                            *
C *   SquaredStepSize: Squared large step size.                        *
C *   IntOrder:        Maximum integration order.                      *
C *                                                                    *
C * Output parameters:                                                 *
C *   Psi:             Extrapolated vector.                            *
C *   DeltaPsi:        Error estimate of extrapolated vector.          *
C *   Tableau:         Updated extrapolation tableau.                  *
C *                                                                    *
C * Other parameters:                                                  *
C *   Difference:      Auxiliary vector of length PsiDim.              *
C *                                                                    *
C * V6.0 MB                                                            *
C *                                                                    *
C **********************************************************************

      Subroutine PolyExtrapol (Column,Psi,DeltaPsi,Difference,Tableau,
     +                         PsiDim,SquaredStepSize,IntOrder)

      Implicit None

      Integer   MaxOrder
      Parameter (MaxOrder = 16)

      Integer    Column,PsiDim,IntOrder
      Real*8     SquaredStepSize
      Complex*16 Psi(PsiDim),DeltaPsi(PsiDim),Difference(PsiDim),
     +           Tableau(PsiDim,IntOrder-1)

      Integer    D,P
      Real*8     Factor,OldFactor,Help,OldSquaredStepSize(MaxOrder)
      Complex*16 Delta,InterimTableau

      Save       OldSquaredStepSize

C --- INITIALIZE VARIABLES ---

      OldSquaredStepSize(Column) = SquaredStepSize
      Do D = 1,PsiDim
         DeltaPsi(D) = Psi(D)
      EndDo
      
C --- FILL IN TABLEAU ---

      If (Column .Eq. 1) Then

C    --- STORE FIRST ESTIMATE IN FIRST COLUMN ---

         Do D = 1,PsiDim
            Tableau(D,1) = Psi(D)
         EndDo
      Else
         
C    --- INITIALIZE COLUMN DIFFERENCE ---

         Do D = 1,PsiDim
            Difference(D) = Psi(D)
         EndDo

C    --- LOOP OVER EACH PREVIOUS COLUMN ---

         Do P = 1,Column-1
            
C       --- INITIALIZE VARIABLES ---

            Help = 1.0D0/(OldSquaredStepSize(Column-P)-SquaredStepSize)
            Factor = SquaredStepSize*Help
            OldFactor = OldSquaredStepSize(Column-P)*Help

C       --- PROPAGATE TABLEAU ONE DIAGONAL FURTHER ---

            Do D = 1,PsiDim
               InterimTableau = Tableau(D,P)
               Tableau(D,P) = DeltaPsi(D)
               Delta = Difference(D)-InterimTableau
               DeltaPsi(D) = Factor*Delta
               Difference(D) = OldFactor*Delta
               Psi(D) = Psi(D)+DeltaPsi(D)
            EndDo
         EndDo

C    --- FILL FINAL COLUMN OF THE TABLEAU ---

         If (Column .LT. IntOrder) Then
            Do D = 1,PsiDim
               Tableau(D,Column) = DeltaPsi(D)
            EndDo
         EndIf
      EndIf

      Return
      End

C **********************************************************************
C *                                                                    *
C *                         SUBROUTINE WRITEBSSTEP                     *
C *                                                                    *
C * Writes for the Bulirsch-Stoer extrapolation integrator the         *
C * stepsize and error of the current integration step to a file.      *
C * (In the current form WriteBSStep is a just dummy routine doing     *
C * absolutely nothing; it is included here for formal reasons, but    *
C * can (if desired) be extended easily such that it writes the size   *
C * and error of the current integration step to a file.)              *
C *                                                                    *
C * Input parameters:                                                  *
C *   Iteration:  Number of the current iteration.                     *
C *   SmallSteps: Number of small integration steps used in the        *
C *               current iteration.                                   *
C *   Stepsize:   Size of current large integration step.              *
C *   Error:      Error of current large integration step.             *
C *                                                                    *
C * Output parameters:                                                 *
C *   none                                                             *
C *                                                                    *
C * V6.0 MB                                                            *
C *                                                                    *
C **********************************************************************

      Subroutine WriteBSStep(Iteration,SmallSteps,Stepsize,Error, 
     . AbsTime)

      Implicit None

      Integer Iteration,SmallSteps
      Real*8  Stepsize,Error,AbsTime

C --- WRITE STEPSIZE AND ERROR ---

      Return
      End

C **********************************************************************
C *                                                                    *
C *                        SUBROUTINE BSERRORMSG                       *
C *                                                                    *
C * Generates for a given error number returned by "BSStep" a          *
C * corresponding error message.                                       *
C *                                                                    *
C * Input parameters:                                                  *
C *   Error: Error code returned by BSStep.                            *
C *                                                                    *
C * Output parameters:                                                 *
C *   Msg:   Error message.                                            *
C *                                                                    *
C * V7.0 MB                                                            *
C *                                                                    *
C **********************************************************************

      Subroutine BSErrorMsg (Error,Msg)

      Implicit None

      Integer      Error
      Character*(*) Msg

C --- GENERATE ERROR MESSAGE ---

      If (Error .Eq. 1) Then
         Msg = 'Illegal integration order'
      ElseIf (Error .Eq. 2) Then
         Msg = 'Stepsize underflow'
      Else
         Msg = 'Unknown error occurred'
      EndIf

      Return
      End
