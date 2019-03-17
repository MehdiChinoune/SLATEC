!DECK STOD
SUBROUTINE STOD(Neq,Y,Yh,Nyh,Yh1,Ewt,Savf,Acor,Wm,Iwm,F,JAC,Rpar,Ipar)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  STOD
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to DEBDF
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (STOD-S, DSTOD-D)
  !***AUTHOR  Watts, H. A., (SNLA)
  !***DESCRIPTION
  !
  !   STOD integrates a system of first order odes over one step in the
  !   integrator package DEBDF.
  ! ----------------------------------------------------------------------
  ! STOD  performs one step of the integration of an initial value
  ! problem for a system of ordinary differential equations.
  ! Note.. STOD  is independent of the value of the iteration method
  ! indicator MITER, when this is .NE. 0, and hence is independent
  ! of the type of chord method used, or the Jacobian structure.
  ! Communication with STOD  is done with the following variables..
  !
  ! Y      = An array of length .GE. n used as the Y argument in
  !          all calls to F and JAC.
  ! NEQ    = Integer array containing problem size in NEQ(1), and
  !          passed as the NEQ argument in all calls to F and JAC.
  ! YH     = An NYH by LMAX array containing the dependent variables
  !          and their approximate scaled derivatives, where
  !          LMAX = MAXORD + 1.  YH(I,J+1) contains the approximate
  !          J-th derivative of Y(I), scaled by H**J/Factorial(j)
  !          (J = 0,1,...,NQ).  On entry for the first step, the first
  !          two columns of YH must be set from the initial values.
  ! NYH    = A constant integer .GE. N, the first dimension of YH.
  ! YH1    = A one-dimensional array occupying the same space as YH.
  ! EWT    = An array of N elements with which the estimated local
  !          errors in YH are compared.
  ! SAVF   = An array of working storage, of length N.
  ! ACOR   = A work array of length N, used for the accumulated
  !          corrections.  On a successful return, ACOR(I) contains
  !          the estimated one-step local error in Y(I).
  ! WM,IWM = Real and integer work arrays associated with matrix
  !          operations in chord iteration (MITER .NE. 0).
  ! PJAC   = Name of routine to evaluate and preprocess Jacobian matrix
  !          if a chord method is being used.
  ! SLVS   = Name of routine to solve linear system in chord iteration.
  ! H      = The step size to be attempted on the next step.
  !          H is altered by the error control algorithm during the
  !          problem.  H can be either positive or negative, but its
  !          sign must remain constant throughout the problem.
  ! HMIN   = The minimum absolute value of the step size H to be used.
  ! HMXI   = Inverse of the maximum absolute value of H to be used.
  !          HMXI = 0.0 is allowed and corresponds to an infinite HMAX.
  !          HMIN and HMXI may be changed at any time, but will not
  !          take effect until the next change of H is considered.
  ! TN     = The independent variable. TN is updated on each step taken.
  ! JSTART = An integer used for input only, with the following
  !          values and meanings..
  !               0  Perform the first step.
  !           .GT.0  Take a new step continuing from the last.
  !              -1  Take the next step with a new value of H, MAXORD,
  !                    N, METH, MITER, and/or matrix parameters.
  !              -2  Take the next step with a new value of H,
  !                    but with other inputs unchanged.
  !          On return, JSTART is set to 1 to facilitate continuation.
  ! KFLAG  = a completion code with the following meanings..
  !               0  The step was successful.
  !              -1  The requested error could not be achieved.
  !              -2  Corrector convergence could not be achieved.
  !          A return with KFLAG = -1 or -2 means either
  !          ABS(H) = HMIN or 10 consecutive failures occurred.
  !          On a return with KFLAG negative, the values of TN and
  !          the YH array are as of the beginning of the last
  !          step, and H is the last step size attempted.
  ! MAXORD = The maximum order of integration method to be allowed.
  ! METH/MITER = The method flags.  See description in driver.
  ! N      = The number of first-order differential equations.
  ! ----------------------------------------------------------------------
  !
  !***SEE ALSO  DEBDF
  !***ROUTINES CALLED  CFOD, PJAC, SLVS, VNWRMS
  !***COMMON BLOCKS    DEBDF1
  !***REVISION HISTORY  (YYMMDD)
  !   800901  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  !   920422  Changed DIMENSION statement.  (WRB)
  !***END PROLOGUE  STOD
  INTEGER IOD, Ipar, KSTeps
  REAL Rpar
  EXTERNAL F, JAC
  !
  !LLL. OPTIMIZE
  INTEGER Neq, Nyh, Iwm, i, i1, IALth, IER, IOWnd, iredo, iret, &
    IPUp, j, jb, JSTart, KFLag, L, LMAx, m, MAXord, MEO, &
    METh, MITer, N, ncf, newq, NFE, NJE, NQ, NQNyh, NQU, &
    NST, NSTepj
  REAL Y, Yh, Yh1, Ewt, Savf, Acor, Wm, ROWnd, CONit, CRAte, EL, &
    ELCo, HOLd, RC, RMAx, TESco, EL0, H, HMIn, HMXi, HU, TN, &
    UROund, dcon, ddn, del, delp, dsm, dup, exdn, exsm, exup, &
    r, rh, rhdn, rhsm, rhup, told, VNWRMS
  DIMENSION Y(*), Yh(Nyh,*), Yh1(*), Ewt(*), Savf(*), Acor(*), Wm(*), &
    Iwm(*), Rpar(*), Ipar(*)
  COMMON /DEBDF1/ ROWnd, CONit, CRAte, EL(13), ELCo(13,12), HOLd, RC, &
    RMAx, TESco(3,12), EL0, H, HMIn, HMXi, HU, TN, &
    UROund, IOWnd(7), KSTeps, IOD(6), IALth, IPUp, &
    LMAx, MEO, NQNyh, NSTepj, IER, JSTart, KFLag, L, &
    METh, MITer, MAXord, N, NQ, NST, NFE, NJE, NQU
  !
  !
  !***FIRST EXECUTABLE STATEMENT  STOD
  KFLag = 0
  told = TN
  ncf = 0
  IF ( JSTart>0 ) GOTO 400
  IF ( JSTart==-1 ) THEN
    !-----------------------------------------------------------------------
    ! THE FOLLOWING BLOCK HANDLES PRELIMINARIES NEEDED WHEN JSTART = -1.
    ! IPUP IS SET TO MITER TO FORCE A MATRIX UPDATE.
    ! IF AN ORDER INCREASE IS ABOUT TO BE CONSIDERED (IALTH = 1),
    ! IALTH IS RESET TO 2 TO POSTPONE CONSIDERATION ONE MORE STEP.
    ! IF THE CALLER HAS CHANGED METH, CFOD  IS CALLED TO RESET
    ! THE COEFFICIENTS OF THE METHOD.
    ! IF THE CALLER HAS CHANGED MAXORD TO A VALUE LESS THAN THE CURRENT
    ! ORDER NQ, NQ IS REDUCED TO MAXORD, AND A NEW H CHOSEN ACCORDINGLY.
    ! IF H IS TO BE CHANGED, YH MUST BE RESCALED.
    ! IF H OR METH IS BEING CHANGED, IALTH IS RESET TO L = NQ + 1
    ! TO PREVENT FURTHER CHANGES IN H FOR THAT MANY STEPS.
    !-----------------------------------------------------------------------
    IPUp = MITer
    LMAx = MAXord + 1
    IF ( IALth==1 ) IALth = 2
    IF ( METh/=MEO ) THEN
      CALL CFOD(METh,ELCo,TESco)
      MEO = METh
      IF ( NQ<=MAXord ) THEN
        IALth = L
        iret = 1
        GOTO 100
      ENDIF
    ELSEIF ( NQ<=MAXord ) THEN
      GOTO 200
    ENDIF
    NQ = MAXord
    L = LMAx
    DO i = 1, L
      EL(i) = ELCo(i,NQ)
    ENDDO
    NQNyh = NQ*Nyh
    RC = RC*EL(1)/EL0
    EL0 = EL(1)
    CONit = 0.5E0/(NQ+2)
    ddn = VNWRMS(N,Savf,Ewt)/TESco(1,L)
    exdn = 1.0E0/L
    rhdn = 1.0E0/(1.3E0*ddn**exdn+0.0000013E0)
    rh = MIN(rhdn,1.0E0)
    iredo = 3
    IF ( H==HOLd ) THEN
      rh = MAX(rh,HMIn/ABS(H))
    ELSE
      rh = MIN(rh,ABS(H/HOLd))
      H = HOLd
    ENDIF
    GOTO 300
  ELSE
    IF ( JSTart==-2 ) GOTO 200
    !-----------------------------------------------------------------------
    ! ON THE FIRST CALL, THE ORDER IS SET TO 1, AND OTHER VARIABLES ARE
    ! INITIALIZED.  RMAX IS THE MAXIMUM RATIO BY WHICH H CAN BE INCREASED
    ! IN A SINGLE STEP.  IT IS INITIALLY 1.E4 TO COMPENSATE FOR THE SMALL
    ! INITIAL H, BUT THEN IS NORMALLY EQUAL TO 10.  IF A FAILURE
    ! OCCURS (IN CORRECTOR CONVERGENCE OR ERROR TEST), RMAX IS SET AT 2
    ! FOR THE NEXT INCREASE.
    !-----------------------------------------------------------------------
    LMAx = MAXord + 1
    NQ = 1
    L = 2
    IALth = 2
    RMAx = 10000.0E0
    RC = 0.0E0
    EL0 = 1.0E0
    CRAte = 0.7E0
    delp = 0.0E0
    HOLd = H
    MEO = METh
    NSTepj = 0
    iret = 3
    !-----------------------------------------------------------------------
    ! CFOD  IS CALLED TO GET ALL THE INTEGRATION COEFFICIENTS FOR THE
    ! CURRENT METH.  THEN THE EL VECTOR AND RELATED CONSTANTS ARE RESET
    ! WHENEVER THE ORDER NQ IS CHANGED, OR AT THE START OF THE PROBLEM.
    !-----------------------------------------------------------------------
    CALL CFOD(METh,ELCo,TESco)
  ENDIF
  100 CONTINUE
  DO i = 1, L
    EL(i) = ELCo(i,NQ)
  ENDDO
  NQNyh = NQ*Nyh
  RC = RC*EL(1)/EL0
  EL0 = EL(1)
  CONit = 0.5E0/(NQ+2)
  SELECT CASE (iret)
    CASE (2)
      rh = MAX(rh,HMIn/ABS(H))
      GOTO 300
    CASE (3)
      GOTO 400
    CASE DEFAULT
  END SELECT
  !-----------------------------------------------------------------------
  ! IF H IS BEING CHANGED, THE H RATIO RH IS CHECKED AGAINST
  ! RMAX, HMIN, AND HMXI, AND THE YH ARRAY RESCALED.  IALTH IS SET TO
  ! L = NQ + 1 TO PREVENT A CHANGE OF H FOR THAT MANY STEPS, UNLESS
  ! FORCED BY A CONVERGENCE OR ERROR TEST FAILURE.
  !-----------------------------------------------------------------------
  200 CONTINUE
  IF ( H==HOLd ) GOTO 400
  rh = H/HOLd
  H = HOLd
  iredo = 3
  300  rh = MIN(rh,RMAx)
  rh = rh/MAX(1.0E0,ABS(H)*HMXi*rh)
  r = 1.0E0
  DO j = 2, L
    r = r*rh
    DO i = 1, N
      Yh(i,j) = Yh(i,j)*r
    ENDDO
  ENDDO
  H = H*rh
  RC = RC*rh
  IALth = L
  IF ( iredo==0 ) THEN
    RMAx = 10.0E0
    GOTO 1200
  ENDIF
  !-----------------------------------------------------------------------
  ! THIS SECTION COMPUTES THE PREDICTED VALUES BY EFFECTIVELY
  ! MULTIPLYING THE YH ARRAY BY THE PASCAL TRIANGLE MATRIX.
  ! RC IS THE RATIO OF NEW TO OLD VALUES OF THE COEFFICIENT  H*EL(1).
  ! WHEN RC DIFFERS FROM 1 BY MORE THAN 30 PERCENT, IPUP IS SET TO MITER
  ! TO FORCE PJAC TO BE CALLED, IF A JACOBIAN IS INVOLVED.
  ! IN ANY CASE, PJAC IS CALLED AT LEAST EVERY 20-TH STEP.
  !-----------------------------------------------------------------------
  400 CONTINUE
  IF ( ABS(RC-1.0E0)>0.3E0 ) IPUp = MITer
  IF ( NST>=NSTepj+20 ) IPUp = MITer
  TN = TN + H
  i1 = NQNyh + 1
  DO jb = 1, NQ
    i1 = i1 - Nyh
    DO i = i1, NQNyh
      Yh1(i) = Yh1(i) + Yh1(i+Nyh)
    ENDDO
  ENDDO
  KSTeps = KSTeps + 1
  !-----------------------------------------------------------------------
  ! UP TO 3 CORRECTOR ITERATIONS ARE TAKEN.  A CONVERGENCE TEST IS
  ! MADE ON THE R.M.S. NORM OF EACH CORRECTION, WEIGHTED BY THE ERROR
  ! WEIGHT VECTOR EWT.  THE SUM OF THE CORRECTIONS IS ACCUMULATED IN THE
  ! VECTOR ACOR(I).  THE YH ARRAY IS NOT ALTERED IN THE CORRECTOR LOOP.
  !-----------------------------------------------------------------------
  500  m = 0
  DO i = 1, N
    Y(i) = Yh(i,1)
  ENDDO
  CALL F(TN,Y,Savf,Rpar,Ipar)
  NFE = NFE + 1
  IF ( IPUp>0 ) THEN
    !-----------------------------------------------------------------------
    ! IF INDICATED, THE MATRIX P = I - H*EL(1)*J IS REEVALUATED AND
    ! PREPROCESSED BEFORE STARTING THE CORRECTOR ITERATION.  IPUP IS SET
    ! TO 0 AS AN INDICATOR THAT THIS HAS BEEN DONE.
    !-----------------------------------------------------------------------
    IPUp = 0
    RC = 1.0E0
    NSTepj = NST
    CRAte = 0.7E0
    CALL PJAC(Neq,Y,Yh,Nyh,Ewt,Acor,Savf,Wm,Iwm,F,JAC,Rpar,Ipar)
    IF ( IER/=0 ) GOTO 800
  ENDIF
  DO i = 1, N
    Acor(i) = 0.0E0
  ENDDO
  600 CONTINUE
  IF ( MITer/=0 ) THEN
    !-----------------------------------------------------------------------
    ! IN THE CASE OF THE CHORD METHOD, COMPUTE THE CORRECTOR ERROR,
    ! AND SOLVE THE LINEAR SYSTEM WITH THAT AS RIGHT-HAND SIDE AND
    ! P AS COEFFICIENT MATRIX.
    !-----------------------------------------------------------------------
    DO i = 1, N
      Y(i) = H*Savf(i) - (Yh(i,2)+Acor(i))
    ENDDO
    CALL SLVS(Wm,Iwm,Y,Savf)
    IF ( IER/=0 ) GOTO 700
    del = VNWRMS(N,Y,Ewt)
    DO i = 1, N
      Acor(i) = Acor(i) + Y(i)
      Y(i) = Yh(i,1) + EL(1)*Acor(i)
    ENDDO
  ELSE
    !-----------------------------------------------------------------------
    ! IN THE CASE OF FUNCTIONAL ITERATION, UPDATE Y DIRECTLY FROM
    ! THE RESULT OF THE LAST FUNCTION EVALUATION.
    !-----------------------------------------------------------------------
    DO i = 1, N
      Savf(i) = H*Savf(i) - Yh(i,2)
      Y(i) = Savf(i) - Acor(i)
    ENDDO
    del = VNWRMS(N,Y,Ewt)
    DO i = 1, N
      Y(i) = Yh(i,1) + EL(1)*Savf(i)
      Acor(i) = Savf(i)
    ENDDO
  ENDIF
  !-----------------------------------------------------------------------
  ! TEST FOR CONVERGENCE.  IF M.GT.0, AN ESTIMATE OF THE CONVERGENCE
  ! RATE CONSTANT IS STORED IN CRATE, AND THIS IS USED IN THE TEST.
  !-----------------------------------------------------------------------
  IF ( m/=0 ) CRAte = MAX(0.2E0*CRAte,del/delp)
  dcon = del*MIN(1.0E0,1.5E0*CRAte)/(TESco(2,NQ)*CONit)
  IF ( dcon<=1.0E0 ) THEN
    !-----------------------------------------------------------------------
    ! THE CORRECTOR HAS CONVERGED.  IPUP IS SET TO -1 IF MITER .NE. 0,
    ! TO SIGNAL THAT THE JACOBIAN INVOLVED MAY NEED UPDATING LATER.
    ! THE LOCAL ERROR TEST IS MADE AND CONTROL PASSES TO STATEMENT 500
    ! IF IT FAILS.
    !-----------------------------------------------------------------------
    IF ( MITer/=0 ) IPUp = -1
    IF ( m==0 ) dsm = del/TESco(2,NQ)
    IF ( m>0 ) dsm = VNWRMS(N,Acor,Ewt)/TESco(2,NQ)
    IF ( dsm>1.0E0 ) THEN
      !-----------------------------------------------------------------------
      ! THE ERROR TEST FAILED.  KFLAG KEEPS TRACK OF MULTIPLE FAILURES.
      ! RESTORE TN AND THE YH ARRAY TO THEIR PREVIOUS VALUES, AND PREPARE
      ! TO TRY THE STEP AGAIN.  COMPUTE THE OPTIMUM STEP SIZE FOR THIS OR
      ! ONE LOWER ORDER.  AFTER 2 OR MORE FAILURES, H IS FORCED TO DECREASE
      ! BY A FACTOR OF 0.2 OR LESS.
      !-----------------------------------------------------------------------
      KFLag = KFLag - 1
      TN = told
      i1 = NQNyh + 1
      DO jb = 1, NQ
        i1 = i1 - Nyh
        DO i = i1, NQNyh
          Yh1(i) = Yh1(i) - Yh1(i+Nyh)
        ENDDO
      ENDDO
      RMAx = 2.0E0
      IF ( ABS(H)<=HMIn*1.00001E0 ) THEN
        !-----------------------------------------------------------------------
        ! ALL RETURNS ARE MADE THROUGH THIS SECTION.  H IS SAVED IN HOLD
        ! TO ALLOW THE CALLER TO CHANGE H ON THE NEXT STEP.
        !-----------------------------------------------------------------------
        KFLag = -1
        GOTO 1300
      ELSEIF ( KFLag<=-3 ) THEN
        !-----------------------------------------------------------------------
        ! CONTROL REACHES THIS SECTION IF 3 OR MORE FAILURES HAVE OCCURRED.
        ! IF 10 FAILURES HAVE OCCURRED, EXIT WITH KFLAG = -1.
        ! IT IS ASSUMED THAT THE DERIVATIVES THAT HAVE ACCUMULATED IN THE
        ! YH ARRAY HAVE ERRORS OF THE WRONG ORDER.  HENCE THE FIRST
        ! DERIVATIVE IS RECOMPUTED, AND THE ORDER IS SET TO 1.  THEN
        ! H IS REDUCED BY A FACTOR OF 10, AND THE STEP IS RETRIED,
        ! UNTIL IT SUCCEEDS OR H REACHES HMIN.
        !-----------------------------------------------------------------------
        IF ( KFLag==-10 ) THEN
          KFLag = -1
          GOTO 1300
        ELSE
          rh = 0.1E0
          rh = MAX(HMIn/ABS(H),rh)
          H = H*rh
          DO i = 1, N
            Y(i) = Yh(i,1)
          ENDDO
          CALL F(TN,Y,Savf,Rpar,Ipar)
          NFE = NFE + 1
          DO i = 1, N
            Yh(i,2) = H*Savf(i)
          ENDDO
          IPUp = MITer
          IALth = 5
          IF ( NQ==1 ) GOTO 400
          NQ = 1
          L = 2
          iret = 3
          GOTO 100
        ENDIF
      ELSE
        iredo = 2
        rhup = 0.0E0
        GOTO 900
      ENDIF
    ELSE
      !-----------------------------------------------------------------------
      ! AFTER A SUCCESSFUL STEP, UPDATE THE YH ARRAY.
      ! CONSIDER CHANGING H IF IALTH = 1.  OTHERWISE DECREASE IALTH BY 1.
      ! IF IALTH IS THEN 1 AND NQ .LT. MAXORD, THEN ACOR IS SAVED FOR
      ! USE IN A POSSIBLE ORDER INCREASE ON THE NEXT STEP.
      ! IF A CHANGE IN H IS CONSIDERED, AN INCREASE OR DECREASE IN ORDER
      ! BY ONE IS CONSIDERED ALSO.  A CHANGE IN H IS MADE ONLY IF IT IS BY A
      ! FACTOR OF AT LEAST 1.1.  IF NOT, IALTH IS SET TO 3 TO PREVENT
      ! TESTING FOR THAT MANY STEPS.
      !-----------------------------------------------------------------------
      KFLag = 0
      iredo = 0
      NST = NST + 1
      HU = H
      NQU = NQ
      DO j = 1, L
        DO i = 1, N
          Yh(i,j) = Yh(i,j) + EL(j)*Acor(i)
        ENDDO
      ENDDO
      IALth = IALth - 1
      IF ( IALth==0 ) THEN
        !-----------------------------------------------------------------------
        ! REGARDLESS OF THE SUCCESS OR FAILURE OF THE STEP, FACTORS
        ! RHDN, RHSM, AND RHUP ARE COMPUTED, BY WHICH H COULD BE MULTIPLIED
        ! AT ORDER NQ - 1, ORDER NQ, OR ORDER NQ + 1, RESPECTIVELY.
        ! IN THE CASE OF FAILURE, RHUP = 0.0 TO AVOID AN ORDER INCREASE.
        ! THE LARGEST OF THESE IS DETERMINED AND THE NEW ORDER CHOSEN
        ! ACCORDINGLY.  IF THE ORDER IS TO BE INCREASED, WE COMPUTE ONE
        ! ADDITIONAL SCALED DERIVATIVE.
        !-----------------------------------------------------------------------
        rhup = 0.0E0
        IF ( L/=LMAx ) THEN
          DO i = 1, N
            Savf(i) = Acor(i) - Yh(i,LMAx)
          ENDDO
          dup = VNWRMS(N,Savf,Ewt)/TESco(3,NQ)
          exup = 1.0E0/(L+1)
          rhup = 1.0E0/(1.4E0*dup**exup+0.0000014E0)
        ENDIF
        GOTO 900
      ELSE
        IF ( IALth<=1 ) THEN
          IF ( L/=LMAx ) THEN
            DO i = 1, N
              Yh(i,LMAx) = Acor(i)
            ENDDO
          ENDIF
        ENDIF
        GOTO 1200
      ENDIF
    ENDIF
  ELSE
    m = m + 1
    IF ( m/=3 ) THEN
      IF ( m<2.OR.del<=2.0E0*delp ) THEN
        delp = del
        CALL F(TN,Y,Savf,Rpar,Ipar)
        NFE = NFE + 1
        GOTO 600
      ENDIF
    ENDIF
  ENDIF
  !-----------------------------------------------------------------------
  ! THE CORRECTOR ITERATION FAILED TO CONVERGE IN 3 TRIES.
  ! IF MITER .NE. 0 AND THE JACOBIAN IS OUT OF DATE, PJAC IS CALLED FOR
  ! THE NEXT TRY.  OTHERWISE THE YH ARRAY IS RETRACTED TO ITS VALUES
  ! BEFORE PREDICTION, AND H IS REDUCED, IF POSSIBLE.  IF H CANNOT BE
  ! REDUCED OR 10 FAILURES HAVE OCCURRED, EXIT WITH KFLAG = -2.
  !-----------------------------------------------------------------------
  700 CONTINUE
  IF ( IPUp/=0 ) THEN
    IPUp = MITer
    GOTO 500
  ENDIF
  800  TN = told
  ncf = ncf + 1
  RMAx = 2.0E0
  i1 = NQNyh + 1
  DO jb = 1, NQ
    i1 = i1 - Nyh
    DO i = i1, NQNyh
      Yh1(i) = Yh1(i) - Yh1(i+Nyh)
    ENDDO
  ENDDO
  IF ( ABS(H)<=HMIn*1.00001E0 ) THEN
    KFLag = -2
    GOTO 1300
  ELSEIF ( ncf==10 ) THEN
    KFLag = -2
    GOTO 1300
  ELSE
    rh = 0.25E0
    IPUp = MITer
    iredo = 1
    rh = MAX(rh,HMIn/ABS(H))
    GOTO 300
  ENDIF
  900  exsm = 1.0E0/L
  rhsm = 1.0E0/(1.2E0*dsm**exsm+0.0000012E0)
  rhdn = 0.0E0
  IF ( NQ/=1 ) THEN
    ddn = VNWRMS(N,Yh(1,L),Ewt)/TESco(1,NQ)
    exdn = 1.0E0/NQ
    rhdn = 1.0E0/(1.3E0*ddn**exdn+0.0000013E0)
  ENDIF
  IF ( rhsm>=rhup ) THEN
    IF ( rhsm>=rhdn ) THEN
      newq = NQ
      rh = rhsm
      GOTO 1000
    ENDIF
  ELSEIF ( rhup>rhdn ) THEN
    newq = L
    rh = rhup
    IF ( rh<1.1E0 ) THEN
      IALth = 3
      GOTO 1200
    ELSE
      r = EL(L)/L
      DO i = 1, N
        Yh(i,newq+1) = Acor(i)*r
      ENDDO
      GOTO 1100
    ENDIF
  ENDIF
  newq = NQ - 1
  rh = rhdn
  IF ( KFLag<0.AND.rh>1.0E0 ) rh = 1.0E0
  1000 CONTINUE
  IF ( (KFLag==0).AND.(rh<1.1E0) ) THEN
    IALth = 3
    GOTO 1200
  ELSE
    IF ( KFLag<=-2 ) rh = MIN(rh,0.2E0)
    !-----------------------------------------------------------------------
    ! IF THERE IS A CHANGE OF ORDER, RESET NQ, L, AND THE COEFFICIENTS.
    ! IN ANY CASE H IS RESET ACCORDING TO RH AND THE YH ARRAY IS RESCALED.
    ! THEN EXIT FROM 680 IF THE STEP WAS OK, OR REDO THE STEP OTHERWISE.
    !-----------------------------------------------------------------------
    IF ( newq==NQ ) THEN
      rh = MAX(rh,HMIn/ABS(H))
      GOTO 300
    ENDIF
  ENDIF
  1100 NQ = newq
  L = NQ + 1
  iret = 2
  GOTO 100
  1200 r = 1.0E0/TESco(2,NQU)
  DO i = 1, N
    Acor(i) = Acor(i)*r
  ENDDO
  1300 HOLd = H
  JSTart = 1
  !----------------------- END OF SUBROUTINE STOD  -----------------------
END SUBROUTINE STOD
