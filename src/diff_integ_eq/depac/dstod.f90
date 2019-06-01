!** DSTOD
SUBROUTINE DSTOD(Neq,Y,Yh,Nyh,Yh1,Ewt,Savf,Acor,Wm,Iwm,DF,DJAC)
  !>
  !  Subsidiary to DDEBDF
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (STOD-S, DSTOD-D)
  !***
  ! **Author:**  Watts, H. A., (SNLA)
  !***
  ! **Description:**
  !
  !   DSTOD integrates a system of first order odes over one step in the
  !   integrator package DDEBDF.
  ! ----------------------------------------------------------------------
  ! DSTOD  performs one step of the integration of an initial value
  ! problem for a system of ordinary differential equations.
  ! Note.. DSTOD  is independent of the value of the iteration method
  ! indicator MITER, when this is .NE. 0, and hence is independent
  ! of the type of chord method used, or the Jacobian structure.
  ! Communication with DSTOD  is done with the following variables..
  !
  ! Y      = An array of length .GE. N used as the Y argument in
  !          all calls to DF and DJAC.
  ! NEQ    = Integer array containing problem size in NEQ(1), and
  !          passed as the NEQ argument in all calls to DF and DJAC.
  ! YH     = An NYH by LMAX array containing the dependent variables
  !          and their approximate scaled derivatives, where
  !          LMAX = MAXORD + 1.  YH(I,J+1) contains the approximate
  !          J-th derivative of Y(I), scaled by H**J/FACTORIAL(J)
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
  ! WM,IWM = DOUBLE PRECISION and INTEGER work arrays associated with
  !          matrix operations in chord iteration (MITER .NE. 0).
  ! DPJAC   = Name of routine to evaluate and preprocess Jacobian matrix
  !          if a chord method is being used.
  ! DSLVS   = Name of routine to solve linear system in chord iteration.
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
  !***
  ! **See also:**  DDEBDF
  !***
  ! **Routines called:**  DCFOD, DPJAC, DSLVS, DVNRMS
  !***
  ! COMMON BLOCKS    DDEBD1

  !* REVISION HISTORY  (YYMMDD)
  !   820301  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  !   920422  Changed DIMENSION statement.  (WRB)
  USE DDEBD1, ONLY : CONit, CRAte, EL, ELCo, HOLd, RC, RMAx, TESco, EL0, H, HMIn, &
    HMXi, HU, TN, KSTeps, IALth, IPUp, LMAx, MEO, NQNyh, NSTepj, IER, JSTart, &
    KFLag, L, METh, MITer, MAXord, N, NQ, NST, NFE, NQU
  !
  INTERFACE
    SUBROUTINE DF(X,U,Uprime)
      REAL(8) :: X
      REAL(8) :: U(:), Uprime(:)
    END SUBROUTINE DF
    SUBROUTINE DJAC(X,U,Pd,Nrowpd)
      INTEGER :: Nrowpd
      REAL(8) :: X
      REAL(8) :: U(:), Pd(:,:)
    END SUBROUTINE DJAC
  END INTERFACE
  INTEGER :: Neq, Nyh
  INTEGER :: Iwm(:)
  REAL(8) :: Acor(N), Ewt(N), Savf(N), Wm(:), Y(N), Yh(Nyh,MAXord+1), &
    Yh1(Nyh*MAXord+Nyh)
  INTEGER :: i, i1, iredo, iret, j, jb, m, ncf, newq
  REAL(8) :: dcon, ddn, del, delp, dsm, dup, exdn, exsm, exup, r, rh, rhdn, &
    rhsm, rhup, told
  !
  !     BEGIN BLOCK PERMITTING ...EXITS TO 690
  !        BEGIN BLOCK PERMITTING ...EXITS TO 60
  !* FIRST EXECUTABLE STATEMENT  DSTOD
  KFLag = 0
  told = TN
  ncf = 0
  IF ( JSTart>0 ) GOTO 400
  IF ( JSTart==-1 ) THEN
    !              BEGIN BLOCK PERMITTING ...EXITS TO 30
    !                 ------------------------------------------------------
    !                  THE FOLLOWING BLOCK HANDLES PRELIMINARIES NEEDED WHEN
    !                  JSTART = -1.  IPUP IS SET TO MITER TO FORCE A MATRIX
    !                  UPDATE.  IF AN ORDER INCREASE IS ABOUT TO BE
    !                  CONSIDERED (IALTH = 1), IALTH IS RESET TO 2 TO
    !                  POSTPONE CONSIDERATION ONE MORE STEP.  IF THE CALLER
    !                  HAS CHANGED METH, DCFOD  IS CALLED TO RESET THE
    !                  COEFFICIENTS OF THE METHOD.  IF THE CALLER HAS
    !                  CHANGED MAXORD TO A VALUE LESS THAN THE CURRENT
    !                  ORDER NQ, NQ IS REDUCED TO MAXORD, AND A NEW H CHOSEN
    !                  ACCORDINGLY.  IF H IS TO BE CHANGED, YH MUST BE
    !                  RESCALED.  IF H OR METH IS BEING CHANGED, IALTH IS
    !                  RESET TO L = NQ + 1 TO PREVENT FURTHER CHANGES IN H
    !                  FOR THAT MANY STEPS.
    !                 ------------------------------------------------------
    IPUp = MITer
    LMAx = MAXord + 1
    IF ( IALth==1 ) IALth = 2
    IF ( METh/=MEO ) THEN
      CALL DCFOD(METh,ELCo,TESco)
      MEO = METh
      !              ......EXIT
      IF ( NQ<=MAXord ) THEN
        IALth = L
        iret = 1
        !        ............EXIT
        GOTO 100
      END IF
    ELSEIF ( NQ<=MAXord ) THEN
      GOTO 200
    END IF
    NQ = MAXord
    L = LMAx
    DO i = 1, L
      EL(i) = ELCo(i,NQ)
    END DO
    NQNyh = NQ*Nyh
    RC = RC*EL(1)/EL0
    EL0 = EL(1)
    CONit = 0.5D0/(NQ+2)
    ddn = DVNRMS(N,Savf,Ewt)/TESco(1,L)
    exdn = 1.0D0/L
    rhdn = 1.0D0/(1.3D0*ddn**exdn+0.0000013D0)
    rh = MIN(rhdn,1.0D0)
    iredo = 3
    IF ( H==HOLd ) THEN
      rh = MAX(rh,HMIn/ABS(H))
    ELSE
      rh = MIN(rh,ABS(H/HOLd))
      H = HOLd
    END IF
    GOTO 300
  ELSE
    IF ( JSTart==-2 ) GOTO 200
    !              ---------------------------------------------------------
    !               ON THE FIRST CALL, THE ORDER IS SET TO 1, AND OTHER
    !               VARIABLES ARE INITIALIZED.  RMAX IS THE MAXIMUM RATIO BY
    !               WHICH H CAN BE INCREASED IN A SINGLE STEP.  IT IS
    !               INITIALLY 1.E4 TO COMPENSATE FOR THE SMALL INITIAL H,
    !               BUT THEN IS NORMALLY EQUAL TO 10.  IF A FAILURE OCCURS
    !               (IN CORRECTOR CONVERGENCE OR ERROR TEST), RMAX IS SET AT
    !               2 FOR THE NEXT INCREASE.
    !              ---------------------------------------------------------
    LMAx = MAXord + 1
    NQ = 1
    L = 2
    IALth = 2
    RMAx = 10000.0D0
    RC = 0.0D0
    EL0 = 1.0D0
    CRAte = 0.7D0
    delp = 0.0D0
    HOLd = H
    MEO = METh
    NSTepj = 0
    iret = 3
    !           ------------------------------------------------------------
    !            DCFOD  IS CALLED TO GET ALL THE INTEGRATION COEFFICIENTS
    !            FOR THE CURRENT METH.  THEN THE EL VECTOR AND RELATED
    !            CONSTANTS ARE RESET WHENEVER THE ORDER NQ IS CHANGED, OR AT
    !            THE START OF THE PROBLEM.
    !           ------------------------------------------------------------
    CALL DCFOD(METh,ELCo,TESco)
  END IF
  !           BEGIN BLOCK PERMITTING ...EXITS TO 680
  100 CONTINUE
  DO i = 1, L
    EL(i) = ELCo(i,NQ)
  END DO
  NQNyh = NQ*Nyh
  RC = RC*EL(1)/EL0
  EL0 = EL(1)
  CONit = 0.5D0/(NQ+2)
  SELECT CASE (iret)
    CASE (2)
      rh = MAX(rh,HMIn/ABS(H))
      GOTO 300
    CASE (3)
      GOTO 400
    CASE DEFAULT
  END SELECT
  !              ---------------------------------------------------------
  !               IF H IS BEING CHANGED, THE H RATIO RH IS CHECKED AGAINST
  !               RMAX, HMIN, AND HMXI, AND THE YH ARRAY RESCALED.  IALTH
  !               IS SET TO L = NQ + 1 TO PREVENT A CHANGE OF H FOR THAT
  !               MANY STEPS, UNLESS FORCED BY A CONVERGENCE OR ERROR TEST
  !               FAILURE.
  !              ---------------------------------------------------------
  200 CONTINUE
  IF ( H==HOLd ) GOTO 400
  rh = H/HOLd
  H = HOLd
  iredo = 3
  300  rh = MIN(rh,RMAx)
  rh = rh/MAX(1.0D0,ABS(H)*HMXi*rh)
  r = 1.0D0
  DO j = 2, L
    r = r*rh
    DO i = 1, N
      Yh(i,j) = Yh(i,j)*r
    END DO
  END DO
  H = H*rh
  RC = RC*rh
  IALth = L
  IF ( iredo==0 ) THEN
    RMAx = 10.0D0
    r = 1.0D0/TESco(2,NQU)
    DO i = 1, N
      Acor(i) = Acor(i)*r
    END DO
    !     ...............EXIT
    GOTO 1000
  END IF
  !                 ------------------------------------------------------
  !                  THIS SECTION COMPUTES THE PREDICTED VALUES BY
  !                  EFFECTIVELY MULTIPLYING THE YH ARRAY BY THE PASCAL
  !                  TRIANGLE MATRIX.  RC IS THE RATIO OF NEW TO OLD
  !                  VALUES OF THE COEFFICIENT  H*EL(1).  WHEN RC DIFFERS
  !                  FROM 1 BY MORE THAN 30 PERCENT, IPUP IS SET TO MITER
  !                  TO FORCE DPJAC TO BE CALLED, IF A JACOBIAN IS
  !                  INVOLVED.  IN ANY CASE, DPJAC IS CALLED AT LEAST
  !                  EVERY 20-TH STEP.
  !                 ------------------------------------------------------
  !                    BEGIN BLOCK PERMITTING ...EXITS TO 610
  !                       BEGIN BLOCK PERMITTING ...EXITS TO 490
  400 CONTINUE
  IF ( ABS(RC-1.0D0)>0.3D0 ) IPUp = MITer
  IF ( NST>=NSTepj+20 ) IPUp = MITer
  TN = TN + H
  i1 = NQNyh + 1
  DO jb = 1, NQ
    i1 = i1 - Nyh
    DO i = i1, NQNyh
      Yh1(i) = Yh1(i) + Yh1(i+Nyh)
    END DO
  END DO
  KSTeps = KSTeps + 1
  !                          ---------------------------------------------
  !                           UP TO 3 CORRECTOR ITERATIONS ARE TAKEN.  A
  !                           CONVERGENCE TEST IS MADE ON THE R.M.S. NORM
  !                           OF EACH CORRECTION, WEIGHTED BY THE ERROR
  !                           WEIGHT VECTOR EWT.  THE SUM OF THE
  !                           CORRECTIONS IS ACCUMULATED IN THE VECTOR
  !                           ACOR(I).  THE YH ARRAY IS NOT ALTERED IN THE
  !                           CORRECTOR LOOP.
  !                          ---------------------------------------------
  500  m = 0
  DO i = 1, N
    Y(i) = Yh(i,1)
  END DO
  CALL DF(TN,Y,Savf)
  NFE = NFE + 1
  IF ( IPUp>0 ) THEN
    !                                ---------------------------------------
    !                                 IF INDICATED, THE MATRIX P = I -
    !                                 H*EL(1)*J IS REEVALUATED AND
    !                                 PREPROCESSED BEFORE STARTING THE
    !                                 CORRECTOR ITERATION.  IPUP IS SET TO 0
    !                                 AS AN INDICATOR THAT THIS HAS BEEN
    !                                 DONE.
    !                                ---------------------------------------
    IPUp = 0
    RC = 1.0D0
    NSTepj = NST
    CRAte = 0.7D0
    CALL DPJAC(Neq,Y,Yh,Nyh,Ewt,Acor,Savf,Wm,Iwm,DF,DJAC)
    !                          ......EXIT
    IF ( IER/=0 ) GOTO 800
  END IF
  DO i = 1, N
    Acor(i) = 0.0D0
  END DO
  600 CONTINUE
  IF ( MITer/=0 ) THEN
    !                                   ------------------------------------
    !                                    IN THE CASE OF THE CHORD METHOD,
    !                                    COMPUTE THE CORRECTOR ERROR, AND
    !                                    SOLVE THE LINEAR SYSTEM WITH THAT
    !                                    AS RIGHT-HAND SIDE AND P AS
    !                                    COEFFICIENT MATRIX.
    !                                   ------------------------------------
    DO i = 1, N
      Y(i) = H*Savf(i) - (Yh(i,2)+Acor(i))
    END DO
    CALL DSLVS(Wm,Iwm,Y)
    !                             ......EXIT
    IF ( IER/=0 ) GOTO 700
    del = DVNRMS(N,Y,Ewt)
    DO i = 1, N
      Acor(i) = Acor(i) + Y(i)
      Y(i) = Yh(i,1) + EL(1)*Acor(i)
    END DO
  ELSE
    !                                   ------------------------------------
    !                                    IN THE CASE OF FUNCTIONAL
    !                                    ITERATION, UPDATE Y DIRECTLY FROM
    !                                    THE RESULT OF THE LAST FUNCTION
    !                                    EVALUATION.
    !                                   ------------------------------------
    DO i = 1, N
      Savf(i) = H*Savf(i) - Yh(i,2)
      Y(i) = Savf(i) - Acor(i)
    END DO
    del = DVNRMS(N,Y,Ewt)
    DO i = 1, N
      Y(i) = Yh(i,1) + EL(1)*Savf(i)
      Acor(i) = Savf(i)
    END DO
  END IF
  !                                ---------------------------------------
  !                                 TEST FOR CONVERGENCE.  IF M.GT.0, AN
  !                                 ESTIMATE OF THE CONVERGENCE RATE
  !                                 CONSTANT IS STORED IN CRATE, AND THIS
  !                                 IS USED IN THE TEST.
  !                                ---------------------------------------
  IF ( m/=0 ) CRAte = MAX(0.2D0*CRAte,del/delp)
  dcon = del*MIN(1.0D0,1.5D0*CRAte)/(TESco(2,NQ)*CONit)
  IF ( dcon>1.0D0 ) THEN
    m = m + 1
    !                             ...EXIT
    IF ( m/=3 ) THEN
      !                             ...EXIT
      IF ( m<2.OR.del<=2.0D0*delp ) THEN
        delp = del
        CALL DF(TN,Y,Savf)
        NFE = NFE + 1
        GOTO 600
      END IF
    END IF
  ELSE
    !                                   ------------------------------------
    !                                    THE CORRECTOR HAS CONVERGED.  IPUP
    !                                    IS SET TO -1 IF MITER .NE. 0, TO
    !                                    SIGNAL THAT THE JACOBIAN INVOLVED
    !                                    MAY NEED UPDATING LATER.  THE LOCAL
    !                                    ERROR TEST IS MADE AND CONTROL
    !                                    PASSES TO STATEMENT 500 IF IT
    !                                    FAILS.
    !                                   ------------------------------------
    IF ( MITer/=0 ) IPUp = -1
    IF ( m==0 ) dsm = del/TESco(2,NQ)
    IF ( m>0 ) dsm = DVNRMS(N,Acor,Ewt)/TESco(2,NQ)
    IF ( dsm>1.0D0 ) THEN
      !                                   ------------------------------------
      !                                    THE ERROR TEST FAILED.  KFLAG KEEPS
      !                                    TRACK OF MULTIPLE FAILURES.
      !                                    RESTORE TN AND THE YH ARRAY TO
      !                                    THEIR PREVIOUS VALUES, AND PREPARE
      !                                    TO TRY THE STEP AGAIN.  COMPUTE THE
      !                                    OPTIMUM STEP SIZE FOR THIS OR ONE
      !                                    LOWER ORDER.  AFTER 2 OR MORE
      !                                    FAILURES, H IS FORCED TO DECREASE
      !                                    BY A FACTOR OF 0.2 OR LESS.
      !                                   ------------------------------------
      KFLag = KFLag - 1
      TN = told
      i1 = NQNyh + 1
      DO jb = 1, NQ
        i1 = i1 - Nyh
        DO i = i1, NQNyh
          Yh1(i) = Yh1(i) - Yh1(i+Nyh)
        END DO
      END DO
      RMAx = 2.0D0
      IF ( ABS(H)<=HMIn*1.00001D0 ) THEN
        !                                      ---------------------------------
        !                                       ALL RETURNS ARE MADE THROUGH
        !                                       THIS SECTION.  H IS SAVED IN
        !                                       HOLD TO ALLOW THE CALLER TO
        !                                       CHANGE H ON THE NEXT STEP.
        !                                      ---------------------------------
        KFLag = -1
        !     .................................EXIT
        GOTO 1000
        !                    ...............EXIT
      ELSEIF ( KFLag>-3 ) THEN
        iredo = 2
        rhup = 0.0D0
        !                       ............EXIT
        GOTO 900
        !                    ---------------------------------------------------
        !                     CONTROL REACHES THIS SECTION IF 3 OR MORE FAILURES
        !                     HAVE OCCURRED.  IF 10 FAILURES HAVE OCCURRED, EXIT
        !                     WITH KFLAG = -1.  IT IS ASSUMED THAT THE
        !                     DERIVATIVES THAT HAVE ACCUMULATED IN THE YH ARRAY
        !                     HAVE ERRORS OF THE WRONG ORDER.  HENCE THE FIRST
        !                     DERIVATIVE IS RECOMPUTED, AND THE ORDER IS SET TO
        !                     1.  THEN H IS REDUCED BY A FACTOR OF 10, AND THE
        !                     STEP IS RETRIED, UNTIL IT SUCCEEDS OR H REACHES
        !                     HMIN.
        !                    ---------------------------------------------------
      ELSEIF ( KFLag/=-10 ) THEN
        rh = 0.1D0
        rh = MAX(HMIn/ABS(H),rh)
        H = H*rh
        DO i = 1, N
          Y(i) = Yh(i,1)
        END DO
        CALL DF(TN,Y,Savf)
        NFE = NFE + 1
        DO i = 1, N
          Yh(i,2) = H*Savf(i)
        END DO
        IPUp = MITer
        IALth = 5
        !              ......EXIT
        IF ( NQ==1 ) GOTO 400
        NQ = 1
        L = 2
        iret = 3
        GOTO 100
      ELSE
        !                       ------------------------------------------------
        !                        ALL RETURNS ARE MADE THROUGH THIS SECTION.  H
        !                        IS SAVED IN HOLD TO ALLOW THE CALLER TO CHANGE
        !                        H ON THE NEXT STEP.
        !                       ------------------------------------------------
        KFLag = -1
        !     ..................EXIT
        GOTO 1000
      END IF
    ELSE
      !                                      BEGIN BLOCK
      !                                      PERMITTING ...EXITS TO 360
      !                                         ------------------------------
      !                                          AFTER A SUCCESSFUL STEP,
      !                                          UPDATE THE YH ARRAY.
      !                                          CONSIDER CHANGING H IF IALTH
      !                                          = 1.  OTHERWISE DECREASE
      !                                          IALTH BY 1.  IF IALTH IS THEN
      !                                          1 AND NQ .LT. MAXORD, THEN
      !                                          ACOR IS SAVED FOR USE IN A
      !                                          POSSIBLE ORDER INCREASE ON
      !                                          THE NEXT STEP.  IF A CHANGE
      !                                          IN H IS CONSIDERED, AN
      !                                          INCREASE OR DECREASE IN ORDER
      !                                          BY ONE IS CONSIDERED ALSO.  A
      !                                          CHANGE IN H IS MADE ONLY IF
      !                                          IT IS BY A FACTOR OF AT LEAST
      !                                          1.1.  IF NOT, IALTH IS SET TO
      !                                          3 TO PREVENT TESTING FOR THAT
      !                                          MANY STEPS.
      !                                         ------------------------------
      KFLag = 0
      iredo = 0
      NST = NST + 1
      HU = H
      NQU = NQ
      DO j = 1, L
        DO i = 1, N
          Yh(i,j) = Yh(i,j) + EL(j)*Acor(i)
        END DO
      END DO
      IALth = IALth - 1
      IF ( IALth/=0 ) THEN
        !                                      ...EXIT
        IF ( IALth<=1 ) THEN
          !                                      ...EXIT
          IF ( L/=LMAx ) THEN
            DO i = 1, N
              Yh(i,LMAx) = Acor(i)
            END DO
          END IF
        END IF
        r = 1.0D0/TESco(2,NQU)
        DO i = 1, N
          Acor(i) = Acor(i)*r
        END DO
        !     .................................EXIT
        GOTO 1000
      ELSE
        !                                            ---------------------------
        !                                             REGARDLESS OF THE SUCCESS
        !                                             OR FAILURE OF THE STEP,
        !                                             FACTORS RHDN, RHSM, AND
        !                                             RHUP ARE COMPUTED, BY
        !                                             WHICH H COULD BE
        !                                             MULTIPLIED AT ORDER NQ -
        !                                             1, ORDER NQ, OR ORDER NQ +
        !                                             1, RESPECTIVELY.  IN THE
        !                                             CASE OF FAILURE, RHUP =
        !                                             0.0 TO AVOID AN ORDER
        !                                             INCREASE.  THE LARGEST OF
        !                                             THESE IS DETERMINED AND
        !                                             THE NEW ORDER CHOSEN
        !                                             ACCORDINGLY.  IF THE ORDER
        !                                             IS TO BE INCREASED, WE
        !                                             COMPUTE ONE ADDITIONAL
        !                                             SCALED DERIVATIVE.
        !                                            ---------------------------
        rhup = 0.0D0
        !                       .....................EXIT
        IF ( L/=LMAx ) THEN
          DO i = 1, N
            Savf(i) = Acor(i) - Yh(i,LMAx)
          END DO
          dup = DVNRMS(N,Savf,Ewt)/TESco(3,NQ)
          exup = 1.0D0/(L+1)
          !                       .....................EXIT
          rhup = 1.0D0/(1.4D0*dup**exup+0.0000014D0)
        END IF
        GOTO 900
      END IF
    END IF
  END IF
  !                             ------------------------------------------
  !                              THE CORRECTOR ITERATION FAILED TO
  !                              CONVERGE IN 3 TRIES.  IF MITER .NE. 0 AND
  !                              THE JACOBIAN IS OUT OF DATE, DPJAC IS
  !                              CALLED FOR THE NEXT TRY.  OTHERWISE THE
  !                              YH ARRAY IS RETRACTED TO ITS VALUES
  !                              BEFORE PREDICTION, AND H IS REDUCED, IF
  !                              POSSIBLE.  IF H CANNOT BE REDUCED OR 10
  !                              FAILURES HAVE OCCURRED, EXIT WITH KFLAG =
  !                              -2.
  !                             ------------------------------------------
  !                          ...EXIT
  700 CONTINUE
  IF ( IPUp/=0 ) THEN
    IPUp = MITer
    GOTO 500
  END IF
  800  TN = told
  ncf = ncf + 1
  RMAx = 2.0D0
  i1 = NQNyh + 1
  DO jb = 1, NQ
    i1 = i1 - Nyh
    DO i = i1, NQNyh
      Yh1(i) = Yh1(i) - Yh1(i+Nyh)
    END DO
  END DO
  IF ( ABS(H)<=HMIn*1.00001D0 ) THEN
    KFLag = -2
    !     ........................EXIT
    GOTO 1000
  ELSEIF ( ncf/=10 ) THEN
    rh = 0.25D0
    IPUp = MITer
    iredo = 1
    !                 .........EXIT
    rh = MAX(rh,HMIn/ABS(H))
    GOTO 300
  ELSE
    KFLag = -2
    !     ........................EXIT
    GOTO 1000
  END IF
  900  exsm = 1.0D0/L
  rhsm = 1.0D0/(1.2D0*dsm**exsm+0.0000012D0)
  rhdn = 0.0D0
  IF ( NQ/=1 ) THEN
    ddn = DVNRMS(N,Yh(:,L),Ewt)/TESco(1,NQ)
    exdn = 1.0D0/NQ
    rhdn = 1.0D0/(1.3D0*ddn**exdn+0.0000013D0)
  END IF
  IF ( rhsm>=rhup ) THEN
    IF ( rhsm>=rhdn ) THEN
      newq = NQ
      rh = rhsm
      IF ( KFLag==0.AND.rh<1.1D0 ) THEN
        IALth = 3
        r = 1.0D0/TESco(2,NQU)
        DO i = 1, N
          Acor(i) = Acor(i)*r
        END DO
        !     .....................EXIT
        GOTO 1000
      ELSE
        IF ( KFLag<=-2 ) rh = MIN(rh,0.2D0)
        !                             ------------------------------------------
        !                              IF THERE IS A CHANGE OF ORDER, RESET NQ,
        !                              L, AND THE COEFFICIENTS.  IN ANY CASE H
        !                              IS RESET ACCORDING TO RH AND THE YH ARRAY
        !                              IS RESCALED.  THEN EXIT FROM 680 IF THE
        !                              STEP WAS OK, OR REDO THE STEP OTHERWISE.
        !                             ------------------------------------------
        !                 ............EXIT
        IF ( newq==NQ ) THEN
          rh = MAX(rh,HMIn/ABS(H))
          GOTO 300
        ELSE
          NQ = newq
          L = NQ + 1
          iret = 2
          !           ..................EXIT
          GOTO 100
        END IF
      END IF
    END IF
  ELSEIF ( rhup>rhdn ) THEN
    newq = L
    rh = rhup
    IF ( rh>=1.1D0 ) THEN
      r = EL(L)/L
      DO i = 1, N
        Yh(i,newq+1) = Acor(i)*r
      END DO
      NQ = newq
      L = NQ + 1
      iret = 2
      !           ..................EXIT
      GOTO 100
    ELSE
      IALth = 3
      r = 1.0D0/TESco(2,NQU)
      DO i = 1, N
        Acor(i) = Acor(i)*r
      END DO
      !     ...........................EXIT
      GOTO 1000
    END IF
  END IF
  newq = NQ - 1
  rh = rhdn
  IF ( KFLag<0.AND.rh>1.0D0 ) rh = 1.0D0
  IF ( KFLag==0.AND.rh<1.1D0 ) THEN
    IALth = 3
    r = 1.0D0/TESco(2,NQU)
    DO i = 1, N
      Acor(i) = Acor(i)*r
      !     ..................EXIT
    END DO
  ELSE
    IF ( KFLag<=-2 ) rh = MIN(rh,0.2D0)
    !                          ---------------------------------------------
    !                           IF THERE IS A CHANGE OF ORDER, RESET NQ, L,
    !                           AND THE COEFFICIENTS.  IN ANY CASE H IS
    !                           RESET ACCORDING TO RH AND THE YH ARRAY IS
    !                           RESCALED.  THEN EXIT FROM 680 IF THE STEP
    !                           WAS OK, OR REDO THE STEP OTHERWISE.
    !                          ---------------------------------------------
    !                 .........EXIT
    IF ( newq==NQ ) THEN
      rh = MAX(rh,HMIn/ABS(H))
      GOTO 300
    ELSE
      NQ = newq
      L = NQ + 1
      iret = 2
      !           ...............EXIT
      GOTO 100
    END IF
  END IF
  1000 HOLd = H
  JSTart = 1
  !     ----------------------- END OF SUBROUTINE DSTOD
  !     -----------------------
END SUBROUTINE DSTOD
