!** DGAMIT
REAL(8) FUNCTION DGAMIT(A,X)
  !>
  !  Calculate Tricomi's form of the incomplete Gamma function.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C7E
  !***
  ! **Type:**      DOUBLE PRECISION (GAMIT-S, DGAMIT-D)
  !***
  ! **Keywords:**  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB,
  !             SPECIAL FUNCTIONS, TRICOMI
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  !   Evaluate Tricomi's incomplete Gamma function defined by
  !
  !   DGAMIT = X**(-A)/GAMMA(A) * integral from 0 to X of EXP(-T) *
  !              T**(A-1.)
  !
  !   for A .GT. 0.0 and by analytic continuation for A .LE. 0.0.
  !   GAMMA(X) is the complete gamma function of X.
  !
  !   DGAMIT is evaluated for arbitrary real values of A and for non-
  !   negative values of X (even though DGAMIT is defined for X .LT.
  !   0.0), except that for X = 0 and A .LE. 0.0, DGAMIT is infinite,
  !   which is a fatal error.
  !
  !   The function and both arguments are DOUBLE PRECISION.
  !
  !   A slight deterioration of 2 or 3 digits accuracy will occur when
  !   DGAMIT is very large or very small in absolute value, because log-
  !   arithmic variables are used.  Also, if the parameter  A  is very
  !   close to a negative integer (but not a negative integer), there is
  !   a loss of accuracy, which is reported if the result is less than
  !   half machine precision.
  !
  !***
  ! **References:**  W. Gautschi, A computational procedure for incomplete
  !                 gamma functions, ACM Transactions on Mathematical
  !                 Software 5, 4 (December 1979), pp. 466-481.
  !               W. Gautschi, Incomplete gamma functions, Algorithm 542,
  !                 ACM Transactions on Mathematical Software 5, 4
  !                 (December 1979), pp. 482-489.
  !***
  ! **Routines called:**  D1MACH, D9GMIT, D9LGIC, D9LGIT, DGAMR, DLGAMS,
  !                    DLNGAM, XERCLR, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   920528  DESCRIPTION and REFERENCES sections revised.  (WRB)
  USE service, ONLY : XERMSG, XERCLR, D1MACH
  REAL(8) :: A, X, aeps, ainta, algap1, alng, alx, h, sga, sgngam, t
  REAL(8), SAVE :: alneps, sqeps, bot
  LOGICAL :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  DGAMIT
  IF ( first ) THEN
    alneps = -LOG(D1MACH(3))
    sqeps = SQRT(D1MACH(4))
    bot = LOG(D1MACH(1))
    first = .FALSE.
  END IF
  !
  IF ( X<0.D0 ) CALL XERMSG('SLATEC','DGAMIT','X IS NEGATIVE',2,2)
  !
  IF ( X/=0.D0 ) alx = LOG(X)
  sga = 1.0D0
  IF ( A/=0.D0 ) sga = SIGN(1.0D0,A)
  ainta = AINT(A+0.5D0*sga)
  aeps = A - ainta
  !
  IF ( X<=0.D0 ) THEN
    DGAMIT = 0.0D0
    IF ( ainta>0.D0.OR.aeps/=0.D0 ) DGAMIT = DGAMR(A+1.0D0)
    RETURN
    !
  ELSEIF ( X<=1.D0 ) THEN
    IF ( A>=(-0.5D0).OR.aeps/=0.D0 ) CALL DLGAMS(A+1.0D0,algap1,sgngam)
    DGAMIT = D9GMIT(A,X,algap1,sgngam,alx)
    RETURN
    !
  ELSEIF ( A<X ) THEN
    !
    alng = D9LGIC(A,X,alx)
    !
    ! EVALUATE DGAMIT IN TERMS OF LOG (DGAMIC (A, X))
    !
    h = 1.0D0
    IF ( aeps/=0.D0.OR.ainta>0.D0 ) THEN
      !
      CALL DLGAMS(A+1.0D0,algap1,sgngam)
      t = LOG(ABS(A)) + alng - algap1
      IF ( t>alneps ) THEN
        !
        t = t - A*alx
        IF ( t<bot ) CALL XERCLR
        DGAMIT = -sga*sgngam*EXP(t)
        RETURN
      ELSE
        !
        IF ( t>(-alneps) ) h = 1.0D0 - sga*sgngam*EXP(t)
        IF ( ABS(h)<=sqeps ) THEN
          !
          CALL XERCLR
          CALL XERMSG('SLATEC','DGAMIT','RESULT LT HALF PRECISION',1,1)
        END IF
      END IF
    END IF
  ELSE
    t = D9LGIT(A,X,LOG_GAMMA(A+1.0D0))
    IF ( t<bot ) CALL XERCLR
    DGAMIT = EXP(t)
    RETURN
  END IF
  !
  t = -A*alx + LOG(ABS(h))
  IF ( t<bot ) CALL XERCLR
  DGAMIT = SIGN(EXP(t),h)
  RETURN
END FUNCTION DGAMIT
