!DECK GAMIT
REAL FUNCTION GAMIT(A,X)
  IMPLICIT NONE
  REAL A, aeps, ainta, algap1, alneps, alng, ALNGAM, alx, bot, &
    GAMR, h, R1MACH, R9GMIT, R9LGIC, R9LGIT, sga, sgngam, sqeps, &
    t, X
  !***BEGIN PROLOGUE  GAMIT
  !***PURPOSE  Calculate Tricomi's form of the incomplete Gamma function.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C7E
  !***TYPE      SINGLE PRECISION (GAMIT-S, DGAMIT-D)
  !***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB,
  !             SPECIAL FUNCTIONS, TRICOMI
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  !   Evaluate Tricomi's incomplete gamma function defined by
  !
  !   GAMIT = X**(-A)/GAMMA(A) * integral from 0 to X of EXP(-T) *
  !             T**(A-1.)
  !
  !   for A .GT. 0.0 and by analytic continuation for A .LE. 0.0.
  !   GAMMA(X) is the complete gamma function of X.
  !
  !   GAMIT is evaluated for arbitrary real values of A and for non-
  !   negative values of X (even though GAMIT is defined for X .LT.
  !   0.0), except that for X = 0 and A .LE. 0.0, GAMIT is infinite,
  !   which is a fatal error.
  !
  !   The function and both arguments are REAL.
  !
  !   A slight deterioration of 2 or 3 digits accuracy will occur when
  !   GAMIT is very large or very small in absolute value, because log-
  !   arithmic variables are used.  Also, if the parameter  A  is very
  !   close to a negative integer (but not a negative integer), there is
  !   a loss of accuracy, which is reported if the result is less than
  !   half machine precision.
  !
  !***REFERENCES  W. Gautschi, A computational procedure for incomplete
  !                 gamma functions, ACM Transactions on Mathematical
  !                 Software 5, 4 (December 1979), pp. 466-481.
  !               W. Gautschi, Incomplete gamma functions, Algorithm 542,
  !                 ACM Transactions on Mathematical Software 5, 4
  !                 (December 1979), pp. 482-489.
  !***ROUTINES CALLED  ALGAMS, ALNGAM, GAMR, R1MACH, R9GMIT, R9LGIC,
  !                    R9LGIT, XERCLR, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   920528  DESCRIPTION and REFERENCES sections revised.  (WRB)
  !***END PROLOGUE  GAMIT
  LOGICAL first
  SAVE alneps, sqeps, bot, first
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  GAMIT
  IF ( first ) THEN
    alneps = -LOG(R1MACH(3))
    sqeps = SQRT(R1MACH(4))
    bot = LOG(R1MACH(1))
  ENDIF
  first = .FALSE.
  !
  IF ( X<0.0 ) CALL XERMSG('SLATEC','GAMIT','X IS NEGATIVE',2,2)
  !
  IF ( X/=0.0 ) alx = LOG(X)
  sga = 1.0
  IF ( A/=0.0 ) sga = SIGN(1.0,A)
  ainta = AINT(A+0.5*sga)
  aeps = A - ainta
  !
  IF ( X<=0.0 ) THEN
    GAMIT = 0.0
    IF ( ainta>0.0.OR.aeps/=0.0 ) GAMIT = GAMR(A+1.0)
    RETURN
    !
  ELSEIF ( X<=1.0 ) THEN
    IF ( A>=(-0.5).OR.aeps/=0.0 ) CALL ALGAMS(A+1.0,algap1,sgngam)
    GAMIT = R9GMIT(A,X,algap1,sgngam,alx)
    RETURN
    !
  ELSEIF ( A<X ) THEN
    !
    alng = R9LGIC(A,X,alx)
    !
    ! EVALUATE GAMIT IN TERMS OF LOG(GAMIC(A,X))
    !
    h = 1.0
    IF ( aeps/=0.0.OR.ainta>0.0 ) THEN
      CALL ALGAMS(A+1.0,algap1,sgngam)
      t = LOG(ABS(A)) + alng - algap1
      IF ( t>alneps ) THEN
        !
        t = t - A*alx
        IF ( t<bot ) CALL XERCLR
        GAMIT = -sga*sgngam*EXP(t)
        RETURN
      ELSE
        IF ( t>(-alneps) ) h = 1.0 - sga*sgngam*EXP(t)
        IF ( ABS(h)<=sqeps ) THEN
          CALL XERCLR
          CALL XERMSG('SLATEC','GAMIT','RESULT LT HALF PRECISION',1,1)
        ENDIF
      ENDIF
    ENDIF
  ELSE
    t = R9LGIT(A,X,ALNGAM(A+1.0))
    IF ( t<bot ) CALL XERCLR
    GAMIT = EXP(t)
    RETURN
  ENDIF
  !
  t = -A*alx + LOG(ABS(h))
  IF ( t<bot ) CALL XERCLR
  GAMIT = SIGN(EXP(t),h)
  RETURN
END FUNCTION GAMIT
