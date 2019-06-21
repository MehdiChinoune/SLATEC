!** GAMIT
REAL(SP) FUNCTION GAMIT(A,X)
  !> Calculate Tricomi's form of the incomplete Gamma function.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C7E
  !***
  ! **Type:**      SINGLE PRECISION (GAMIT-S, DGAMIT-D)
  !***
  ! **Keywords:**  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB,
  !             SPECIAL FUNCTIONS, TRICOMI
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  !   Evaluate Tricomi's incomplete gamma function defined by
  !
  !   GAMIT = X**(-A)/GAMMA(A) * integral from 0 to X of EXP(-T) *
  !             T**(A-1.)
  !
  !   for A > 0.0 and by analytic continuation for A <= 0.0.
  !   GAMMA(X) is the complete gamma function of X.
  !
  !   GAMIT is evaluated for arbitrary real values of A and for non-
  !   negative values of X (even though GAMIT is defined for X <
  !   0.0), except that for X = 0 and A <= 0.0, GAMIT is infinite,
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
  !***
  ! **References:**  W. Gautschi, A computational procedure for incomplete
  !                 gamma functions, ACM Transactions on Mathematical
  !                 Software 5, 4 (December 1979), pp. 466-481.
  !               W. Gautschi, Incomplete gamma functions, Algorithm 542,
  !                 ACM Transactions on Mathematical Software 5, 4
  !                 (December 1979), pp. 482-489.
  !***
  ! **Routines called:**  ALGAMS, GAMR, R1MACH, R9GMIT, R9LGIC,
  !                    R9LGIT, XERCLR, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   920528  DESCRIPTION and REFERENCES sections revised.  (WRB)
  USE service, ONLY : XERMSG, XERCLR, R1MACH
  REAL(SP) :: A, X
  REAL(SP) :: aeps, ainta, algap1, alng, alx, h, sga, sgngam, t
  REAL(SP), PARAMETER :: alneps = -LOG(R1MACH(3)), sqeps = SQRT(R1MACH(4)), &
    bot = LOG(R1MACH(1))
  !* FIRST EXECUTABLE STATEMENT  GAMIT
  !
  IF( X<0._SP ) CALL XERMSG('GAMIT','X IS NEGATIVE',2,2)
  !
  IF( X/=0._SP ) alx = LOG(X)
  sga = 1._SP
  IF( A/=0._SP ) sga = SIGN(1._SP,A)
  ainta = AINT(A+0.5_SP*sga)
  aeps = A - ainta
  !
  IF( X<=0._SP ) THEN
    GAMIT = 0._SP
    IF( ainta>0._SP .OR. aeps/=0._SP ) GAMIT = GAMR(A+1._SP)
    RETURN
    !
  ELSEIF( X<=1._SP ) THEN
    IF( A>=(-0.5_SP) .OR. aeps/=0._SP ) CALL ALGAMS(A+1._SP,algap1,sgngam)
    GAMIT = R9GMIT(A,X,algap1,sgngam)
    RETURN
    !
  ELSEIF( A<X ) THEN
    !
    alng = R9LGIC(A,X,alx)
    !
    ! EVALUATE GAMIT IN TERMS OF LOG(GAMIC(A,X))
    !
    h = 1._SP
    IF( aeps/=0._SP .OR. ainta>0._SP ) THEN
      CALL ALGAMS(A+1._SP,algap1,sgngam)
      t = LOG(ABS(A)) + alng - algap1
      IF( t>alneps ) THEN
        !
        t = t - A*alx
        IF( t<bot ) CALL XERCLR
        GAMIT = -sga*sgngam*EXP(t)
        RETURN
      ELSE
        IF( t>(-alneps) ) h = 1._SP - sga*sgngam*EXP(t)
        IF( ABS(h)<=sqeps ) THEN
          CALL XERCLR
          CALL XERMSG('GAMIT','RESULT LT HALF PRECISION',1,1)
        END IF
      END IF
    END IF
  ELSE
    t = R9LGIT(A,X,LOG_GAMMA(A+1._SP))
    IF( t<bot ) CALL XERCLR
    GAMIT = EXP(t)
    RETURN
  END IF
  !
  t = -A*alx + LOG(ABS(h))
  IF( t<bot ) CALL XERCLR
  GAMIT = SIGN(EXP(t),h)
  RETURN
END FUNCTION GAMIT
