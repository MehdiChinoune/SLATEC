!** DGAMIT
REAL(DP) ELEMENTAL FUNCTION DGAMIT(A,X)
  !> Calculate Tricomi's form of the incomplete Gamma function.
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
  !   for A > 0.0 and by analytic continuation for A <= 0.0.
  !   GAMMA(X) is the complete gamma function of X.
  !
  !   DGAMIT is evaluated for arbitrary real values of A and for non-
  !   negative values of X (even though DGAMIT is defined for X <
  !   0.0), except that for X = 0 and A <= 0.0, DGAMIT is infinite,
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
  USE service, ONLY : eps_2_dp
  !
  REAL(DP), INTENT(IN) :: A, X
  !
  REAL(DP) :: aeps, ainta, algap1, alng, alx, h, sga, sgngam, t
  REAL(DP), PARAMETER :: alneps = -LOG(eps_2_dp)
  !* FIRST EXECUTABLE STATEMENT  DGAMIT
  IF( X<0._DP ) ERROR STOP 'DGAMIT : X IS NEGATIVE'
  !
  IF( X/=0._DP ) alx = LOG(X)
  sga = 1._DP
  IF( A/=0._DP ) sga = SIGN(1._DP,A)
  ainta = AINT(A+0.5_DP*sga)
  aeps = A - ainta
  !
  IF( X==0._DP ) THEN
    DGAMIT = 0._DP
    IF( ainta>0._DP .OR. aeps/=0._DP ) DGAMIT = DGAMR(A+1._DP)
  ELSEIF( X<=1._DP ) THEN
    IF( A>=(-0.5_DP) .OR. aeps/=0._DP ) CALL DLGAMS(A+1._DP,algap1,sgngam)
    DGAMIT = D9GMIT(A,X,algap1,sgngam)
  ELSEIF( A<X ) THEN
    alng = D9LGIC(A,X,alx)
    ! EVALUATE DGAMIT IN TERMS OF LOG (DGAMIC (A, X))
    h = 1._DP
    IF( aeps/=0._DP .OR. ainta>0._DP ) THEN
      CALL DLGAMS(A+1._DP,algap1,sgngam)
      t = LOG(ABS(A)) + alng - algap1
      IF( t>alneps ) THEN
        t = t - A*alx
        DGAMIT = -sga*sgngam*EXP(t)
        RETURN
      ELSE
        !
        IF( t>(-alneps) ) h = 1._DP - sga*sgngam*EXP(t)
        ! IF( ABS(h)<=sqeps ) THEN
          ! CALL XERMSG('DGAMIT','RESULT LT HALF PRECISION',1,1)
        ! END IF
      END IF
    END IF
    t = -A*alx + LOG(ABS(h))
    DGAMIT = SIGN(EXP(t),h)
  ELSE
    t = D9LGIT(A,X,LOG_GAMMA(A+1._DP))
    DGAMIT = EXP(t)
  END IF

  RETURN
END FUNCTION DGAMIT