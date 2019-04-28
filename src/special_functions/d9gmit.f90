!** D9GMIT
REAL(8) FUNCTION D9GMIT(A,X,Algap1,Sgngam)
  !>
  !  Compute Tricomi's incomplete Gamma function for small
  !            arguments.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C7E
  !***
  ! **Type:**      DOUBLE PRECISION (R9GMIT-S, D9GMIT-D)
  !***
  ! **Keywords:**  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB, SMALL X,
  !             SPECIAL FUNCTIONS, TRICOMI
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! Compute Tricomi's incomplete gamma function for small X.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  D1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890911  Removed unnecessary intrinsics.  (WRB)
  !   890911  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900720  Routine changed from user-callable to subsidiary.  (WRB)
  USE service, ONLY : XERMSG, D1MACH
  INTEGER k, m, ma
  REAL(8) :: A, X, Algap1, Sgngam, ae, aeps, algs, alg2, fk, s, sgng2, t, te
  REAL(8), SAVE :: eps, bot
  LOGICAL :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  D9GMIT
  IF ( first ) THEN
    eps = 0.5D0*D1MACH(3)
    bot = LOG(D1MACH(1))
    first = .FALSE.
  END IF
  !
  IF ( X<=0.D0 ) CALL XERMSG('SLATEC','D9GMIT','X SHOULD BE GT 0',1,2)
  !
  ma = INT( A + 0.5D0 )
  IF ( A<0.D0 ) ma = INT( A - 0.5D0 )
  aeps = A - ma
  !
  ae = A
  IF ( A<(-0.5D0) ) ae = aeps
  !
  t = 1.D0
  te = ae
  s = t
  DO k = 1, 200
    fk = k
    te = -X*te/fk
    t = te/(ae+fk)
    s = s + t
    IF ( ABS(t)<eps*ABS(s) ) GOTO 100
  END DO
  CALL XERMSG('SLATEC','D9GMIT',&
    'NO CONVERGENCE IN 200 TERMS OF TAYLOR-S SERIES',2,2)
  !
  100 CONTINUE
  IF ( A>=(-0.5D0) ) algs = -Algap1 + LOG(s)
  IF ( A>=(-0.5D0) ) THEN
    !
    D9GMIT = EXP(algs)
  ELSE
    !
    algs = -LOG_GAMMA(1.D0+aeps) + LOG(s)
    s = 1.0D0
    m = -ma - 1
    IF ( m/=0 ) THEN
      t = 1.0D0
      DO k = 1, m
        t = X*t/(aeps-(m+1-k))
        s = s + t
        IF ( ABS(t)<eps*ABS(s) ) EXIT
      END DO
    END IF
    !
    D9GMIT = 0.0D0
    algs = -ma*LOG(X) + algs
    IF ( s==0.D0.OR.aeps==0.D0 ) THEN
      D9GMIT = EXP(algs)
    ELSE
      !
      sgng2 = Sgngam*SIGN(1.0D0,s)
      alg2 = -X - Algap1 + LOG(ABS(s))
      !
      IF ( alg2>bot ) D9GMIT = sgng2*EXP(alg2)
      IF ( algs>bot ) D9GMIT = D9GMIT + EXP(algs)
      RETURN
    END IF
  END IF
  !
END FUNCTION D9GMIT
