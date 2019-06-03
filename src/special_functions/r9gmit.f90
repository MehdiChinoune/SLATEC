!** R9GMIT
REAL FUNCTION R9GMIT(A,X,Algap1,Sgngam)
  !>
  !  Compute Tricomi's incomplete Gamma function for small
  !            arguments.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C7E
  !***
  ! **Type:**      SINGLE PRECISION (R9GMIT-S, D9GMIT-D)
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
  ! **Routines called:**  R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900720  Routine changed from user-callable to subsidiary.  (WRB)
  USE service, ONLY : XERMSG, R1MACH
  REAL :: A, Algap1, Sgngam, X
  INTEGER :: k, m, ma
  REAL :: ae, aeps, alg2, algs, fk, s, sgng2, t, te
  REAL, PARAMETER :: eps = 0.5*R1MACH(3), bot = LOG(R1MACH(1))
  !* FIRST EXECUTABLE STATEMENT  R9GMIT
  !
  IF ( X<=0.0 ) CALL XERMSG('R9GMIT','X SHOULD BE GT 0',1,2)
  !
  ma = INT( A + 0.5 )
  IF ( A<0.0 ) ma = INT( A - 0.5 )
  aeps = A - ma
  !
  ae = A
  IF ( A<(-0.5) ) ae = aeps
  !
  t = 1.0
  te = ae
  s = t
  DO k = 1, 200
    fk = k
    te = -X*te/fk
    t = te/(ae+fk)
    s = s + t
    IF ( ABS(t)<eps*ABS(s) ) GOTO 100
  END DO
  CALL XERMSG('R9GMIT',&
    'NO CONVERGENCE IN 200 TERMS OF TAYLOR-S SERIES',2,2)
  !
  100 CONTINUE
  IF ( A>=(-0.5) ) algs = -Algap1 + LOG(s)
  IF ( A>=(-0.5) ) THEN
    !
    R9GMIT = EXP(algs)
  ELSE
    !
    algs = -LOG_GAMMA(1.0+aeps) + LOG(s)
    s = 1.0
    m = -ma - 1
    IF ( m/=0 ) THEN
      t = 1.0
      DO k = 1, m
        t = X*t/(aeps-m-1+k)
        s = s + t
        IF ( ABS(t)<eps*ABS(s) ) EXIT
      END DO
    END IF
    !
    R9GMIT = 0.0
    algs = -ma*LOG(X) + algs
    IF ( s==0.0.OR.aeps==0.0 ) THEN
      R9GMIT = EXP(algs)
    ELSE
      !
      sgng2 = Sgngam*SIGN(1.0,s)
      alg2 = -X - Algap1 + LOG(ABS(s))
      !
      IF ( alg2>bot ) R9GMIT = sgng2*EXP(alg2)
      IF ( algs>bot ) R9GMIT = R9GMIT + EXP(algs)
      RETURN
    END IF
  END IF
  !
END FUNCTION R9GMIT
