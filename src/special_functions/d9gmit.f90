!** D9GMIT
REAL(DP) ELEMENTAL FUNCTION D9GMIT(A,X,Algap1,Sgngam)
  !> Compute Tricomi's incomplete Gamma function for small
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
  USE service, ONLY : D1MACH
  REAL(DP), INTENT(IN) :: A, X, Algap1, Sgngam
  INTEGER :: k, m, ma
  REAL(DP) :: ae, aeps, algs, alg2, fk, s, sgng2, t, te
  REAL(DP), PARAMETER :: eps = 0.5_DP*D1MACH(3), bot = LOG(D1MACH(1))
  !* FIRST EXECUTABLE STATEMENT  D9GMIT
  !
  IF( X<=0._DP ) ERROR STOP 'D9GMIT : X SHOULD BE > 0'
  !
  ma = INT( A + 0.5_DP )
  IF( A<0._DP ) ma = INT( A - 0.5_DP )
  aeps = A - ma
  !
  ae = A
  IF( A<(-0.5_DP) ) ae = aeps
  !
  t = 1._DP
  te = ae
  s = t
  DO k = 1, 200
    fk = k
    te = -X*te/fk
    t = te/(ae+fk)
    s = s + t
    IF( ABS(t)<eps*ABS(s) ) GOTO 100
  END DO
  ERROR STOP 'D9GMIT : NO CONVERGENCE IN 200 TERMS OF TAYLOR-S SERIES'
  !
  100 CONTINUE
  IF( A>=(-0.5_DP) ) THEN
    algs = -Algap1 + LOG(s)
    D9GMIT = EXP(algs)
  ELSE
    !
    algs = -LOG_GAMMA(1._DP+aeps) + LOG(s)
    s = 1._DP
    m = -ma - 1
    IF( m/=0 ) THEN
      t = 1._DP
      DO k = 1, m
        t = X*t/(aeps-(m+1-k))
        s = s + t
        IF( ABS(t)<eps*ABS(s) ) EXIT
      END DO
    END IF
    !
    D9GMIT = 0._DP
    algs = -ma*LOG(X) + algs
    IF( s==0._DP .OR. aeps==0._DP ) THEN
      D9GMIT = EXP(algs)
    ELSE
      !
      sgng2 = Sgngam*SIGN(1._DP,s)
      alg2 = -X - Algap1 + LOG(ABS(s))
      !
      IF( alg2>bot ) THEN
        D9GMIT = sgng2*EXP(alg2)
        D9GMIT = D9GMIT + EXP(algs)
      END IF
      RETURN
    END IF
  END IF
  !
END FUNCTION D9GMIT