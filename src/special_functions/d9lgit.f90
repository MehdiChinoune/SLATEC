!** D9LGIT
REAL(DP) ELEMENTAL FUNCTION D9LGIT(A,X,Algap1)
  !> Compute the logarithm of Tricomi's incomplete Gamma function with Perron's
  !  continued fraction for large X and A >= X.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C7E
  !***
  ! **Type:**      DOUBLE PRECISION (R9LGIT-S, D9LGIT-D)
  !***
  ! **Keywords:**  FNLIB, INCOMPLETE GAMMA FUNCTION, LOGARITHM,
  !             PERRON'S CONTINUED FRACTION, SPECIAL FUNCTIONS, TRICOMI
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! Compute the log of Tricomi's incomplete gamma function with Perron's
  ! continued fraction for large X and for A >= X.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  D1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900720  Routine changed from user-callable to subsidiary.  (WRB)
  USE service, ONLY : D1MACH
  REAL(DP), INTENT(IN) :: A, X, Algap1
  INTEGER :: k
  REAL(DP) :: ax, a1x, fk, hstar, p, r, s, t
  REAL(DP), PARAMETER :: eps = 0.5_DP*D1MACH(3), sqeps = SQRT(D1MACH(4))
  !* FIRST EXECUTABLE STATEMENT  D9LGIT
  !
  IF( X<=0._DP .OR. A<X ) ERROR STOP 'D9LGIT : X SHOULD BE > 0.0 AND <= A'
  !
  ax = A + X
  a1x = ax + 1._DP
  r = 0._DP
  p = 1._DP
  s = p
  DO k = 1, 200
    fk = k
    t = (A+fk)*X*(1._DP+r)
    r = t/((ax+fk)*(a1x+fk)-t)
    p = r*p
    s = s + p
    IF( ABS(p)<eps*s ) GOTO 100
  END DO
  ERROR STOP 'D9LGIT : NO CONVERGENCE IN 200 TERMS OF CONTINUED FRACTION'
  !
  100  hstar = 1._DP - X*s/a1x
  !IF( hstar<sqeps ) CALL XERMSG('D9LGIT','RESULT LESS THAN HALF PRECISION',1,1)
  !
  D9LGIT = -X - Algap1 - LOG(hstar)
  !
END FUNCTION D9LGIT