!** R9LGIT
REAL(SP) ELEMENTAL FUNCTION R9LGIT(A,X,Algap1)
  !> Compute the logarithm of Tricomi's incomplete Gamma
  !  function with Perron's continued fraction for large X and A >= X.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C7E
  !***
  ! **Type:**      SINGLE PRECISION (R9LGIT-S, D9LGIT-D)
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
  ! **Routines called:**  R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900720  Routine changed from user-callable to subsidiary.  (WRB)
  USE service, ONLY : eps_2_sp, eps_sp
  !
  REAL(SP), INTENT(IN) :: A, Algap1, X
  !
  INTEGER :: k
  REAL(SP) :: a1x, ax, fk, hstar, p, r, s, t
  REAL(SP), PARAMETER :: eps = 0.5_SP*eps_2_sp, sqeps = SQRT(eps_sp)
  !* FIRST EXECUTABLE STATEMENT  R9LGIT
  !
  IF( X<=0._SP .OR. A<X ) ERROR STOP 'R9LGIT : X SHOULD BE > 0.0 AND <= A'
  !
  ax = A + X
  a1x = ax + 1._SP
  r = 0._SP
  p = 1._SP
  s = p
  DO k = 1, 200
    fk = k
    t = (A+fk)*X*(1._SP+r)
    r = t/((ax+fk)*(a1x+fk)-t)
    p = r*p
    s = s + p
    IF( ABS(p)<eps*s ) GOTO 100
  END DO
  ERROR STOP 'R9LGIT : NO CONVERGENCE IN 200 TERMS OF CONTINUED FRACTION'
  !
  100  hstar = 1._SP - X*s/a1x
  ! IF( hstar<sqeps ) CALL XERMSG('R9LGIT','RESULT LESS THAN HALF PRECISION',1,1)
  !
  R9LGIT = -X - Algap1 - LOG(hstar)

END FUNCTION R9LGIT