!** R9LGIC
REAL(SP) ELEMENTAL FUNCTION R9LGIC(A,X,Alx)
  !> Compute the log complementary incomplete Gamma function for large X and for A <= X.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C7E
  !***
  ! **Type:**      SINGLE PRECISION (R9LGIC-S, D9LGIC-D)
  !***
  ! **Keywords:**  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB, LARGE X,
  !             LOGARITHM, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! Compute the log complementary incomplete gamma function for large X
  ! and for A <= X.
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
  USE service, ONLY : R1MACH
  REAL(SP), INTENT(IN) :: A, Alx, X
  INTEGER :: k
  REAL(SP) :: fk, p, r, s, t, xma, xpa
  REAL(SP), PARAMETER :: eps = 0.5_SP*R1MACH(3)
  !* FIRST EXECUTABLE STATEMENT  R9LGIC
  !
  xpa = X + 1._SP - A
  xma = X - 1._SP - A
  !
  r = 0._SP
  p = 1._SP
  s = p
  DO k = 1, 200
    fk = k
    t = fk*(A-fk)*(1._SP+r)
    r = -t/((xma+2._SP*fk)*(xpa+2._SP*fk)+t)
    p = r*p
    s = s + p
    IF( ABS(p)<eps*s ) GOTO 100
  END DO
  ERROR STOP 'R9LGIC : NO CONVERGENCE IN 200 TERMS OF CONTINUED FRACTION'
  !
  100  R9LGIC = A*Alx - X + LOG(s/xpa)

END FUNCTION R9LGIC