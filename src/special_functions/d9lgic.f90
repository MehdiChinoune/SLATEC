!** D9LGIC
REAL(DP) ELEMENTAL FUNCTION D9LGIC(A,X,Alx)
  !> Compute the log complementary incomplete Gamma function for large X and for A <= X.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C7E
  !***
  ! **Type:**      DOUBLE PRECISION (R9LGIC-S, D9LGIC-D)
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
  ! **Routines called:**  D1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900720  Routine changed from user-callable to subsidiary.  (WRB)
  USE service, ONLY : eps_2_dp
  !
  REAL(DP), INTENT(IN) :: A, X, Alx
  !
  INTEGER :: k
  REAL(DP) :: fk, p, r, s, t, xma, xpa
  REAL(DP), PARAMETER :: eps = 0.5_DP*eps_2_dp
  !* FIRST EXECUTABLE STATEMENT  D9LGIC
  !
  xpa = X + 1._DP - A
  xma = X - 1._DP - A
  !
  r = 0._DP
  p = 1._DP
  s = p
  DO k = 1, 300
    fk = k
    t = fk*(A-fk)*(1._DP+r)
    r = -t/((xma+2._DP*fk)*(xpa+2._DP*fk)+t)
    p = r*p
    s = s + p
    IF( ABS(p)<eps*s ) GOTO 100
  END DO
  ERROR STOP 'D9LGIC : NO CONVERGENCE IN 300 TERMS OF CONTINUED FRACTION'
  !
  100  D9LGIC = A*Alx - X + LOG(s/xpa)
  !
END FUNCTION D9LGIC