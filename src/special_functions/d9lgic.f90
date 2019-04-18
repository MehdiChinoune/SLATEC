!** D9LGIC
REAL(8) FUNCTION D9LGIC(A,X,Alx)
  !>
  !***
  !  Compute the log complementary incomplete Gamma function
  !            for large X and for A .LE. X.
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
  ! and for A .LE. X.
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
  USE service, ONLY : XERMSG, D1MACH
  INTEGER k
  REAL(8) :: A, X, Alx, fk, p, r, s, t, xma, xpa
  REAL(8) :: eps = 0.D0
  !* FIRST EXECUTABLE STATEMENT  D9LGIC
  IF ( eps==0.D0 ) eps = 0.5D0*D1MACH(3)
  !
  xpa = X + 1.0D0 - A
  xma = X - 1.D0 - A
  !
  r = 0.D0
  p = 1.D0
  s = p
  DO k = 1, 300
    fk = k
    t = fk*(A-fk)*(1.D0+r)
    r = -t/((xma+2.D0*fk)*(xpa+2.D0*fk)+t)
    p = r*p
    s = s + p
    IF ( ABS(p)<eps*s ) GOTO 100
  END DO
  CALL XERMSG('SLATEC','D9LGIC',&
    'NO CONVERGENCE IN 300 TERMS OF CONTINUED FRACTION',1,2)
  !
  100  D9LGIC = A*Alx - X + LOG(s/xpa)
  !
END FUNCTION D9LGIC
