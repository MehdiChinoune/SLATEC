!*==R9LGIC.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK R9LGIC
FUNCTION R9LGIC(A,X,Alx)
  IMPLICIT NONE
  !*--R9LGIC5
  !*** Start of declarations inserted by SPAG
  REAL A, Alx, eps, fk, p, r, R1MACH, R9LGIC, s, t, X, xma, xpa
  INTEGER k
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  R9LGIC
  !***SUBSIDIARY
  !***PURPOSE  Compute the log complementary incomplete Gamma function
  !            for large X and for A .LE. X.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C7E
  !***TYPE      SINGLE PRECISION (R9LGIC-S, D9LGIC-D)
  !***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB, LARGE X,
  !             LOGARITHM, SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! Compute the log complementary incomplete gamma function for large X
  ! and for A .LE. X.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  R1MACH, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900720  Routine changed from user-callable to subsidiary.  (WRB)
  !***END PROLOGUE  R9LGIC
  SAVE eps
  DATA eps/0.0/
  !***FIRST EXECUTABLE STATEMENT  R9LGIC
  IF ( eps==0.0 ) eps = 0.5*R1MACH(3)
  !
  xpa = X + 1.0 - A
  xma = X - 1.0 - A
  !
  r = 0.0
  p = 1.0
  s = p
  DO k = 1, 200
    fk = k
    t = fk*(A-fk)*(1.0+r)
    r = -t/((xma+2.0*fk)*(xpa+2.0*fk)+t)
    p = r*p
    s = s + p
    IF ( ABS(p)<eps*s ) GOTO 100
  ENDDO
  CALL XERMSG('SLATEC','R9LGIC',&
    'NO CONVERGENCE IN 200 TERMS OF CONTINUED FRACTION',1,2)
  !
  100  R9LGIC = A*Alx - X + LOG(s/xpa)
  !
END FUNCTION R9LGIC
