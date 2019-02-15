!DECK D9LGIT
REAL(8) FUNCTION D9LGIT(A,X,Algap1)
  IMPLICIT NONE
  INTEGER k
  !***BEGIN PROLOGUE  D9LGIT
  !***SUBSIDIARY
  !***PURPOSE  Compute the logarithm of Tricomi's incomplete Gamma
  !            function with Perron's continued fraction for large X and
  !            A .GE. X.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C7E
  !***TYPE      DOUBLE PRECISION (R9LGIT-S, D9LGIT-D)
  !***KEYWORDS  FNLIB, INCOMPLETE GAMMA FUNCTION, LOGARITHM,
  !             PERRON'S CONTINUED FRACTION, SPECIAL FUNCTIONS, TRICOMI
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! Compute the log of Tricomi's incomplete gamma function with Perron's
  ! continued fraction for large X and for A .GE. X.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  D1MACH, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900720  Routine changed from user-callable to subsidiary.  (WRB)
  !***END PROLOGUE  D9LGIT
  REAL(8) :: A, X, Algap1, ax, a1x, eps, fk, hstar, p, r, &
    s, sqeps, t, D1MACH
  LOGICAL first
  SAVE eps, sqeps, first
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  D9LGIT
  IF ( first ) THEN
    eps = 0.5D0*D1MACH(3)
    sqeps = SQRT(D1MACH(4))
  ENDIF
  first = .FALSE.
  !
  IF ( X<=0.D0.OR.A<X ) CALL XERMSG('SLATEC','D9LGIT',&
    'X SHOULD BE GT 0.0 AND LE A',2,2)
  !
  ax = A + X
  a1x = ax + 1.0D0
  r = 0.D0
  p = 1.D0
  s = p
  DO k = 1, 200
    fk = k
    t = (A+fk)*X*(1.D0+r)
    r = t/((ax+fk)*(a1x+fk)-t)
    p = r*p
    s = s + p
    IF ( ABS(p)<eps*s ) GOTO 100
  ENDDO
  CALL XERMSG('SLATEC','D9LGIT',&
    'NO CONVERGENCE IN 200 TERMS OF CONTINUED FRACTION',3,2)
  !
  100  hstar = 1.0D0 - X*s/a1x
  IF ( hstar<sqeps ) CALL XERMSG('SLATEC','D9LGIT',&
    'RESULT LESS THAN HALF PRECISION',1,1)
  !
  D9LGIT = -X - Algap1 - LOG(hstar)
  !
END FUNCTION D9LGIT
