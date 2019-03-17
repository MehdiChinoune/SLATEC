!DECK CCOT
COMPLEX FUNCTION CCOT(Z)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  CCOT
  !***PURPOSE  Compute the cotangent.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C4A
  !***TYPE      COMPLEX (COT-S, DCOT-D, CCOT-C)
  !***KEYWORDS  COTANGENT, ELEMENTARY FUNCTIONS, FNLIB, TRIGONOMETRIC
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! CCOT(Z) calculates the complex trigonometric cotangent of Z.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  R1MACH, XERCLR, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !***END PROLOGUE  CCOT
  REAL den, R1MACH, sn2x, sqeps, x2, y2
  COMPLEX Z
  SAVE sqeps
  DATA sqeps/0./
  !***FIRST EXECUTABLE STATEMENT  CCOT
  IF ( sqeps==0. ) sqeps = SQRT(R1MACH(4))
  !
  x2 = 2.0*REAL(Z)
  y2 = 2.0*AIMAG(Z)
  !
  sn2x = SIN(x2)
  CALL XERCLR
  !
  den = COSH(y2) - COS(x2)
  IF ( den==0. ) CALL XERMSG('SLATEC','CCOT',&
    'COT IS SINGULAR FOR INPUT Z (X IS 0 OR PI AND Y IS 0)',2,2)
  !
  IF ( ABS(den)<=MAX(ABS(x2),1.)*sqeps ) THEN
    CALL XERCLR
    CALL XERMSG('SLATEC','CCOT',&
      'ANSWER LT HALF PRECISION, ABS(X) TOO BIG OR X TOO NEAR 0 OR PI',1,1)
  ENDIF
  !
  CCOT = CMPLX(sn2x/den,-SINH(y2)/den)
  !
END FUNCTION CCOT