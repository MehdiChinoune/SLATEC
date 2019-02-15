!DECK CSINH
COMPLEX FUNCTION CSINH(Z)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  CSINH
  !***PURPOSE  Compute the complex hyperbolic sine.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C4C
  !***TYPE      COMPLEX (CSINH-C)
  !***KEYWORDS  ELEMENTARY FUNCTIONS, FNLIB, HYPERBOLIC SINE
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! CSINH(Z) calculates the complex hyperbolic sine of complex
  ! argument Z.  Z is in units of radians.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  CSINH
  COMPLEX Z, ci
  SAVE ci
  DATA ci/(0.,1.)/
  !***FIRST EXECUTABLE STATEMENT  CSINH
  CSINH = -ci*SIN(ci*Z)
  !
END FUNCTION CSINH
