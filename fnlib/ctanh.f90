!DECK CTANH
COMPLEX FUNCTION CTANH(Z)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  CTANH
  !***PURPOSE  Compute the complex hyperbolic tangent.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C4C
  !***TYPE      COMPLEX (CTANH-C)
  !***KEYWORDS  ELEMENTARY FUNCTIONS, FNLIB, HYPERBOLIC TANGENT
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! CTANH(Z) calculates the complex hyperbolic tangent of complex
  ! argument Z.  Z is in units of radians.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  CTAN
  !***REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   861211  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  CTANH
  COMPLEX Z, ci, CTAN
  SAVE ci
  DATA ci/(0.,1.)/
  !***FIRST EXECUTABLE STATEMENT  CTANH
  CTANH = -ci*CTAN(ci*Z)
  !
END FUNCTION CTANH
