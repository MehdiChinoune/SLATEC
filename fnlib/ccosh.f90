!DECK CCOSH
COMPLEX FUNCTION CCOSH(Z)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  CCOSH
  !***PURPOSE  Compute the complex hyperbolic cosine.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C4C
  !***TYPE      COMPLEX (CCOSH-C)
  !***KEYWORDS  ELEMENTARY FUNCTIONS, FNLIB, HYPERBOLIC COSINE
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! CCOSH(Z) calculates the complex hyperbolic cosine of Z.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  CCOSH
  COMPLEX Z, ci
  SAVE ci
  DATA ci/(0.,1.)/
  !***FIRST EXECUTABLE STATEMENT  CCOSH
  CCOSH = COS(ci*Z)
  !
END FUNCTION CCOSH
