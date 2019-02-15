!DECK CGAMR
COMPLEX FUNCTION CGAMR(Z)
  IMPLICIT NONE
  INTEGER irold
  REAL x
  !***BEGIN PROLOGUE  CGAMR
  !***PURPOSE  Compute the reciprocal of the Gamma function.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C7A
  !***TYPE      COMPLEX (GAMR-S, DGAMR-D, CGAMR-C)
  !***KEYWORDS  FNLIB, RECIPROCAL GAMMA FUNCTION, SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! CGAMR(Z) calculates the reciprocal gamma function for COMPLEX
  ! argument Z.  This is a preliminary version that is not accurate.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  CLNGAM, XERCLR, XGETF, XSETF
  !***REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   861211  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  CGAMR
  COMPLEX Z, CLNGAM
  !***FIRST EXECUTABLE STATEMENT  CGAMR
  CGAMR = (0.0,0.0)
  x = REAL(Z)
  IF ( x<=0.0.AND.AINT(x)==x.AND.AIMAG(Z)==0.0 ) RETURN
  !
  CALL XGETF(irold)
  CALL XSETF(1)
  CGAMR = CLNGAM(Z)
  CALL XERCLR
  CALL XSETF(irold)
  CGAMR = EXP(-CGAMR)
  !
END FUNCTION CGAMR
