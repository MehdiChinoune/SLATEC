!DECK CACOS
COMPLEX FUNCTION CACOS(Z)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  CACOS
  !***PURPOSE  Compute the complex arc cosine.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C4A
  !***TYPE      COMPLEX (CACOS-C)
  !***KEYWORDS  ARC COSINE, ELEMENTARY FUNCTIONS, FNLIB, TRIGONOMETRIC
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! CACOS(Z) calculates the complex trigonometric arc cosine of Z.
  ! The result is in units of radians, and the real part is in the
  ! first or second quadrant.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  CASIN
  !***REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   861211  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  CACOS
  REAL pi2
  COMPLEX Z, CASIN
  SAVE pi2
  DATA pi2/1.57079632679489661923E0/
  !***FIRST EXECUTABLE STATEMENT  CACOS
  CACOS = pi2 - CASIN(Z)
  !
END FUNCTION CACOS
