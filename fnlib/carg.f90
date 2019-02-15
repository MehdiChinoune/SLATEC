!DECK CARG
FUNCTION CARG(Z)
  IMPLICIT NONE
  REAL CARG
  !***BEGIN PROLOGUE  CARG
  !***PURPOSE  Compute the argument of a complex number.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  A4A
  !***TYPE      COMPLEX (CARG-C)
  !***KEYWORDS  ARGUMENT OF A COMPLEX NUMBER, ELEMENTARY FUNCTIONS, FNLIB
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! CARG(Z) calculates the argument of the complex number Z.  Note
  ! that CARG returns a real result.  If Z = X+iY, then CARG is ATAN(Y/X),
  ! except when both X and Y are zero, in which case the result
  ! will be zero.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   861211  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  CARG
  COMPLEX Z
  !***FIRST EXECUTABLE STATEMENT  CARG
  CARG = 0.0
  IF ( REAL(Z)/=0..OR.AIMAG(Z)/=0. ) CARG = ATAN2(AIMAG(Z),REAL(Z))
  !
END FUNCTION CARG
