!** CARG
REAL FUNCTION CARG(Z)
  !>
  !  Compute the argument of a complex number.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  A4A
  !***
  ! **Type:**      COMPLEX (CARG-C)
  !***
  ! **Keywords:**  ARGUMENT OF A COMPLEX NUMBER, ELEMENTARY FUNCTIONS, FNLIB
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! CARG(Z) calculates the argument of the complex number Z.  Note
  ! that CARG returns a real result.  If Z = X+iY, then CARG is ATAN(Y/X),
  ! except when both X and Y are zero, in which case the result
  ! will be zero.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   861211  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  
  COMPLEX Z
  !* FIRST EXECUTABLE STATEMENT  CARG
  CARG = 0.0
  IF ( REAL(Z)/=0..OR.AIMAG(Z)/=0. ) CARG = ATAN2(AIMAG(Z),REAL(Z))
  !
END FUNCTION CARG
