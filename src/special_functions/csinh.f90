!** CSINH
COMPLEX FUNCTION CSINH(Z)
  IMPLICIT NONE
  !>
  !***
  !  Compute the complex hyperbolic sine.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C4C
  !***
  ! **Type:**      COMPLEX (CSINH-C)
  !***
  ! **Keywords:**  ELEMENTARY FUNCTIONS, FNLIB, HYPERBOLIC SINE
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! CSINH(Z) calculates the complex hyperbolic sine of complex
  ! argument Z.  Z is in units of radians.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)

  COMPLEX Z
  COMPLEX, PARAMETER :: ci = (0.,1.)
  !* FIRST EXECUTABLE STATEMENT  CSINH
  CSINH = -ci*SIN(ci*Z)
  !
END FUNCTION CSINH
