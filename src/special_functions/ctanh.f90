!** CTANH
COMPLEX FUNCTION CTANH(Z)
  IMPLICIT NONE
  !>
  !***
  !  Compute the complex hyperbolic tangent.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C4C
  !***
  ! **Type:**      COMPLEX (CTANH-C)
  !***
  ! **Keywords:**  ELEMENTARY FUNCTIONS, FNLIB, HYPERBOLIC TANGENT
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! CTANH(Z) calculates the complex hyperbolic tangent of complex
  ! argument Z.  Z is in units of radians.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  CTAN

  !* REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   861211  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)

  COMPLEX Z, CTAN
  COMPLEX, PARAMETER :: ci = (0.,1.)
  !* FIRST EXECUTABLE STATEMENT  CTANH
  CTANH = -ci*CTAN(ci*Z)
  !
END FUNCTION CTANH
