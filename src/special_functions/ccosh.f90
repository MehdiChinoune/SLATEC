!** CCOSH
COMPLEX FUNCTION CCOSH(Z)
  IMPLICIT NONE
  !>
  !***
  !  Compute the complex hyperbolic cosine.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C4C
  !***
  ! **Type:**      COMPLEX (CCOSH-C)
  !***
  ! **Keywords:**  ELEMENTARY FUNCTIONS, FNLIB, HYPERBOLIC COSINE
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! CCOSH(Z) calculates the complex hyperbolic cosine of Z.
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
  COMPLEX, PARAMETER :: ci  = (0.,1.)
  !* FIRST EXECUTABLE STATEMENT  CCOSH
  CCOSH = COS(ci*Z)
  !
END FUNCTION CCOSH
