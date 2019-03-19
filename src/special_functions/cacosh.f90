!** CACOSH
COMPLEX FUNCTION CACOSH(Z)
  IMPLICIT NONE
  !>
  !***
  !  Compute the arc hyperbolic cosine.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C4C
  !***
  ! **Type:**      COMPLEX (ACOSH-S, DACOSH-D, CACOSH-C)
  !***
  ! **Keywords:**  ACOSH, ARC HYPERBOLIC COSINE, ELEMENTARY FUNCTIONS, FNLIB,
  !             INVERSE HYPERBOLIC COSINE
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! CACOSH(Z) calculates the complex arc hyperbolic cosine of Z.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  CACOS

  !* REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   861211  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  
  COMPLEX Z, ci, CACOS
  SAVE ci
  DATA ci/(0.,1.)/
  !* FIRST EXECUTABLE STATEMENT  CACOSH
  CACOSH = ci*CACOS(Z)
  !
END FUNCTION CACOSH
