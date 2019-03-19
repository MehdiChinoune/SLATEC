!** CASINH
COMPLEX FUNCTION CASINH(Z)
  IMPLICIT NONE
  !>
  !***
  !  Compute the arc hyperbolic sine.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C4C
  !***
  ! **Type:**      COMPLEX (ASINH-S, DASINH-D, CASINH-C)
  !***
  ! **Keywords:**  ARC HYPERBOLIC SINE, ASINH, ELEMENTARY FUNCTIONS, FNLIB,
  !             INVERSE HYPERBOLIC SINE
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! CASINH(Z) calculates the complex arc hyperbolic sine of Z.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  CASIN

  !* REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   861211  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  
  COMPLEX Z, ci, CASIN
  SAVE ci
  DATA ci/(0.,1.)/
  !* FIRST EXECUTABLE STATEMENT  CASINH
  CASINH = -ci*CASIN(ci*Z)
  !
END FUNCTION CASINH
