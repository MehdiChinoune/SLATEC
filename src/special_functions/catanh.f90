!** CATANH
COMPLEX FUNCTION CATANH(Z)
  IMPLICIT NONE
  !>
  !***
  !  Compute the arc hyperbolic tangent.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C4C
  !***
  ! **Type:**      COMPLEX (ATANH-S, DATANH-D, CATANH-C)
  !***
  ! **Keywords:**  ARC HYPERBOLIC TANGENT, ATANH, ELEMENTARY FUNCTIONS,
  !             FNLIB, INVERSE HYPERBOLIC TANGENT
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! CATANH(Z) calculates the complex arc hyperbolic tangent of Z.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  CATAN

  !* REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   861211  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  
  COMPLEX Z, ci, CATAN
  SAVE ci
  DATA ci/(0.,1.)/
  !* FIRST EXECUTABLE STATEMENT  CATANH
  CATANH = -ci*CATAN(ci*Z)
  !
END FUNCTION CATANH
