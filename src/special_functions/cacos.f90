!** CACOS
COMPLEX FUNCTION CACOS(Z)
  IMPLICIT NONE
  !>
  !***
  !  Compute the complex arc cosine.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C4A
  !***
  ! **Type:**      COMPLEX (CACOS-C)
  !***
  ! **Keywords:**  ARC COSINE, ELEMENTARY FUNCTIONS, FNLIB, TRIGONOMETRIC
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! CACOS(Z) calculates the complex trigonometric arc cosine of Z.
  ! The result is in units of radians, and the real part is in the
  ! first or second quadrant.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  CASIN

  !* REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   861211  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  
  REAL pi2
  COMPLEX Z, CASIN
  SAVE pi2
  DATA pi2/1.57079632679489661923E0/
  !* FIRST EXECUTABLE STATEMENT  CACOS
  CACOS = pi2 - CASIN(Z)
  !
END FUNCTION CACOS
