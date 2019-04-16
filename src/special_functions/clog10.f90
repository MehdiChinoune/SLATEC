!** CLOG10
COMPLEX FUNCTION CLOG10(Z)
  !>
  !***
  !  Compute the principal value of the complex base 10
  !            logarithm.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C4B
  !***
  ! **Type:**      COMPLEX (CLOG10-C)
  !***
  ! **Keywords:**  BASE TEN LOGARITHM, ELEMENTARY FUNCTIONS, FNLIB
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! CLOG10(Z) calculates the principal value of the complex common
  ! or base 10 logarithm of Z for -PI .LT. arg(Z) .LE. +PI.
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
  REAL, PARAMETER :: aloge = 0.43429448190325182765E0
  !* FIRST EXECUTABLE STATEMENT  CLOG10
  CLOG10 = aloge*LOG(Z)
  !
END FUNCTION CLOG10
