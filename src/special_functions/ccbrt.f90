!** CCBRT
COMPLEX(SP) ELEMENTAL FUNCTION CCBRT(Z)
  !> Compute the cube root.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C2
  !***
  ! **Type:**      COMPLEX (CBRT-S, DCBRT-D, CCBRT-C)
  !***
  ! **Keywords:**  CUBE ROOT, ELEMENTARY FUNCTIONS, FNLIB, ROOTS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! CCBRT(Z) calculates the complex cube root of Z.  The principal root
  ! for which -PI < arg(Z) <= +PI is returned.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  CARG, CBRT

  !* REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)

  COMPLEX(SP), INTENT(IN) :: Z
  REAL(SP) :: r, theta
  !* FIRST EXECUTABLE STATEMENT  CCBRT
  theta = CARG(Z)/3._SP
  r = CBRT(ABS(Z))
  !
  CCBRT = CMPLX(r*COS(theta),r*SIN(theta),SP)
  !
END FUNCTION CCBRT