!** CGAMMA
COMPLEX(SP) FUNCTION CGAMMA(Z)
  !> Compute the complete Gamma function.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C7A
  !***
  ! **Type:**      COMPLEX (GAMMA-S, DGAMMA-D, CGAMMA-C)
  !***
  ! **Keywords:**  COMPLETE GAMMA FUNCTION, FNLIB, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! CGAMMA(Z) calculates the complete gamma function for COMPLEX
  ! argument Z.  This is a preliminary version that is portable
  ! but not accurate.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  CLNGAM

  !* REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   861211  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)

  COMPLEX(SP) :: Z
  !* FIRST EXECUTABLE STATEMENT  CGAMMA
  CGAMMA = EXP(CLNGAM(Z))
  !
END FUNCTION CGAMMA
