!** CGAMR
COMPLEX FUNCTION CGAMR(Z)
  !>
  !***
  !  Compute the reciprocal of the Gamma function.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C7A
  !***
  ! **Type:**      COMPLEX (GAMR-S, DGAMR-D, CGAMR-C)
  !***
  ! **Keywords:**  FNLIB, RECIPROCAL GAMMA FUNCTION, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! CGAMR(Z) calculates the reciprocal gamma function for COMPLEX
  ! argument Z.  This is a preliminary version that is not accurate.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  CLNGAM, XERCLR, XGETF, XSETF

  !* REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   861211  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)

  INTEGER irold
  REAL x
  COMPLEX Z
  !* FIRST EXECUTABLE STATEMENT  CGAMR
  CGAMR = (0.0,0.0)
  x = REAL(Z)
  IF ( x<=0.0.AND.AINT(x)==x.AND.AIMAG(Z)==0.0 ) RETURN
  !
  CALL XGETF(irold)
  CALL XSETF(1)
  CGAMR = CLNGAM(Z)
  CALL XERCLR
  CALL XSETF(irold)
  CGAMR = EXP(-CGAMR)
  !
END FUNCTION CGAMR
