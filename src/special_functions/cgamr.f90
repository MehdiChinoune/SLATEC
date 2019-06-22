!** CGAMR
COMPLEX(SP) FUNCTION CGAMR(Z)
  !> Compute the reciprocal of the Gamma function.
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
  USE service, ONLY : control_xer, num_xer
  INTEGER :: irold
  REAL(SP) :: x
  COMPLEX(SP) :: Z
  !* FIRST EXECUTABLE STATEMENT  CGAMR
  CGAMR = (0._SP,0._SP)
  x = REAL(Z)
  IF( x<=0._SP .AND. AINT(x)==x .AND. AIMAG(Z)==0._SP ) RETURN
  !
  irold = control_xer
  control_xer = 1
  CGAMR = CLNGAM(Z)
  num_xer = 0
  control_xer = irold
  CGAMR = EXP(-CGAMR)
  !
END FUNCTION CGAMR
