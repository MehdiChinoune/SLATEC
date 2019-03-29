!** CTAN
COMPLEX FUNCTION CTAN(Z)
  IMPLICIT NONE
  !>
  !***
  !  Compute the complex tangent.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C4A
  !***
  ! **Type:**      COMPLEX (CTAN-C)
  !***
  ! **Keywords:**  ELEMENTARY FUNCTIONS, FNLIB, TANGENT, TRIGONOMETRIC
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! CTAN(Z) calculates the complex trigonometric tangent of complex
  ! argument Z.  Z is in units of radians.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  R1MACH, XERCLR, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)

  REAL den, R1MACH, sn2x, x2, y2
  COMPLEX Z
  REAL :: sqeps = 0.
  !* FIRST EXECUTABLE STATEMENT  CTAN
  IF ( sqeps==0. ) sqeps = SQRT(R1MACH(4))
  !
  x2 = 2.0*REAL(Z)
  y2 = 2.0*AIMAG(Z)
  !
  sn2x = SIN(x2)
  CALL XERCLR
  !
  den = COS(x2) + COSH(y2)
  IF ( den==0. ) CALL XERMSG('SLATEC','CTAN',&
    'TAN IS SINGULAR FOR INPUT Z (X IS PI/2 OR 3*PI/2 AND Y IS 0)',2,2)
  !
  IF ( ABS(den)<=MAX(ABS(x2),1.)*sqeps ) THEN
    CALL XERCLR
    CALL XERMSG('SLATEC','CTAN',&
      'ANSWER LT HALF PRECISION, ABS(X) TOO BIG OR X TOO NEAR PI/2 OR 3*PI/2',1,1)
  ENDIF
  !
  CTAN = CMPLX(sn2x/den,SINH(y2)/den)
  !
END FUNCTION CTAN
