!** CCOT
COMPLEX(SP) FUNCTION CCOT(Z)
  !> Compute the cotangent.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C4A
  !***
  ! **Type:**      COMPLEX (COT-S, DCOT-D, CCOT-C)
  !***
  ! **Keywords:**  COTANGENT, ELEMENTARY FUNCTIONS, FNLIB, TRIGONOMETRIC
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! CCOT(Z) calculates the complex trigonometric cotangent of Z.
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
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  USE service, ONLY : XERMSG, XERCLR, R1MACH
  COMPLEX(SP) :: Z
  REAL(SP) :: den, sn2x, x2, y2
  REAL(SP), PARAMETER :: sqeps = SQRT(R1MACH(4))
  !* FIRST EXECUTABLE STATEMENT  CCOT
  !
  x2 = 2._SP*REAL(Z)
  y2 = 2._SP*AIMAG(Z)
  !
  sn2x = SIN(x2)
  CALL XERCLR
  !
  den = COSH(y2) - COS(x2)
  IF( den==0. ) CALL XERMSG('CCOT',&
    'COT IS SINGULAR FOR INPUT Z (X IS 0 OR PI AND Y IS 0)',2,2)
  !
  IF( ABS(den)<=MAX(ABS(x2),1._SP)*sqeps ) THEN
    CALL XERCLR
    CALL XERMSG('CCOT',&
      'ANSWER LT HALF PRECISION, ABS(X) TOO BIG OR X TOO NEAR 0 OR PI',1,1)
  END IF
  !
  CCOT = CMPLX(sn2x/den,-SINH(y2)/den,SP)
  !
END FUNCTION CCOT
