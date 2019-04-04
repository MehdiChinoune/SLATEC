!** ALNGAM
REAL FUNCTION ALNGAM(X)
  IMPLICIT NONE
  !>
  !***
  !  Compute the logarithm of the absolute value of the Gamma
  !            function.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C7A
  !***
  ! **Type:**      SINGLE PRECISION (ALNGAM-S, DLNGAM-D, CLNGAM-C)
  !***
  ! **Keywords:**  ABSOLUTE VALUE, COMPLETE GAMMA FUNCTION, FNLIB, LOGARITHM,
  !             SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! ALNGAM(X) computes the logarithm of the absolute value of the
  ! gamma function at X.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  GAMMA, R1MACH, R9LGMC, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   900727  Added EXTERNAL statement.  (WRB)

  REAL sinpiy, X, y
  REAL, EXTERNAL :: R1MACH, R9LGMC, GAMMA
  REAL, SAVE :: xmax, dxrel
  REAL, PARAMETER :: sq2pil = 0.91893853320467274E0
  REAL, PARAMETER :: sqpi2l = 0.22579135264472743E0
  REAL, PARAMETER :: pi = 3.14159265358979324E0
  LOGICAL :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  ALNGAM
  IF ( first ) THEN
    xmax = R1MACH(2)/LOG(R1MACH(2))
    dxrel = SQRT(R1MACH(4))
    first = .FALSE.
  END IF
  !
  y = ABS(X)
  IF ( y>10.0 ) THEN
    !
    ! LOG (ABS (GAMMA(X))) FOR ABS(X) .GT. 10.0
    !
    IF ( y>xmax ) CALL XERMSG('SLATEC','ALNGAM',&
      'ABS(X) SO BIG ALNGAM OVERFLOWS',2,2)
    !
    IF ( X>0. ) THEN
      ALNGAM = sq2pil + (X-0.5)*LOG(X) - X + R9LGMC(y)
      RETURN
    END IF
    !
    sinpiy = ABS(SIN(pi*y))
    IF ( sinpiy==0. ) CALL XERMSG('SLATEC','ALNGAM',&
      'X IS A NEGATIVE INTEGER',3,2)
    !
    IF ( ABS((X-AINT(X-0.5))/X)<dxrel ) CALL XERMSG('SLATEC','ALNGAM',&
      'ANSWER LT HALF PRECISION BECAUSE X TOO NEAR NEGATIVE INTEGER',1,1)
    !
    ALNGAM = sqpi2l + (X-0.5)*LOG(y) - X - LOG(sinpiy) - R9LGMC(y)
    RETURN
  END IF
  !
  ! LOG (ABS (GAMMA(X))) FOR  ABS(X) .LE. 10.0
  !
  ALNGAM = LOG(ABS(GAMMA(X)))
  RETURN
END FUNCTION ALNGAM
