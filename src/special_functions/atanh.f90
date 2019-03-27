!** ATANH
REAL FUNCTION ATANH(X)
  IMPLICIT NONE
  !>
  !***
  !  Compute the arc hyperbolic tangent.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C4C
  !***
  ! **Type:**      SINGLE PRECISION (ATANH-S, DATANH-D, CATANH-C)
  !***
  ! **Keywords:**  ARC HYPERBOLIC TANGENT, ATANH, ELEMENTARY FUNCTIONS,
  !             FNLIB, INVERSE HYPERBOLIC TANGENT
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! ATANH(X) computes the arc hyperbolic tangent of X.
  !
  ! Series for ATNH       on the interval  0.          to  2.50000D-01
  !                                        with weighted error   6.70E-18
  !                                         log weighted error  17.17
  !                               significant figures required  16.01
  !                                    decimal places required  17.76
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  CSEVL, INITS, R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)

  REAL CSEVL, dxrel, R1MACH, sqeps, X, y
  INTEGER INITS, nterms
  SAVE nterms, dxrel, sqeps
  REAL, PARAMETER :: atnhcs(15) = [ .094395102393195492E0, .049198437055786159E0, &
    .002102593522455432E0, .000107355444977611E0, .000005978267249293E0, &
    .000000350506203088E0, .000000021263743437E0, .000000001321694535E0, &
    .000000000083658755E0, .000000000005370503E0, .000000000000348665E0, &
    .000000000000022845E0, .000000000000001508E0, .000000000000000100E0, &
    .000000000000000006E0 ]
  LOGICAL :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  ATANH
  IF ( first ) THEN
    nterms = INITS(atnhcs,15,0.1*R1MACH(3))
    dxrel = SQRT(R1MACH(4))
    sqeps = SQRT(3.0*R1MACH(3))
    first = .FALSE.
  ENDIF
  !
  y = ABS(X)
  IF ( y>=1.0 ) CALL XERMSG('SLATEC','ATANH','ABS(X) GE 1',2,2)
  !
  IF ( 1.0-y<dxrel ) CALL XERMSG('SLATEC','ATANH',&
    'ANSWER LT HALF PRECISION BECAUSE ABS(X) TOO NEAR 1',1,1)
  !
  ATANH = X
  IF ( y>sqeps.AND.y<=0.5 ) ATANH = X*(1.0+CSEVL(8.*X*X-1.,atnhcs,nterms))
  IF ( y>0.5 ) ATANH = 0.5*LOG((1.0+X)/(1.0-X))
  !
END FUNCTION ATANH
