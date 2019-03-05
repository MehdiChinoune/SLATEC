!DECK ALNGAM
FUNCTION ALNGAM(X)
  IMPLICIT NONE
  REAL ALNGAM, dxrel, GAMMA, pi, R1MACH, R9LGMC, sinpiy, sq2pil, &
    sqpi2l, X, xmax, y
  !***BEGIN PROLOGUE  ALNGAM
  !***PURPOSE  Compute the logarithm of the absolute value of the Gamma
  !            function.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C7A
  !***TYPE      SINGLE PRECISION (ALNGAM-S, DLNGAM-D, CLNGAM-C)
  !***KEYWORDS  ABSOLUTE VALUE, COMPLETE GAMMA FUNCTION, FNLIB, LOGARITHM,
  !             SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! ALNGAM(X) computes the logarithm of the absolute value of the
  ! gamma function at X.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  GAMMA, R1MACH, R9LGMC, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   900727  Added EXTERNAL statement.  (WRB)
  !***END PROLOGUE  ALNGAM
  LOGICAL first
  EXTERNAL GAMMA
  SAVE sq2pil, sqpi2l, pi, xmax, dxrel, first
  DATA sq2pil/0.91893853320467274E0/
  DATA sqpi2l/0.22579135264472743E0/
  DATA pi/3.14159265358979324E0/
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  ALNGAM
  IF ( first ) THEN
    xmax = R1MACH(2)/LOG(R1MACH(2))
    dxrel = SQRT(R1MACH(4))
  ENDIF
  first = .FALSE.
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
    ENDIF
    !
    sinpiy = ABS(SIN(pi*y))
    IF ( sinpiy==0. ) CALL XERMSG('SLATEC','ALNGAM',&
      'X IS A NEGATIVE INTEGER',3,2)
    !
    IF ( ABS((X-AINT(X-0.5))/X)<dxrel ) CALL XERMSG('SLATEC','ALNGAM',&
      'ANSWER LT HALF PRECISION BECAUSE X TOO NEAR NEGATIVE INTEGER',&
      1,1)
    !
    ALNGAM = sqpi2l + (X-0.5)*LOG(y) - X - LOG(sinpiy) - R9LGMC(y)
    RETURN
  ENDIF
  !
  ! LOG (ABS (GAMMA(X))) FOR  ABS(X) .LE. 10.0
  !
  ALNGAM = LOG(ABS(GAMMA(X)))
  RETURN
END FUNCTION ALNGAM
