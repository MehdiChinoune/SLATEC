!*==DLNGAM.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK DLNGAM
REAL(8) FUNCTION DLNGAM(X)
  IMPLICIT NONE
  !*--DLNGAM5
  !***BEGIN PROLOGUE  DLNGAM
  !***PURPOSE  Compute the logarithm of the absolute value of the Gamma
  !            function.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C7A
  !***TYPE      DOUBLE PRECISION (ALNGAM-S, DLNGAM-D, CLNGAM-C)
  !***KEYWORDS  ABSOLUTE VALUE, COMPLETE GAMMA FUNCTION, FNLIB, LOGARITHM,
  !             SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! DLNGAM(X) calculates the double precision logarithm of the
  ! absolute value of the Gamma function for double precision
  ! argument X.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  D1MACH, D9LGMC, DGAMMA, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900727  Added EXTERNAL statement.  (WRB)
  !***END PROLOGUE  DLNGAM
  REAL(8) :: X, dxrel, pi, sinpiy, sqpi2l, sq2pil, xmax, y, &
    DGAMMA, D9LGMC, D1MACH, temp
  LOGICAL first
  EXTERNAL DGAMMA
  SAVE sq2pil, sqpi2l, pi, xmax, dxrel, first
  DATA sq2pil/0.91893853320467274178032973640562D0/
  DATA sqpi2l/ + .225791352644727432363097614947441D+0/
  DATA pi/3.14159265358979323846264338327950D0/
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  DLNGAM
  IF ( first ) THEN
    temp = 1.D0/LOG(D1MACH(2))
    xmax = temp*D1MACH(2)
    dxrel = SQRT(D1MACH(4))
  ENDIF
  first = .FALSE.
  !
  y = ABS(X)
  IF ( y>10.D0 ) THEN
    !
    ! LOG ( ABS (DGAMMA(X)) ) FOR ABS(X) .GT. 10.0
    !
    IF ( y>xmax ) CALL XERMSG('SLATEC','DLNGAM',&
      'ABS(X) SO BIG DLNGAM OVERFLOWS',2,2)
    !
    IF ( X>0.D0 ) DLNGAM = sq2pil + (X-0.5D0)*LOG(X) - X + D9LGMC(y)
    IF ( X>0.D0 ) RETURN
    !
    sinpiy = ABS(SIN(pi*y))
    IF ( sinpiy==0.D0 ) CALL XERMSG('SLATEC','DLNGAM',&
      'X IS A NEGATIVE INTEGER',3,2)
    !
    IF ( ABS((X-AINT(X-0.5D0))/X)<dxrel ) CALL XERMSG('SLATEC','DLNGAM',&
      'ANSWER LT HALF PRECISION BECAUSE X TOO NEAR NEGATIVE INTEGER',1,1)
    !
    DLNGAM = sqpi2l + (X-0.5D0)*LOG(y) - X - LOG(sinpiy) - D9LGMC(y)
    GOTO 99999
  ENDIF
  !
  ! LOG (ABS (DGAMMA(X)) ) FOR ABS(X) .LE. 10.0
  !
  DLNGAM = LOG(ABS(DGAMMA(X)))
  RETURN
  !
  99999 CONTINUE
  END FUNCTION DLNGAM
