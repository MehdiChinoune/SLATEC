!*==DF1P.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DF1P
DOUBLE PRECISION FUNCTION DF1P(X)
  IMPLICIT NONE
  !*--DF1P5
  !***BEGIN PROLOGUE  DF1P
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  DF1P
  DOUBLE PRECISION alfa1 , alfa2 , p1 , p2 , X , d1 , d2
  !***FIRST EXECUTABLE STATEMENT  DF1P
  !  P1 = 1/7, P2 = 2/3
  DATA p1/0.1428571428571428D+00/
  DATA p2/0.6666666666666667D+00/
  alfa1 = -0.25D0
  alfa2 = -0.5D0
  d1 = ABS(X-p1)
  d2 = ABS(X-p2)
  DF1P = 0.0D+00
  IF ( d1/=0.0D+00.AND.d2/=0.0D+00 ) DF1P = d1**alfa1 + d2**alfa2
END FUNCTION DF1P
