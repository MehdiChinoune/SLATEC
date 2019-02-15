!DECK F1P
REAL FUNCTION F1P(X)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  F1P
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  F1P
  REAL alfa1, alfa2, p1, p2, X, d1, d2
  !  P1 = 1/7, P2 = 2/3
  DATA p1/0.1428571428571428E+00/
  DATA p2/0.6666666666666667E+00/
  !***FIRST EXECUTABLE STATEMENT  F1P
  alfa1 = -0.25E0
  alfa2 = -0.5E0
  d1 = ABS(X-p1)
  d2 = ABS(X-p2)
  F1P = 0.0E+00
  IF ( d1/=0.0E+00.AND.d2/=0.0E+00 ) F1P = d1**alfa1 + d2**alfa2
END FUNCTION F1P
