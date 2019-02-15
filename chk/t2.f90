!DECK T2
REAL FUNCTION T2(X)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  T2
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  F2S
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  T2
  REAL a, b, F2S, X, x1, y
  !***FIRST EXECUTABLE STATEMENT  T2
  a = 0.1E+00
  b = 0.1E+01
  x1 = X + 0.1E+01
  y = (b-a)/x1 + a
  T2 = (b-a)*F2S(y)/x1/x1
END FUNCTION T2
