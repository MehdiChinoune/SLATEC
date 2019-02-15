!DECK T0
REAL FUNCTION T0(X)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  T0
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  F0S
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  T0
  REAL a, b, F0S, X, x1, y
  !***FIRST EXECUTABLE STATEMENT  T0
  a = 0.0E+00
  b = 0.1E+01
  x1 = X + 0.1E+01
  y = (b-a)/x1 + a
  T0 = (b-a)*F0S(y)/x1/x1
END FUNCTION T0
