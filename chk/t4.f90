!DECK T4
REAL FUNCTION T4(X)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  T4
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  F4S
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  T4
  REAL a, b, F4S, X, x1, y
  !***FIRST EXECUTABLE STATEMENT  T4
  a = 0.0E+00
  b = 0.1E+01
  x1 = X + 0.1E+01
  y = (b-a)/x1 + a
  T4 = (b-a)*F4S(y)/x1/x1
END FUNCTION T4
