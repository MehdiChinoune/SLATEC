!DECK T5
REAL FUNCTION T5(X)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  T5
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  F5S
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  T5
  REAL a, b, F5S, X, x1, y
  !***FIRST EXECUTABLE STATEMENT  T5
  a = 0.0E+00
  b = 0.1E+01
  x1 = X + 0.1E+01
  y = (b-a)/x1 + a
  T5 = (b-a)*F5S(y)/x1/x1
END FUNCTION T5
