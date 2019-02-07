!*==DT3.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DT3
REAL(8) FUNCTION DT3(X)
  IMPLICIT NONE
  !*--DT35
  !***BEGIN PROLOGUE  DT3
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  DF3S
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  DT3
  REAL(8) :: a, b, DF3S, X, x1, y
  !***FIRST EXECUTABLE STATEMENT  DT3
  a = 0.0D+00
  b = 0.5D+01
  x1 = X + 0.1D+01
  y = (b-a)/x1 + a
  DT3 = (b-a)*DF3S(y)/x1/x1
END FUNCTION DT3
