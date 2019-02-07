!*==DT4.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DT4
REAL(8) FUNCTION DT4(X)
  IMPLICIT NONE
  !*--DT45
  !***BEGIN PROLOGUE  DT4
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  DF4S
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  DT4
  REAL(8) :: a, b, DF4S, X, x1, y
  !***FIRST EXECUTABLE STATEMENT  DT4
  a = 0.0D+00
  b = 0.1D+01
  x1 = X + 0.1D+01
  y = (b-a)/x1 + a
  DT4 = (b-a)*DF4S(y)/x1/x1
END FUNCTION DT4
