!*==DT5.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DT5
REAL(8) FUNCTION DT5(X)
  IMPLICIT NONE
  !*--DT55
  !***BEGIN PROLOGUE  DT5
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  DF5S
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  DT5
  REAL(8) :: a, b, DF5S, X, x1, y
  !***FIRST EXECUTABLE STATEMENT  DT5
  a = 0.0D+00
  b = 0.1D+01
  x1 = X + 0.1D+01
  y = (b-a)/x1 + a
  DT5 = (b-a)*DF5S(y)/x1/x1
END FUNCTION DT5
