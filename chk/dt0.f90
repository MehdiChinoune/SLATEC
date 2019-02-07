!*==DT0.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DT0
DOUBLE PRECISION FUNCTION DT0(X)
  IMPLICIT NONE
  !*--DT05
  !***BEGIN PROLOGUE  DT0
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  DF0S
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  DT0
  DOUBLE PRECISION a , b , DF0S , X , x1 , y
  !***FIRST EXECUTABLE STATEMENT  DT0
  a = 0.0D+00
  b = 0.1D+01
  x1 = X + 0.1D+01
  y = (b-a)/x1 + a
  DT0 = (b-a)*DF0S(y)/x1/x1
END FUNCTION DT0
