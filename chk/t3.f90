!*==T3.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK T3
REAL FUNCTION T3(X)
  IMPLICIT NONE
  !*--T35
  !***BEGIN PROLOGUE  T3
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  F3S
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  T3
  REAL a , b , F3S , X , x1 , y
  !***FIRST EXECUTABLE STATEMENT  T3
  a = 0.0E+00
  b = 0.5E+01
  x1 = X + 0.1E+01
  y = (b-a)/x1 + a
  T3 = (b-a)*F3S(y)/x1/x1
END FUNCTION T3
