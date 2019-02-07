!*==T1.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK T1
REAL FUNCTION T1(X)
  IMPLICIT NONE
  !*--T15
  !***BEGIN PROLOGUE  T1
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  F1S
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  T1
  REAL a , b , F1S , X , x1 , y
  !***FIRST EXECUTABLE STATEMENT  T1
  a = 0.0E+00
  b = 0.1E+01
  x1 = X + 0.1E+01
  y = (b-a)/x1 + a
  T1 = (b-a)*F1S(y)/x1/x1
END FUNCTION T1
