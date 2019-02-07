!*==DT1.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DT1
REAL(8) FUNCTION DT1(X)
  IMPLICIT NONE
  !*--DT15
  !***BEGIN PROLOGUE  DT1
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  DF1S
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  DT1
  REAL(8) :: a , b , DF1S , X , x1 , y
  !***FIRST EXECUTABLE STATEMENT  DT1
  a = 0.0D+00
  b = 0.1D+01
  x1 = X + 0.1D+01
  y = (b-a)/x1 + a
  DT1 = (b-a)*DF1S(y)/x1/x1
END FUNCTION DT1
