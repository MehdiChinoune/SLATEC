!*==DT2.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DT2
DOUBLE PRECISION FUNCTION DT2(X)
  IMPLICIT NONE
  !*--DT25
  !***BEGIN PROLOGUE  DT2
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  DF2S
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  DT2
  DOUBLE PRECISION a , b , DF2S , X , x1 , y
  !***FIRST EXECUTABLE STATEMENT  DT2
  a = 0.1D+00
  b = 0.1D+01
  x1 = X + 0.1D+01
  y = (b-a)/x1 + a
  DT2 = (b-a)*DF2S(y)/x1/x1
END FUNCTION DT2
