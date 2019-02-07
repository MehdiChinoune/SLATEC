!*==F2P.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK F2P
REAL FUNCTION F2P(X)
  IMPLICIT NONE
  !*--F2P5
  !***BEGIN PROLOGUE  F2P
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  F2P
  REAL X
  !***FIRST EXECUTABLE STATEMENT  F2P
  F2P = SIN(0.314159E+03*X)/(0.314159E+01*X)
END FUNCTION F2P
