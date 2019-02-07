!*==F2S.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK F2S
REAL FUNCTION F2S(X)
  IMPLICIT NONE
  !*--F2S5
  !***BEGIN PROLOGUE  F2S
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  F2S
  REAL X
  !***FIRST EXECUTABLE STATEMENT  F2S
  F2S = 100.0
  IF ( X/=0.0 ) F2S = SIN(0.314159E+03*X)/(0.314159E+01*X)
END FUNCTION F2S
