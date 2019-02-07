!*==F1S.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK F1S
REAL FUNCTION F1S(X)
  IMPLICIT NONE
  !*--F1S5
  !***BEGIN PROLOGUE  F1S
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  F1S
  REAL X
  !***FIRST EXECUTABLE STATEMENT  F1S
  F1S = 0.2E+01/(0.2E+01+SIN(0.314159E+02*X))
END FUNCTION F1S
