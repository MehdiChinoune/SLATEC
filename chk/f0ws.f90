!*==F0WS.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK F0WS
REAL FUNCTION F0WS(X)
  IMPLICIT NONE
  !*--F0WS5
  !***BEGIN PROLOGUE  F0WS
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  F0WS
  REAL X
  !***FIRST EXECUTABLE STATEMENT  F0WS
  F0WS = SIN(10.0*X)
END FUNCTION F0WS
