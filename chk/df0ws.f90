!*==DF0WS.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DF0WS
DOUBLE PRECISION FUNCTION DF0WS(X)
  IMPLICIT NONE
  !*--DF0WS5
  !***BEGIN PROLOGUE  DF0WS
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  DF0WS
  DOUBLE PRECISION X
  !***FIRST EXECUTABLE STATEMENT  DF0WS
  DF0WS = SIN(0.1D+02*X)
END FUNCTION DF0WS
