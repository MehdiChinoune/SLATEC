!*==DF0O.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DF0O
DOUBLE PRECISION FUNCTION DF0O(X)
  IMPLICIT NONE
  !*--DF0O5
  !***BEGIN PROLOGUE  DF0O
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  DF0O
  DOUBLE PRECISION X
  !***FIRST EXECUTABLE STATEMENT  DF0O
  DF0O = (0.2D+01*SIN(X))**14
END FUNCTION DF0O
