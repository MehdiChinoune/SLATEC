!*==F0O.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK F0O
REAL FUNCTION F0O(X)
  IMPLICIT NONE
  !*--F0O5
  !***BEGIN PROLOGUE  F0O
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  F0O
  REAL X
  !***FIRST EXECUTABLE STATEMENT  F0O
  F0O = (2.0E0*SIN(X))**14
END FUNCTION F0O
