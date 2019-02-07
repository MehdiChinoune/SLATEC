!*==F1N.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK F1N
REAL FUNCTION F1N(X)
  IMPLICIT NONE
  !*--F1N5
  !***BEGIN PROLOGUE  F1N
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  F1N
  REAL X
  !***FIRST EXECUTABLE STATEMENT  F1N
  F1N = 1.0E0/(X**4+X**2+1.0E0)
END FUNCTION F1N
