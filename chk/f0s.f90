!*==F0S.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK F0S
REAL FUNCTION F0S(X)
  IMPLICIT NONE
  !*--F0S5
  !***BEGIN PROLOGUE  F0S
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  F0S
  REAL X
  !***FIRST EXECUTABLE STATEMENT  F0S
  F0S = 0.0
  IF ( X/=0.0 ) F0S = 1.0/SQRT(X)
END FUNCTION F0S
