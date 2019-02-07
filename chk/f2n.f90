!*==F2N.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK F2N
REAL FUNCTION F2N(X)
  IMPLICIT NONE
  !*--F2N5
  !***BEGIN PROLOGUE  F2N
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  F2N
  REAL X
  !***FIRST EXECUTABLE STATEMENT  F2N
  F2N = X**(-0.9E+00)
END FUNCTION F2N
