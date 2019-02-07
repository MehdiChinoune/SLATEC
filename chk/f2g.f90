!*==F2G.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK F2G
REAL FUNCTION F2G(X)
  IMPLICIT NONE
  !*--F2G5
  !***BEGIN PROLOGUE  F2G
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  F2G
  REAL X
  !***FIRST EXECUTABLE STATEMENT  F2G
  F2G = X*SIN(0.3E+02*X)*COS(0.5E+02*X)
END FUNCTION F2G
