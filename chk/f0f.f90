!*==F0F.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK F0F
REAL FUNCTION F0F(X)
  IMPLICIT NONE
  !*--F0F5
  !***BEGIN PROLOGUE  F0F
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  F0F
  REAL X
  !***FIRST EXECUTABLE STATEMENT  F0F
  F0F = 0.0
  IF ( X/=0.0 ) F0F = SIN(0.5E+02*X)/(X*SQRT(X))
END FUNCTION F0F
