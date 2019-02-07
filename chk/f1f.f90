!*==F1F.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK F1F
REAL FUNCTION F1F(X)
  IMPLICIT NONE
  !*--F1F5
  !***BEGIN PROLOGUE  F1F
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  F1F
  REAL X, x1, y
  !***FIRST EXECUTABLE STATEMENT  F1F
  x1 = X + 1.0
  F1F = 5.0/x1/x1
  y = 5.0/x1
  IF ( y>3.1415926535897932 ) F1F = 0.0
END FUNCTION F1F
