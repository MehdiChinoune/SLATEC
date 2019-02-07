!*==F1G.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK F1G
REAL FUNCTION F1G(X)
  IMPLICIT NONE
  !*--F1G5
  !***BEGIN PROLOGUE  F1G
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  F1G
  REAL pi , X
  DATA pi/3.1415926535897932/
  !***FIRST EXECUTABLE STATEMENT  F1G
  F1G = 2.0/(2.0+SIN(10.0*pi*X))
END FUNCTION F1G
