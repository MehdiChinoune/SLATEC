!*==DF3G.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DF3G
REAL(8) FUNCTION DF3G(X)
  IMPLICIT NONE
  !*--DF3G5
  !***BEGIN PROLOGUE  DF3G
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  DF3G
  REAL(8) :: X
  !***FIRST EXECUTABLE STATEMENT  DF3G
  DF3G = ABS(X-0.33D+00)**(-.90D+00)
END FUNCTION DF3G
