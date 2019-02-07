!*==DF0S.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DF0S
REAL(8) FUNCTION DF0S(X)
  IMPLICIT NONE
  !*--DF0S5
  !***BEGIN PROLOGUE  DF0S
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  DF0S
  REAL(8) :: X
  !***FIRST EXECUTABLE STATEMENT  DF0S
  DF0S = 0.0D+00
  IF ( X/=0.0D+00 ) DF0S = 0.1D+01/SQRT(X)
END FUNCTION DF0S
