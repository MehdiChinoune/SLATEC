!*==DF4S.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DF4S
REAL(8) FUNCTION DF4S(X)
  IMPLICIT NONE
  !*--DF4S5
  !***BEGIN PROLOGUE  DF4S
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  DF4S
  REAL(8) :: X
  !***FIRST EXECUTABLE STATEMENT  DF4S
  DF4S = 0.00D+00
  IF ( X-0.33D+00/=0.00D+00 ) DF4S = ABS(X-0.33D+00)**(-0.999D+00)
END FUNCTION DF4S
