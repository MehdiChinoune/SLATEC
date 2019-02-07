!*==DF1C.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DF1C
REAL(8) FUNCTION DF1C(X)
  IMPLICIT NONE
  !*--DF1C5
  !***BEGIN PROLOGUE  DF1C
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  DF1C
  REAL(8) :: X
  !***FIRST EXECUTABLE STATEMENT  DF1C
  DF1C = 0.0D+00
  IF ( X/=0.33D+00 ) DF1C = (X-0.5D+00)*ABS(X-0.33D+00)**(-0.9D+00)
END FUNCTION DF1C
