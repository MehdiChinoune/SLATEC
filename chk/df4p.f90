!*==DF4P.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DF4P
REAL(8) FUNCTION DF4P(X)
  IMPLICIT NONE
  !*--DF4P5
  !***BEGIN PROLOGUE  DF4P
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  DF4P
  REAL(8) :: X
  !***FIRST EXECUTABLE STATEMENT  DF4P
  DF4P = 0.0D+00
  IF ( X>0.0D+00 ) DF4P = 0.1D+01/(X*SQRT(X))
END FUNCTION DF4P
