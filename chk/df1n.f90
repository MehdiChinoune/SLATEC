!*==DF1N.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DF1N
REAL(8) FUNCTION DF1N(X)
  IMPLICIT NONE
  !*--DF1N5
  !***BEGIN PROLOGUE  DF1N
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  DF1N
  REAL(8) :: X
  !***FIRST EXECUTABLE STATEMENT  DF1N
  DF1N = 1.0D0/(X**4+X**2+1.0D0)
END FUNCTION DF1N
