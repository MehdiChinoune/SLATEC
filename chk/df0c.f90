!*==DF0C.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DF0C
REAL(8) FUNCTION DF0C(X)
  IMPLICIT NONE
  !*--DF0C5
  !***BEGIN PROLOGUE  DF0C
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  DF0C
  REAL(8) :: X
  !***FIRST EXECUTABLE STATEMENT  DF0C
  DF0C = 1.D0/(X*X+1.D-4)
END FUNCTION DF0C
