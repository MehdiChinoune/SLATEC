!*==DF2N.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DF2N
REAL(8) FUNCTION DF2N(X)
  IMPLICIT NONE
  !*--DF2N5
  !***BEGIN PROLOGUE  DF2N
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  DF2N
  REAL(8) :: X
  !***FIRST EXECUTABLE STATEMENT  DF2N
  DF2N = X**(-0.9D+00)
END FUNCTION DF2N
