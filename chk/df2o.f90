!*==DF2O.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DF2O
REAL(8) FUNCTION DF2O(X)
  IMPLICIT NONE
  !*--DF2O5
  !***BEGIN PROLOGUE  DF2O
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  DF2O
  REAL(8) :: X
  !***FIRST EXECUTABLE STATEMENT  DF2O
  DF2O = 0.0D+00
  IF ( X/=0.0D+00 ) DF2O = 0.1D+01/(X*X*SQRT(X))
END FUNCTION DF2O
