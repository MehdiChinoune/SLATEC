!*==DF5S.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DF5S
REAL(8) FUNCTION DF5S(X)
  IMPLICIT NONE
  !*--DF5S5
  !***BEGIN PROLOGUE  DF5S
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  DF5S
  REAL(8) :: X
  !***FIRST EXECUTABLE STATEMENT  DF5S
  DF5S = 0.0D+00
  IF ( X/=0.0D+00 ) DF5S = 1.0D+00/(X*SQRT(X))
END FUNCTION DF5S
