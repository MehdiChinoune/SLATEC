!*==DF2G.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DF2G
REAL(8) FUNCTION DF2G(X)
  IMPLICIT NONE
  !*--DF2G5
  !***BEGIN PROLOGUE  DF2G
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  DF2G
  REAL(8) :: X
  !***FIRST EXECUTABLE STATEMENT  DF2G
  DF2G = X*SIN(0.3D+02*X)*COS(0.5D+02*X)
END FUNCTION DF2G
