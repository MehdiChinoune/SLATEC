!*==DF2S.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DF2S
REAL(8) FUNCTION DF2S(X)
  IMPLICIT NONE
  !*--DF2S5
  !***BEGIN PROLOGUE  DF2S
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  DF2S
  REAL(8) :: X
  !***FIRST EXECUTABLE STATEMENT  DF2S
  DF2S = 0.1D+03
  IF ( X/=0.0D+00 ) DF2S = SIN(0.314159D+03*X)/(0.314159D+01*X)
END FUNCTION DF2S
