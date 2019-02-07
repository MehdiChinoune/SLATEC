!*==F4S.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK F4S
REAL FUNCTION F4S(X)
  IMPLICIT NONE
  !*--F4S5
  !***BEGIN PROLOGUE  F4S
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  F4S
  REAL X
  !***FIRST EXECUTABLE STATEMENT  F4S
  IF ( X==.33E+00 ) THEN
    F4S = 0.0
    GOTO 99999
  ENDIF
  F4S = ABS(X-0.33E+00)**(-0.999E+00)
  RETURN
  99999 CONTINUE
  END FUNCTION F4S
