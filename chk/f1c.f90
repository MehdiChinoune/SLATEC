!DECK F1C
REAL FUNCTION F1C(X)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  F1C
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  F1C
  REAL X
  !***FIRST EXECUTABLE STATEMENT  F1C
  F1C = 0.0
  IF ( X/=0.33 ) F1C = (X-0.5)*ABS(X-0.33)**(-0.9)
END FUNCTION F1C
