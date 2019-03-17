!DECK PPSPF
REAL FUNCTION PPSPF(X,Iz,C,A,Bh)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  PPSPF
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to BLKTRI
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (PPSPF-S)
  !***AUTHOR  (UNKNOWN)
  !***SEE ALSO  BLKTRI
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !***END PROLOGUE  PPSPF
  REAL A, Bh, C, sum, X
  INTEGER Iz, j
  DIMENSION A(*), C(*), Bh(*)
  !***FIRST EXECUTABLE STATEMENT  PPSPF
  sum = 0.
  DO j = 1, Iz
    sum = sum + 1./(X-Bh(j))
  ENDDO
  PPSPF = sum
END FUNCTION PPSPF
