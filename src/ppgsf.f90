!DECK PPGSF
FUNCTION PPGSF(X,Iz,C,A,Bh)
  IMPLICIT NONE
  REAL A, Bh, C, PPGSF, sum, X
  INTEGER Iz, j
  !***BEGIN PROLOGUE  PPGSF
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to CBLKTR
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (PPGSF-S)
  !***AUTHOR  (UNKNOWN)
  !***SEE ALSO  CBLKTR
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !***END PROLOGUE  PPGSF
  DIMENSION A(*), C(*), Bh(*)
  !***FIRST EXECUTABLE STATEMENT  PPGSF
  sum = 0.
  DO j = 1, Iz
    sum = sum - 1./(X-Bh(j))**2
  ENDDO
  PPGSF = sum
END FUNCTION PPGSF
