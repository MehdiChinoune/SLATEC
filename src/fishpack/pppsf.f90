!DECK PPPSF
REAL FUNCTION PPPSF(X,Iz,C,A,Bh)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  PPPSF
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to CBLKTR
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (PPPSF-S)
  !***AUTHOR  (UNKNOWN)
  !***SEE ALSO  CBLKTR
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !***END PROLOGUE  PPPSF
  REAL A, Bh, C, sum, X
  INTEGER Iz, j
  DIMENSION A(*), C(*), Bh(*)
  !***FIRST EXECUTABLE STATEMENT  PPPSF
  sum = 0.
  DO j = 1, Iz
    sum = sum + 1./(X-Bh(j))
  ENDDO
  PPPSF = sum
END FUNCTION PPPSF
