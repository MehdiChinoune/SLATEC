!DECK PGSF
REAL FUNCTION PGSF(X,Iz,C,A,Bh)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  PGSF
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to CBLKTR
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (PGSF-S)
  !***AUTHOR  (UNKNOWN)
  !***SEE ALSO  CBLKTR
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !***END PROLOGUE  PGSF
  REAL A, Bh, C, dd, fsg, hsg, X
  INTEGER Iz, j
  DIMENSION A(*), C(*), Bh(*)
  !***FIRST EXECUTABLE STATEMENT  PGSF
  fsg = 1.
  hsg = 1.
  DO j = 1, Iz
    dd = 1./(X-Bh(j))
    fsg = fsg*A(j)*dd
    hsg = hsg*C(j)*dd
  ENDDO
  IF ( MOD(Iz,2)/=0 ) THEN
    PGSF = 1. + fsg + hsg
    RETURN
  ENDIF
  PGSF = 1. - fsg - hsg
  RETURN
END FUNCTION PGSF
