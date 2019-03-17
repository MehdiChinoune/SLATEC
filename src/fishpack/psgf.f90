!DECK PSGF
REAL FUNCTION PSGF(X,Iz,C,A,Bh)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  PSGF
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to BLKTRI
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (PSGF-S)
  !***AUTHOR  (UNKNOWN)
  !***SEE ALSO  BLKTRI
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !***END PROLOGUE  PSGF
  REAL A, Bh, C, dd, fsg, hsg, X
  INTEGER Iz, j
  DIMENSION A(*), C(*), Bh(*)
  !***FIRST EXECUTABLE STATEMENT  PSGF
  fsg = 1.
  hsg = 1.
  DO j = 1, Iz
    dd = 1./(X-Bh(j))
    fsg = fsg*A(j)*dd
    hsg = hsg*C(j)*dd
  ENDDO
  IF ( MOD(Iz,2)/=0 ) THEN
    PSGF = 1. + fsg + hsg
    RETURN
  ENDIF
  PSGF = 1. - fsg - hsg
  RETURN
END FUNCTION PSGF