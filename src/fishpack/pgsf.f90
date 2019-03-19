!** PGSF
REAL FUNCTION PGSF(X,Iz,C,A,Bh)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to CBLKTR
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (PGSF-S)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **See also:**  CBLKTR
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  
  REAL A, Bh, C, dd, fsg, hsg, X
  INTEGER Iz, j
  DIMENSION A(*), C(*), Bh(*)
  !* FIRST EXECUTABLE STATEMENT  PGSF
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
