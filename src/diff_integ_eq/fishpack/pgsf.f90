!** PGSF
REAL(SP) PURE FUNCTION PGSF(X,Iz,C,A,Bh)
  !> Subsidiary to CBLKTR
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

  INTEGER, INTENT(IN) :: Iz
  REAL(SP), INTENT(IN) :: X, A(Iz), Bh(Iz), C(Iz)
  !
  INTEGER :: j
  REAL(SP) :: dd, fsg, hsg
  !* FIRST EXECUTABLE STATEMENT  PGSF
  fsg = 1._SP
  hsg = 1._SP
  DO j = 1, Iz
    dd = 1._SP/(X-Bh(j))
    fsg = fsg*A(j)*dd
    hsg = hsg*C(j)*dd
  END DO
  IF( MOD(Iz,2)/=0 ) THEN
    PGSF = 1._SP + fsg + hsg
    RETURN
  END IF
  PGSF = 1._SP - fsg - hsg
  !
  RETURN
END FUNCTION PGSF