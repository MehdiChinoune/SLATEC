!** PSGF
REAL FUNCTION PSGF(X,Iz,C,A,Bh)
  !>
  !  Subsidiary to BLKTRI
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (PSGF-S)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **See also:**  BLKTRI
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)

  INTEGER :: Iz
  REAL :: X, A(Iz), Bh(Iz), C(Iz)
  INTEGER :: j
  REAL :: dd, fsg, hsg
  !* FIRST EXECUTABLE STATEMENT  PSGF
  fsg = 1.
  hsg = 1.
  DO j = 1, Iz
    dd = 1./(X-Bh(j))
    fsg = fsg*A(j)*dd
    hsg = hsg*C(j)*dd
  END DO
  IF ( MOD(Iz,2)/=0 ) THEN
    PSGF = 1. + fsg + hsg
    RETURN
  END IF
  PSGF = 1. - fsg - hsg
  RETURN
END FUNCTION PSGF
