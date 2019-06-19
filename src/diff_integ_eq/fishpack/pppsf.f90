!** PPPSF
REAL(SP) FUNCTION PPPSF(X,Iz,C,A,Bh)
  !> Subsidiary to CBLKTR
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (PPPSF-S)
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

  INTEGER :: Iz
  REAL(SP) :: X, A(Iz), Bh(Iz), C(Iz)
  INTEGER :: j
  REAL(SP) :: summ
  !* FIRST EXECUTABLE STATEMENT  PPPSF
  summ = 0.
  DO j = 1, Iz
    summ = summ + 1./(X-Bh(j))
  END DO
  PPPSF = summ
END FUNCTION PPPSF
