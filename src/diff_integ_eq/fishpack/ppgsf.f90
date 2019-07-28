!** PPGSF
REAL(SP) PURE FUNCTION PPGSF(X,Iz,C,A,Bh)
  !> Subsidiary to CBLKTR
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (PPGSF-S)
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
  INTEGER :: j
  REAL(SP) :: summ
  !* FIRST EXECUTABLE STATEMENT  PPGSF
  summ = 0._SP
  DO j = 1, Iz
    summ = summ - 1._SP/(X-Bh(j))**2
  END DO
  PPGSF = summ
  !
END FUNCTION PPGSF