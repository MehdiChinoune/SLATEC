!** PPSPF
REAL FUNCTION PPSPF(X,Iz,C,A,Bh)
  !>
  !  Subsidiary to BLKTRI
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (PPSPF-S)
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
  REAL :: summ
  !* FIRST EXECUTABLE STATEMENT  PPSPF
  summ = 0.
  DO j = 1, Iz
    summ = summ + 1./(X-Bh(j))
  END DO
  PPSPF = summ
END FUNCTION PPSPF
