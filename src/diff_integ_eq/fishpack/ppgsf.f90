!** PPGSF
REAL FUNCTION PPGSF(X,Iz,C,A,Bh)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to CBLKTR
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
  
  REAL A(*), Bh(*), C(*), summ, X
  INTEGER Iz, j
  !* FIRST EXECUTABLE STATEMENT  PPGSF
  summ = 0.
  DO j = 1, Iz
    summ = summ - 1./(X-Bh(j))**2
  END DO
  PPGSF = summ
END FUNCTION PPGSF
