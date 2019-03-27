!** PPPSF
REAL FUNCTION PPPSF(X,Iz,C,A,Bh)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to CBLKTR
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
  
  REAL A(*), Bh(*), C(*), sum, X
  INTEGER Iz, j
  !* FIRST EXECUTABLE STATEMENT  PPPSF
  sum = 0.
  DO j = 1, Iz
    sum = sum + 1./(X-Bh(j))
  ENDDO
  PPPSF = sum
END FUNCTION PPPSF
