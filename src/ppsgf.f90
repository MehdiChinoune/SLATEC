!*==PPSGF.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK PPSGF
FUNCTION PPSGF(X,Iz,C,A,Bh)
  IMPLICIT NONE
  !*--PPSGF5
  !*** Start of declarations inserted by SPAG
  REAL A, Bh, C, PPSGF, sum, X
  INTEGER Iz, j
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  PPSGF
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to BLKTRI
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (PPSGF-S)
  !***AUTHOR  (UNKNOWN)
  !***SEE ALSO  BLKTRI
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !***END PROLOGUE  PPSGF
  DIMENSION A(*), C(*), Bh(*)
  !***FIRST EXECUTABLE STATEMENT  PPSGF
  sum = 0.
  DO j = 1, Iz
    sum = sum - 1./(X-Bh(j))**2
  ENDDO
  PPSGF = sum
END FUNCTION PPSGF
