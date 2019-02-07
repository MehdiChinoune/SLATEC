!*==DGVEC.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DGVEC
SUBROUTINE DGVEC(X,G)
  IMPLICIT NONE
  !*--DGVEC5
  !***BEGIN PROLOGUE  DGVEC
  !***PURPOSE  Subsidiary to
  !***LIBRARY   SLATEC
  !***AUTHOR  (UNKNOWN)
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  DGVEC
  REAL(8) :: X , G
  DIMENSION G(*)
  !***FIRST EXECUTABLE STATEMENT  DGVEC
  G(1) = 0.0D0
  G(2) = 1.0D0 + COS(X)
END SUBROUTINE DGVEC
