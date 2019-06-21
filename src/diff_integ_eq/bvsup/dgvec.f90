!** DGVEC
SUBROUTINE DGVEC(X,G)
  !> Subsidiary to
  !***
  ! **Library:**   SLATEC
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)

  REAL(DP) :: X, G(2)
  !* FIRST EXECUTABLE STATEMENT  DGVEC
  G(1) = 0._DP
  G(2) = 1._DP + COS(X)
END SUBROUTINE DGVEC
