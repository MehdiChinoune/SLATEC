!** GVEC
SUBROUTINE GVEC(X,G)
  !>
  !***
  !  Subsidiary to
  !***
  ! **Library:**   SLATEC
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   ??????  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  
  REAL G(*), X
  !* FIRST EXECUTABLE STATEMENT  GVEC
  G(1) = 0.0
  G(2) = 1.0 + COS(X)
END SUBROUTINE GVEC
