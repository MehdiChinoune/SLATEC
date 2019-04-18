!** PRVEC
REAL FUNCTION PRVEC(M,U,V)
  !>
  !***
  !  Subsidiary to BVSUP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (PRVEC-S, DPRVEC-D)
  !***
  ! **Author:**  Watts, H. A., (SNLA)
  !***
  ! **Description:**
  !
  !  This subroutine computes the inner product of a vector U
  !  with the imaginary product or mate vector corresponding to V
  !
  !***
  ! **See also:**  BVSUP
  !***
  ! **Routines called:**  SDOT

  !* REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  USE linear, ONLY : SDOT
  INTEGER M, n, np
  REAL U(*), V(*), vp
  !* FIRST EXECUTABLE STATEMENT  PRVEC
  n = M/2
  np = n + 1
  vp = SDOT(n,U(1),1,V(np),1)
  PRVEC = SDOT(n,U(np),1,V(1),1) - vp
END FUNCTION PRVEC
