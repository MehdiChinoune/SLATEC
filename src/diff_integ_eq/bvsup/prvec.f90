!** PRVEC
REAL FUNCTION PRVEC(M,U,V)
  !>
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
  INTEGER :: M
  REAL :: U(M), V(M)
  INTEGER :: n
  !* FIRST EXECUTABLE STATEMENT  PRVEC
  n = M/2
  PRVEC = DOT_PRODUCT(U(n+1:2*n),V(1:n)) - DOT_PRODUCT(U(1:n),V(n+1:2*n))
END FUNCTION PRVEC
