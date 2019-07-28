!** DPRVEC
REAL(DP) PURE FUNCTION DPRVEC(M,U,V)
  !> Subsidiary to DBVSUP
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (PRVEC-S, DPRVEC-D)
  !***
  ! **Author:**  Watts, H. A., (SNLA)
  !***
  ! **Description:**
  !
  !  This subroutine computes the inner product of a vector U
  !  with the imaginary product or mate vector corresponding to V.
  !
  !***
  ! **See also:**  DBVSUP
  !***
  ! **Routines called:**  DDOT

  !* REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  !
  INTEGER, INTENT(IN) :: M
  REAL(DP), INTENT(IN) :: U(M), V(M)
  !
  INTEGER :: n
  !* FIRST EXECUTABLE STATEMENT  DPRVEC
  n = M/2
  DPRVEC = DOT_PRODUCT(U(n+1:2*n),V(1:n)) - DOT_PRODUCT(U(1:n),V(n+1:2*n))
  !
END FUNCTION DPRVEC