!** DPRVEC
REAL(8) FUNCTION DPRVEC(M,U,V)
  IMPLICIT NONE
  !>
  !***
  !  Subsidiary to DBVSUP
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
  REAL(8) :: DDOT
  INTEGER M, n, np
  REAL(8) :: U(*), V(*), vp
  !* FIRST EXECUTABLE STATEMENT  DPRVEC
  n = M/2
  np = n + 1
  vp = DDOT(n,U(1),1,V(np),1)
  DPRVEC = DDOT(n,U(np),1,V(1),1) - vp
END FUNCTION DPRVEC
