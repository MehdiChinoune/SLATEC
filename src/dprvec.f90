!*==DPRVEC.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DPRVEC
REAL(8) FUNCTION DPRVEC(M,U,V)
  IMPLICIT NONE
  !*--DPRVEC5
  !***BEGIN PROLOGUE  DPRVEC
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to DBVSUP
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (PRVEC-S, DPRVEC-D)
  !***AUTHOR  Watts, H. A., (SNLA)
  !***DESCRIPTION
  !
  !  This subroutine computes the inner product of a vector U
  !  with the imaginary product or mate vector corresponding to V.
  !
  !***SEE ALSO  DBVSUP
  !***ROUTINES CALLED  DDOT
  !***REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  !***END PROLOGUE  DPRVEC
  !
  REAL(8) :: DDOT
  INTEGER M, n, np
  REAL(8) :: U(*), V(*), vp
  !***FIRST EXECUTABLE STATEMENT  DPRVEC
  n = M/2
  np = n + 1
  vp = DDOT(n,U(1),1,V(np),1)
  DPRVEC = DDOT(n,U(np),1,V(1),1) - vp
END FUNCTION DPRVEC
