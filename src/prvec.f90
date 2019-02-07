!*==PRVEC.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK PRVEC
FUNCTION PRVEC(M,U,V)
  IMPLICIT NONE
  !*--PRVEC5
  !*** Start of declarations inserted by SPAG
  INTEGER M, n, np
  REAL PRVEC, SDOT, U, V, vp
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  PRVEC
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to BVSUP
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (PRVEC-S, DPRVEC-D)
  !***AUTHOR  Watts, H. A., (SNLA)
  !***DESCRIPTION
  !
  !  This subroutine computes the inner product of a vector U
  !  with the imaginary product or mate vector corresponding to V
  !
  !***SEE ALSO  BVSUP
  !***ROUTINES CALLED  SDOT
  !***REVISION HISTORY  (YYMMDD)
  !   750601  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !   910722  Updated AUTHOR section.  (ALS)
  !***END PROLOGUE  PRVEC
  !
  DIMENSION U(*), V(*)
  !***FIRST EXECUTABLE STATEMENT  PRVEC
  n = M/2
  np = n + 1
  vp = SDOT(n,U(1),1,V(np),1)
  PRVEC = SDOT(n,U(np),1,V(1),1) - vp
END FUNCTION PRVEC
