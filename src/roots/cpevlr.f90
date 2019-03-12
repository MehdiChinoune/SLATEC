!DECK CPEVLR
SUBROUTINE CPEVLR(N,M,A,X,C)
  IMPLICIT NONE
  REAL ci, cim1, X
  INTEGER i, j, M, mini, N, np1
  !***BEGIN PROLOGUE  CPEVLR
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to CPZERO
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (CPEVLR-S)
  !***AUTHOR  (UNKNOWN)
  !***SEE ALSO  CPZERO
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   810223  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900402  Added TYPE section.  (WRB)
  !***END PROLOGUE  CPEVLR
  REAL A(*), C(*)
  !***FIRST EXECUTABLE STATEMENT  CPEVLR
  np1 = N + 1
  DO j = 1, np1
    ci = 0.0
    cim1 = A(j)
    mini = MIN(M+1,N+2-j)
    DO i = 1, mini
      IF ( j/=1 ) ci = C(i)
      IF ( i/=1 ) cim1 = C(i-1)
      C(i) = cim1 + X*ci
    ENDDO
  ENDDO
END SUBROUTINE CPEVLR
