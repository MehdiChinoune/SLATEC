!*==DBDIFF.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK DBDIFF
SUBROUTINE DBDIFF(L,V)
  IMPLICIT NONE
  !*--DBDIFF5
  !***BEGIN PROLOGUE  DBDIFF
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to DBSKIN
  !***LIBRARY   SLATEC
  !***TYPE      DOUBLE PRECISION (BDIFF-S, DBDIFF-D)
  !***AUTHOR  Amos, D. E., (SNLA)
  !***DESCRIPTION
  !
  !     DBDIFF computes the sum of B(L,K)*V(K)*(-1)**K where B(L,K)
  !     are the binomial coefficients.  Truncated sums are computed by
  !     setting last part of the V vector to zero. On return, the binomial
  !     sum is in V(L).
  !
  !***SEE ALSO  DBSKIN
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   820601  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)
  !***END PROLOGUE  DBDIFF
  !
  INTEGER i , j , k , L
  REAL(8) :: V
  DIMENSION V(*)
  !***FIRST EXECUTABLE STATEMENT  DBDIFF
  IF ( L==1 ) RETURN
  DO j = 2 , L
    k = L
    DO i = j , L
      V(k) = V(k-1) - V(k)
      k = k - 1
    ENDDO
  ENDDO
END SUBROUTINE DBDIFF
