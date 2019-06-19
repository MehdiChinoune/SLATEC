!** DBDIFF
SUBROUTINE DBDIFF(L,V)
  !> Subsidiary to DBSKIN
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (BDIFF-S, DBDIFF-D)
  !***
  ! **Author:**  Amos, D. E., (SNLA)
  !***
  ! **Description:**
  !
  !     DBDIFF computes the sum of B(L,K)*V(K)*(-1)**K where B(L,K)
  !     are the binomial coefficients.  Truncated sums are computed by
  !     setting last part of the V vector to zero. On return, the binomial
  !     sum is in V(L).
  !
  !***
  ! **See also:**  DBSKIN
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   820601  DATE WRITTEN
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)

  !
  INTEGER :: L
  REAL(DP) :: V(L)
  INTEGER :: i, j, k
  !* FIRST EXECUTABLE STATEMENT  DBDIFF
  IF( L==1 ) RETURN
  DO j = 2, L
    k = L
    DO i = j, L
      V(k) = V(k-1) - V(k)
      k = k - 1
    END DO
  END DO
END SUBROUTINE DBDIFF
