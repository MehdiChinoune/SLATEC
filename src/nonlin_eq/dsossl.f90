!** DSOSSL
SUBROUTINE DSOSSL(K,N,L,X,C,B,M)
  !> Subsidiary to DSOS
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (SOSSOL-S, DSOSSL-D)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     DSOSSL solves an upper triangular type of linear system by back
  !     substitution.
  !
  !     The matrix C is upper trapezoidal and stored as a linear array by
  !     rows. The equations have been normalized so that the diagonal
  !     entries of C are understood to be unity. The off diagonal entries
  !     and the elements of the constant right hand side vector B have
  !     already been stored as the negatives of the corresponding equation
  !     values.
  !     With each call to DSOSSL a (K-1) by (K-1) triangular system is
  !     resolved. For L greater than K, column L of C is included in the
  !     right hand side vector.
  !
  !***
  ! **See also:**  DSOS
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   801001  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900328  Added TYPE section.  (WRB)

  !
  !
  INTEGER :: K, L, M, N
  REAL(DP) :: B(:), C(:), X(:)
  INTEGER :: j, j1, k1, k2, k3, km1, kn, lk
  REAL(DP) :: xmax
  !
  !* FIRST EXECUTABLE STATEMENT  DSOSSL
  km1 = K - 1
  lk = km1
  IF( L==K ) lk = K
  kn = M
  !
  !
  DO k1 = 1, km1
    k2 = K - k1
    k3 = k2 + 1
    xmax = 0.0D0
    kn = kn - N + k2 - 1
    IF( k3<=lk ) THEN
      j1 = kn
      !
      DO j = k3, lk
        j1 = j1 + 1
        xmax = xmax + C(j1)*X(j)
      END DO
    END IF
    !
    IF( L>K ) THEN
      j1 = kn + L - k2
      xmax = xmax + C(j1)*X(L)
    END IF
    X(k2) = xmax + B(k2)
  END DO
  !
END SUBROUTINE DSOSSL
