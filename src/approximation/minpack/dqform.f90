!** DQFORM
PURE SUBROUTINE DQFORM(M,N,Q,Ldq,Wa)
  !> Subsidiary to DNSQ and DNSQE
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      DOUBLE PRECISION (QFORM-S, DQFORM-D)
  !***
  ! **Author:**  (UNKNOWN)
  !***
  ! **Description:**
  !
  !     This subroutine proceeds from the computed QR factorization of
  !     an M by N matrix A to accumulate the M by M orthogonal matrix
  !     Q from its factored form.
  !
  !     The subroutine statement is
  !
  !       SUBROUTINE DQFORM(M,N,Q,LDQ,WA)
  !
  !     where
  !
  !       M is a positive integer input variable set to the number
  !         of rows of A and the order of Q.
  !
  !       N is a positive integer input variable set to the number
  !         of columns of A.
  !
  !       Q is an M by M array. On input the full lower trapezoid in
  !         the first MIN(M,N) columns of Q contains the factored form.
  !         On output Q has been accumulated into a square matrix.
  !
  !       LDQ is a positive integer input variable not less than M
  !         which specifies the leading dimension of the array Q.
  !
  !       WA is a work array of length M.
  !
  !***
  ! **See also:**  DNSQ, DNSQE
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   800301  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTIONsection.  (WRB)
  !   900328  Added TYPE section.  (WRB)

  INTEGER, INTENT(IN) :: Ldq, M, N
  REAL(DP), INTENT(INOUT) :: Q(Ldq,M)
  REAL(DP), INTENT(OUT) :: Wa(M)
  !
  INTEGER :: i, j, jm1, k, l, minmn, np1
  REAL(DP) :: summ, temp
  !
  !     ZERO OUT UPPER TRIANGLE OF Q IN THE FIRST MIN(M,N) COLUMNS.
  !
  !* FIRST EXECUTABLE STATEMENT  DQFORM
  minmn = MIN(M,N)
  IF( minmn>=2 ) THEN
    DO j = 2, minmn
      jm1 = j - 1
      DO i = 1, jm1
        Q(i,j) = 0._DP
      END DO
    END DO
  END IF
  !
  !     INITIALIZE REMAINING COLUMNS TO THOSE OF THE IDENTITY MATRIX.
  !
  np1 = N + 1
  IF( M>=np1 ) THEN
    DO j = np1, M
      DO i = 1, M
        Q(i,j) = 0._DP
      END DO
      Q(j,j) = 1._DP
    END DO
  END IF
  !
  !     ACCUMULATE Q FROM ITS FACTORED FORM.
  !
  DO l = 1, minmn
    k = minmn - l + 1
    DO i = k, M
      Wa(i) = Q(i,k)
      Q(i,k) = 0._DP
    END DO
    Q(k,k) = 1._DP
    IF( Wa(k)/=0._DP ) THEN
      DO j = k, M
        summ = 0._DP
        DO i = k, M
          summ = summ + Q(i,j)*Wa(i)
        END DO
        temp = summ/Wa(k)
        DO i = k, M
          Q(i,j) = Q(i,j) - temp*Wa(i)
        END DO
      END DO
    END IF
  END DO
  !
  !     LAST CARD OF SUBROUTINE DQFORM.
  !
END SUBROUTINE DQFORM