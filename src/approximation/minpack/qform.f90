!** QFORM
SUBROUTINE QFORM(M,N,Q,Ldq,Wa)
  !> Subsidiary to SNSQ and SNSQE
  !***
  ! **Library:**   SLATEC
  !***
  ! **Type:**      SINGLE PRECISION (QFORM-S, DQFORM-D)
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
  !       SUBROUTINE QFORM(M,N,Q,LDQ,WA)
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
  !         the first min(M,N) columns of Q contains the factored form.
  !         On output Q has been accumulated into a square matrix.
  !
  !       LDQ is a positive integer input variable not less than M
  !         which specifies the leading dimension of the array Q.
  !
  !       WA is a work array of length M.
  !
  !***
  ! **See also:**  SNSQ, SNSQE
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   800301  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTIONsection.  (WRB)
  !   900328  Added TYPE section.  (WRB)

  INTEGER :: M, N, Ldq
  REAL(SP) :: Q(Ldq,M), Wa(M)
  INTEGER :: i, j, jm1, k, l, minmn, np1
  REAL(SP) :: summ, temp
  REAL(SP), PARAMETER :: one = 1._SP, zero = 0._SP
  !* FIRST EXECUTABLE STATEMENT  QFORM
  minmn = MIN(M,N)
  IF( minmn>=2 ) THEN
    DO j = 2, minmn
      jm1 = j - 1
      DO i = 1, jm1
        Q(i,j) = zero
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
        Q(i,j) = zero
      END DO
      Q(j,j) = one
    END DO
  END IF
  !
  !     ACCUMULATE Q FROM ITS FACTORED FORM.
  !
  DO l = 1, minmn
    k = minmn - l + 1
    DO i = k, M
      Wa(i) = Q(i,k)
      Q(i,k) = zero
    END DO
    Q(k,k) = one
    IF( Wa(k)/=zero ) THEN
      DO j = k, M
        summ = zero
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
  !     LAST CARD OF SUBROUTINE QFORM.
  !
END SUBROUTINE QFORM
