!*==QFORM.f90  processed by SPAG 6.72Dc at 11:01 on  6 Feb 2019
!DECK QFORM
SUBROUTINE QFORM(M,N,Q,Ldq,Wa)
  IMPLICIT NONE
  !*--QFORM5
  !***BEGIN PROLOGUE  QFORM
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to SNSQ and SNSQE
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (QFORM-S, DQFORM-D)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
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
  !***SEE ALSO  SNSQ, SNSQE
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   800301  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   900328  Added TYPE section.  (WRB)
  !***END PROLOGUE  QFORM
  INTEGER M, N, Ldq
  REAL Q(Ldq,*), Wa(*)
  INTEGER i, j, jm1, k, l, minmn, np1
  REAL one, sum, temp, zero
  SAVE one, zero
  DATA one, zero/1.0E0, 0.0E0/
  !***FIRST EXECUTABLE STATEMENT  QFORM
  minmn = MIN(M,N)
  IF ( minmn>=2 ) THEN
    DO j = 2, minmn
      jm1 = j - 1
      DO i = 1, jm1
        Q(i,j) = zero
      ENDDO
    ENDDO
  ENDIF
  !
  !     INITIALIZE REMAINING COLUMNS TO THOSE OF THE IDENTITY MATRIX.
  !
  np1 = N + 1
  IF ( M>=np1 ) THEN
    DO j = np1, M
      DO i = 1, M
        Q(i,j) = zero
      ENDDO
      Q(j,j) = one
    ENDDO
  ENDIF
  !
  !     ACCUMULATE Q FROM ITS FACTORED FORM.
  !
  DO l = 1, minmn
    k = minmn - l + 1
    DO i = k, M
      Wa(i) = Q(i,k)
      Q(i,k) = zero
    ENDDO
    Q(k,k) = one
    IF ( Wa(k)/=zero ) THEN
      DO j = k, M
        sum = zero
        DO i = k, M
          sum = sum + Q(i,j)*Wa(i)
        ENDDO
        temp = sum/Wa(k)
        DO i = k, M
          Q(i,j) = Q(i,j) - temp*Wa(i)
        ENDDO
      ENDDO
    ENDIF
  ENDDO
  !
  !     LAST CARD OF SUBROUTINE QFORM.
  !
END SUBROUTINE QFORM
