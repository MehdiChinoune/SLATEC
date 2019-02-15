!DECK QRFAC
SUBROUTINE QRFAC(M,N,A,Lda,Pivot,Ipvt,Lipvt,Sigma,Acnorm,Wa)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  QRFAC
  !***SUBSIDIARY
  !***PURPOSE  Subsidiary to SNLS1, SNLS1E, SNSQ and SNSQE
  !***LIBRARY   SLATEC
  !***TYPE      SINGLE PRECISION (QRFAC-S, DQRFAC-D)
  !***AUTHOR  (UNKNOWN)
  !***DESCRIPTION
  !
  !     This subroutine uses Householder transformations with column
  !     pivoting (optional) to compute a QR factorization of the
  !     M by N matrix A. That is, QRFAC determines an orthogonal
  !     matrix Q, a permutation matrix P, and an upper trapezoidal
  !     matrix R with diagonal elements of nonincreasing magnitude,
  !     such that A*P = Q*R. The Householder transformation for
  !     column K, K = 1,2,...,MIN(M,N), is of the form
  !
  !                           T
  !           I - (1/U(K))*U*U
  !
  !     where U has zeros in the first K-1 positions. The form of
  !     this transformation and the method of pivoting first
  !     appeared in the corresponding LINPACK subroutine.
  !
  !     The subroutine statement is
  !
  !       SUBROUTINE QRFAC(M,N,A,LDA,PIVOT,IPVT,LIPVT,SIGMA,ACNORM,WA)
  !
  !     where
  !
  !       M is a positive integer input variable set to the number
  !         of rows of A.
  !
  !       N is a positive integer input variable set to the number
  !         of columns of A.
  !
  !       A is an M by N array. On input A contains the matrix for
  !         which the QR factorization is to be computed. On output
  !         the strict upper trapezoidal part of A contains the strict
  !         upper trapezoidal part of R, and the lower trapezoidal
  !         part of A contains a factored form of Q (the non-trivial
  !         elements of the U vectors described above).
  !
  !       LDA is a positive integer input variable not less than M
  !         which specifies the leading dimension of the array A.
  !
  !       PIVOT is a logical input variable. If pivot is set .TRUE.,
  !         then column pivoting is enforced. If pivot is set .FALSE.,
  !         then no column pivoting is done.
  !
  !       IPVT is an integer output array of length LIPVT. IPVT
  !         defines the permutation matrix P such that A*P = Q*R.
  !         Column J of P is column IPVT(J) of the identity matrix.
  !         If pivot is .FALSE., IPVT is not referenced.
  !
  !       LIPVT is a positive integer input variable. If PIVOT is
  !             .FALSE., then LIPVT may be as small as 1. If PIVOT is
  !             .TRUE., then LIPVT must be at least N.
  !
  !       SIGMA is an output array of length N which contains the
  !         diagonal elements of R.
  !
  !       ACNORM is an output array of length N which contains the
  !         norms of the corresponding columns of the input matrix A.
  !         If this information is not needed, then ACNORM can coincide
  !         with SIGMA.
  !
  !       WA is a work array of length N. If pivot is .FALSE., then WA
  !         can coincide with SIGMA.
  !
  !***SEE ALSO  SNLS1, SNLS1E, SNSQ, SNSQE
  !***ROUTINES CALLED  ENORM, R1MACH
  !***REVISION HISTORY  (YYMMDD)
  !   800301  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   900328  Added TYPE section.  (WRB)
  !***END PROLOGUE  QRFAC
  INTEGER M, N, Lda, Lipvt
  INTEGER Ipvt(*)
  LOGICAL Pivot
  REAL A(Lda,*), Sigma(*), Acnorm(*), Wa(*)
  INTEGER i, j, jp1, k, kmax, minmn
  REAL ajnorm, epsmch, one, p05, sum, temp, zero
  REAL R1MACH, ENORM
  SAVE one, p05, zero
  DATA one, p05, zero/1.0E0, 5.0E-2, 0.0E0/
  !***FIRST EXECUTABLE STATEMENT  QRFAC
  epsmch = R1MACH(4)
  !
  !     COMPUTE THE INITIAL COLUMN NORMS AND INITIALIZE SEVERAL ARRAYS.
  !
  DO j = 1, N
    Acnorm(j) = ENORM(M,A(1,j))
    Sigma(j) = Acnorm(j)
    Wa(j) = Sigma(j)
    IF ( Pivot ) Ipvt(j) = j
  ENDDO
  !
  !     REDUCE A TO R WITH HOUSEHOLDER TRANSFORMATIONS.
  !
  minmn = MIN(M,N)
  DO j = 1, minmn
    IF ( Pivot ) THEN
      !
      !        BRING THE COLUMN OF LARGEST NORM INTO THE PIVOT POSITION.
      !
      kmax = j
      DO k = j, N
        IF ( Sigma(k)>Sigma(kmax) ) kmax = k
      ENDDO
      IF ( kmax/=j ) THEN
        DO i = 1, M
          temp = A(i,j)
          A(i,j) = A(i,kmax)
          A(i,kmax) = temp
        ENDDO
        Sigma(kmax) = Sigma(j)
        Wa(kmax) = Wa(j)
        k = Ipvt(j)
        Ipvt(j) = Ipvt(kmax)
        Ipvt(kmax) = k
      ENDIF
    ENDIF
    !
    !        COMPUTE THE HOUSEHOLDER TRANSFORMATION TO REDUCE THE
    !        J-TH COLUMN OF A TO A MULTIPLE OF THE J-TH UNIT VECTOR.
    !
    ajnorm = ENORM(M-j+1,A(j,j))
    IF ( ajnorm/=zero ) THEN
      IF ( A(j,j)<zero ) ajnorm = -ajnorm
      DO i = j, M
        A(i,j) = A(i,j)/ajnorm
      ENDDO
      A(j,j) = A(j,j) + one
      !
      !        APPLY THE TRANSFORMATION TO THE REMAINING COLUMNS
      !        AND UPDATE THE NORMS.
      !
      jp1 = j + 1
      IF ( N>=jp1 ) THEN
        DO k = jp1, N
          sum = zero
          DO i = j, M
            sum = sum + A(i,j)*A(i,k)
          ENDDO
          temp = sum/A(j,j)
          DO i = j, M
            A(i,k) = A(i,k) - temp*A(i,j)
          ENDDO
          IF ( .NOT.(.NOT.Pivot.OR.Sigma(k)==zero) ) THEN
            temp = A(j,k)/Sigma(k)
            Sigma(k) = Sigma(k)*SQRT(MAX(zero,one-temp**2))
            IF ( p05*(Sigma(k)/Wa(k))**2<=epsmch ) THEN
              Sigma(k) = ENORM(M-j,A(jp1,k))
              Wa(k) = Sigma(k)
            ENDIF
          ENDIF
        ENDDO
      ENDIF
    ENDIF
    Sigma(j) = -ajnorm
  ENDDO
  !
  !     LAST CARD OF SUBROUTINE QRFAC.
  !
END SUBROUTINE QRFAC
