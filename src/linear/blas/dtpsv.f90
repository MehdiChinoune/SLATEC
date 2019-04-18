!** DTPSV
SUBROUTINE DTPSV(Uplo,Trans,Diag,N,Ap,X,Incx)
  !>
  !  Solve one of the systems of equations.
  !***
  ! **Library:**   SLATEC (BLAS)
  !***
  ! **Category:**  D1B4
  !***
  ! **Type:**      DOUBLE PRECISION (STPSV-S, DTPSV-D, CTPSV-C)
  !***
  ! **Keywords:**  LEVEL 2 BLAS, LINEAR ALGEBRA
  !***
  ! **Author:**  Dongarra, J. J., (ANL)
  !           Du Croz, J., (NAG)
  !           Hammarling, S., (NAG)
  !           Hanson, R. J., (SNLA)
  !***
  ! **Description:**
  !
  !  DTPSV  solves one of the systems of equations
  !
  !     A*x = b,   or   A'*x = b,
  !
  !  where b and x are n element vectors and A is an n by n unit, or
  !  non-unit, upper or lower triangular matrix, supplied in packed form.
  !
  !  No test for singularity or near-singularity is included in this
  !  routine. Such tests must be performed before calling this routine.
  !
  !  Parameters
  !  ==========
  !
  !  UPLO   - CHARACTER*1.
  !           On entry, UPLO specifies whether the matrix is an upper or
  !           lower triangular matrix as follows:
  !
  !              UPLO = 'U' or 'u'   A is an upper triangular matrix.
  !
  !              UPLO = 'L' or 'l'   A is a lower triangular matrix.
  !
  !           Unchanged on exit.
  !
  !  TRANS  - CHARACTER*1.
  !           On entry, TRANS specifies the equations to be solved as
  !           follows:
  !
  !              TRANS = 'N' or 'n'   A*x = b.
  !
  !              TRANS = 'T' or 't'   A'*x = b.
  !
  !              TRANS = 'C' or 'c'   A'*x = b.
  !
  !           Unchanged on exit.
  !
  !  DIAG   - CHARACTER*1.
  !           On entry, DIAG specifies whether or not A is unit
  !           triangular as follows:
  !
  !              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
  !
  !              DIAG = 'N' or 'n'   A is not assumed to be unit
  !                                  triangular.
  !
  !           Unchanged on exit.
  !
  !  N      - INTEGER.
  !           On entry, N specifies the order of the matrix A.
  !           N must be at least zero.
  !           Unchanged on exit.
  !
  !  AP     - DOUBLE PRECISION array of DIMENSION at least
  !           ( ( n*( n + 1))/2).
  !           Before entry with  UPLO = 'U' or 'u', the array AP must
  !           contain the upper triangular matrix packed sequentially,
  !           column by column, so that AP( 1 ) contains a( 1, 1 ),
  !           AP( 2 ) and AP( 3 ) contain a( 1, 2 ) and a( 2, 2 )
  !           respectively, and so on.
  !           Before entry with UPLO = 'L' or 'l', the array AP must
  !           contain the lower triangular matrix packed sequentially,
  !           column by column, so that AP( 1 ) contains a( 1, 1 ),
  !           AP( 2 ) and AP( 3 ) contain a( 2, 1 ) and a( 3, 1 )
  !           respectively, and so on.
  !           Note that when  DIAG = 'U' or 'u', the diagonal elements of
  !           A are not referenced, but are assumed to be unity.
  !           Unchanged on exit.
  !
  !  X      - DOUBLE PRECISION array of dimension at least
  !           ( 1 + ( n - 1 )*abs( INCX ) ).
  !           Before entry, the incremented array X must contain the n
  !           element right-hand side vector b. On exit, X is overwritten
  !           with the solution vector x.
  !
  !  INCX   - INTEGER.
  !           On entry, INCX specifies the increment for the elements of
  !           X. INCX must not be zero.
  !           Unchanged on exit.
  !
  !***
  ! **References:**  Dongarra, J. J., Du Croz, J., Hammarling, S., and
  !                 Hanson, R. J.  An extended set of Fortran basic linear
  !                 algebra subprograms.  ACM TOMS, Vol. 14, No. 1,
  !                 pp. 1-17, March 1988.
  !***
  ! **Routines called:**  LSAME, XERBLA

  !* REVISION HISTORY  (YYMMDD)
  !   861022  DATE WRITTEN
  !   910605  Modified to meet SLATEC prologue standards.  Only comment
  !           lines were modified.  (BKS)
  USE service, ONLY : XERBLA
  !     .. Scalar Arguments ..
  INTEGER Incx, N
  CHARACTER :: Diag, Trans, Uplo
  !     .. Array Arguments ..
  REAL(8) :: Ap(*), X(*)
  !     .. Parameters ..
  REAL(8), PARAMETER :: ZERO = 0.0D+0
  !     .. Local Scalars ..
  REAL(8) :: temp
  INTEGER i, info, ix, j, jx, k, kk, kx
  LOGICAL nounit
  !* FIRST EXECUTABLE STATEMENT  DTPSV
  !
  !     Test the input parameters.
  !
  info = 0
  IF ( .NOT.LSAME(Uplo,'U').AND..NOT.LSAME(Uplo,'L') ) THEN
    info = 1
  ELSEIF ( .NOT.LSAME(Trans,'N').AND..NOT.LSAME(Trans,'T').AND.&
      .NOT.LSAME(Trans,'C') ) THEN
    info = 2
  ELSEIF ( .NOT.LSAME(Diag,'U').AND..NOT.LSAME(Diag,'N') ) THEN
    info = 3
  ELSEIF ( N<0 ) THEN
    info = 4
  ELSEIF ( Incx==0 ) THEN
    info = 7
  END IF
  IF ( info/=0 ) THEN
    CALL XERBLA('DTPSV ',info)
    RETURN
  END IF
  !
  !     Quick return if possible.
  !
  IF ( N==0 ) RETURN
  !
  nounit = LSAME(Diag,'N')
  !
  !     Set up the start point in X if the increment is not unity. This
  !     will be  ( N - 1 )*INCX  too small for descending loops.
  !
  IF ( Incx<=0 ) THEN
    kx = 1 - (N-1)*Incx
  ELSEIF ( Incx/=1 ) THEN
    kx = 1
  END IF
  !
  !     Start the operations. In this version the elements of AP are
  !     accessed sequentially with one pass through AP.
  !
  IF ( LSAME(Trans,'N') ) THEN
    !
    !        Form  x := inv( A )*x.
    !
    IF ( LSAME(Uplo,'U') ) THEN
      kk = (N*(N+1))/2
      IF ( Incx==1 ) THEN
        DO j = N, 1, -1
          IF ( X(j)/=ZERO ) THEN
            IF ( nounit ) X(j) = X(j)/Ap(kk)
            temp = X(j)
            k = kk - 1
            DO i = j - 1, 1, -1
              X(i) = X(i) - temp*Ap(k)
              k = k - 1
            END DO
          END IF
          kk = kk - j
        END DO
      ELSE
        jx = kx + (N-1)*Incx
        DO j = N, 1, -1
          IF ( X(jx)/=ZERO ) THEN
            IF ( nounit ) X(jx) = X(jx)/Ap(kk)
            temp = X(jx)
            ix = jx
            DO k = kk - 1, kk - j + 1, -1
              ix = ix - Incx
              X(ix) = X(ix) - temp*Ap(k)
            END DO
          END IF
          jx = jx - Incx
          kk = kk - j
        END DO
      END IF
    ELSE
      kk = 1
      IF ( Incx==1 ) THEN
        DO j = 1, N
          IF ( X(j)/=ZERO ) THEN
            IF ( nounit ) X(j) = X(j)/Ap(kk)
            temp = X(j)
            k = kk + 1
            DO i = j + 1, N
              X(i) = X(i) - temp*Ap(k)
              k = k + 1
            END DO
          END IF
          kk = kk + (N-j+1)
        END DO
      ELSE
        jx = kx
        DO j = 1, N
          IF ( X(jx)/=ZERO ) THEN
            IF ( nounit ) X(jx) = X(jx)/Ap(kk)
            temp = X(jx)
            ix = jx
            DO k = kk + 1, kk + N - j
              ix = ix + Incx
              X(ix) = X(ix) - temp*Ap(k)
            END DO
          END IF
          jx = jx + Incx
          kk = kk + (N-j+1)
        END DO
      END IF
    END IF
    !
    !        Form  x := inv( A' )*x.
    !
  ELSEIF ( LSAME(Uplo,'U') ) THEN
    kk = 1
    IF ( Incx==1 ) THEN
      DO j = 1, N
        temp = X(j)
        k = kk
        DO i = 1, j - 1
          temp = temp - Ap(k)*X(i)
          k = k + 1
        END DO
        IF ( nounit ) temp = temp/Ap(kk+j-1)
        X(j) = temp
        kk = kk + j
      END DO
    ELSE
      jx = kx
      DO j = 1, N
        temp = X(jx)
        ix = kx
        DO k = kk, kk + j - 2
          temp = temp - Ap(k)*X(ix)
          ix = ix + Incx
        END DO
        IF ( nounit ) temp = temp/Ap(kk+j-1)
        X(jx) = temp
        jx = jx + Incx
        kk = kk + j
      END DO
    END IF
  ELSE
    kk = (N*(N+1))/2
    IF ( Incx==1 ) THEN
      DO j = N, 1, -1
        temp = X(j)
        k = kk
        DO i = N, j + 1, -1
          temp = temp - Ap(k)*X(i)
          k = k - 1
        END DO
        IF ( nounit ) temp = temp/Ap(kk-N+j)
        X(j) = temp
        kk = kk - (N-j+1)
      END DO
    ELSE
      kx = kx + (N-1)*Incx
      jx = kx
      DO j = N, 1, -1
        temp = X(jx)
        ix = kx
        DO k = kk, kk - (N-(j+1)), -1
          temp = temp - Ap(k)*X(ix)
          ix = ix - Incx
        END DO
        IF ( nounit ) temp = temp/Ap(kk-N+j)
        X(jx) = temp
        jx = jx - Incx
        kk = kk - (N-j+1)
      END DO
    END IF
  END IF
  !
  !
  !     End of DTPSV .
  !
END SUBROUTINE DTPSV
