!** DTPMV
SUBROUTINE DTPMV(Uplo,Trans,Diag,N,Ap,X,Incx)
  !>
  !  Perform one of the matrix-vector operations.
  !***
  ! **Library:**   SLATEC (BLAS)
  !***
  ! **Category:**  D1B4
  !***
  ! **Type:**      DOUBLE PRECISION (STPMV-S, DTPMV-D, CTPMV-C)
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
  !  DTPMV  performs one of the matrix-vector operations
  !
  !     x := A*x,   or   x := A'*x,
  !
  !  where x is an n element vector and  A is an n by n unit, or non-unit,
  !  upper or lower triangular matrix, supplied in packed form.
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
  !           On entry, TRANS specifies the operation to be performed as
  !           follows:
  !
  !              TRANS = 'N' or 'n'   x := A*x.
  !
  !              TRANS = 'T' or 't'   x := A'*x.
  !
  !              TRANS = 'C' or 'c'   x := A'*x.
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
  !           element vector x. On exit, X is overwritten with the
  !           transformed vector x.
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
  !* FIRST EXECUTABLE STATEMENT  DTPMV
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
    CALL XERBLA('DTPMV ',info)
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
    !        Form  x:= A*x.
    !
    IF ( LSAME(Uplo,'U') ) THEN
      kk = 1
      IF ( Incx==1 ) THEN
        DO j = 1, N
          IF ( X(j)/=ZERO ) THEN
            temp = X(j)
            k = kk
            DO i = 1, j - 1
              X(i) = X(i) + temp*Ap(k)
              k = k + 1
            END DO
            IF ( nounit ) X(j) = X(j)*Ap(kk+j-1)
          END IF
          kk = kk + j
        END DO
      ELSE
        jx = kx
        DO j = 1, N
          IF ( X(jx)/=ZERO ) THEN
            temp = X(jx)
            ix = kx
            DO k = kk, kk + j - 2
              X(ix) = X(ix) + temp*Ap(k)
              ix = ix + Incx
            END DO
            IF ( nounit ) X(jx) = X(jx)*Ap(kk+j-1)
          END IF
          jx = jx + Incx
          kk = kk + j
        END DO
      END IF
    ELSE
      kk = (N*(N+1))/2
      IF ( Incx==1 ) THEN
        DO j = N, 1, -1
          IF ( X(j)/=ZERO ) THEN
            temp = X(j)
            k = kk
            DO i = N, j + 1, -1
              X(i) = X(i) + temp*Ap(k)
              k = k - 1
            END DO
            IF ( nounit ) X(j) = X(j)*Ap(kk-N+j)
          END IF
          kk = kk - (N-j+1)
        END DO
      ELSE
        kx = kx + (N-1)*Incx
        jx = kx
        DO j = N, 1, -1
          IF ( X(jx)/=ZERO ) THEN
            temp = X(jx)
            ix = kx
            DO k = kk, kk - (N-(j+1)), -1
              X(ix) = X(ix) + temp*Ap(k)
              ix = ix - Incx
            END DO
            IF ( nounit ) X(jx) = X(jx)*Ap(kk-N+j)
          END IF
          jx = jx - Incx
          kk = kk - (N-j+1)
        END DO
      END IF
    END IF
    !
    !        Form  x := A'*x.
    !
  ELSEIF ( LSAME(Uplo,'U') ) THEN
    kk = (N*(N+1))/2
    IF ( Incx==1 ) THEN
      DO j = N, 1, -1
        temp = X(j)
        IF ( nounit ) temp = temp*Ap(kk)
        k = kk - 1
        DO i = j - 1, 1, -1
          temp = temp + Ap(k)*X(i)
          k = k - 1
        END DO
        X(j) = temp
        kk = kk - j
      END DO
    ELSE
      jx = kx + (N-1)*Incx
      DO j = N, 1, -1
        temp = X(jx)
        ix = jx
        IF ( nounit ) temp = temp*Ap(kk)
        DO k = kk - 1, kk - j + 1, -1
          ix = ix - Incx
          temp = temp + Ap(k)*X(ix)
        END DO
        X(jx) = temp
        jx = jx - Incx
        kk = kk - j
      END DO
    END IF
  ELSE
    kk = 1
    IF ( Incx==1 ) THEN
      DO j = 1, N
        temp = X(j)
        IF ( nounit ) temp = temp*Ap(kk)
        k = kk + 1
        DO i = j + 1, N
          temp = temp + Ap(k)*X(i)
          k = k + 1
        END DO
        X(j) = temp
        kk = kk + (N-j+1)
      END DO
    ELSE
      jx = kx
      DO j = 1, N
        temp = X(jx)
        ix = jx
        IF ( nounit ) temp = temp*Ap(kk)
        DO k = kk + 1, kk + N - j
          ix = ix + Incx
          temp = temp + Ap(k)*X(ix)
        END DO
        X(jx) = temp
        jx = jx + Incx
        kk = kk + (N-j+1)
      END DO
    END IF
  END IF
  !
  !
  !     End of DTPMV .
  !
END SUBROUTINE DTPMV
