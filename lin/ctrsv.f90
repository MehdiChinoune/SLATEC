!DECK CTRSV
SUBROUTINE CTRSV(Uplo,Trans,Diag,N,A,Lda,X,Incx)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  CTRSV
  !***PURPOSE  Solve a complex triangular system of equations.
  !***LIBRARY   SLATEC (BLAS)
  !***CATEGORY  D1B4
  !***TYPE      COMPLEX (STRSV-S, DTRSV-D, CTRSV-C)
  !***KEYWORDS  LEVEL 2 BLAS, LINEAR ALGEBRA
  !***AUTHOR  Dongarra, J. J., (ANL)
  !           Du Croz, J., (NAG)
  !           Hammarling, S., (NAG)
  !           Hanson, R. J., (SNLA)
  !***DESCRIPTION
  !
  !  CTRSV  solves one of the systems of equations
  !
  !     A*x = b,   or   A'*x = b,   or   conjg( A')*x = b,
  !
  !  where b and x are n element vectors and A is an n by n unit, or
  !  non-unit, upper or lower triangular matrix.
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
  !              TRANS = 'C' or 'c'   conjg( A' )*x = b.
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
  !  A      - COMPLEX          array of DIMENSION ( LDA, n ).
  !           Before entry with  UPLO = 'U' or 'u', the leading n by n
  !           upper triangular part of the array A must contain the upper
  !           triangular matrix and the strictly lower triangular part of
  !           A is not referenced.
  !           Before entry with UPLO = 'L' or 'l', the leading n by n
  !           lower triangular part of the array A must contain the lower
  !           triangular matrix and the strictly upper triangular part of
  !           A is not referenced.
  !           Note that when  DIAG = 'U' or 'u', the diagonal elements of
  !           A are not referenced either, but are assumed to be unity.
  !           Unchanged on exit.
  !
  !  LDA    - INTEGER.
  !           On entry, LDA specifies the first dimension of A as declared
  !           in the calling (sub) program. LDA must be at least
  !           max( 1, n ).
  !           Unchanged on exit.
  !
  !  X      - COMPLEX          array of dimension at least
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
  !***REFERENCES  Dongarra, J. J., Du Croz, J., Hammarling, S., and
  !                 Hanson, R. J.  An extended set of Fortran basic linear
  !                 algebra subprograms.  ACM TOMS, Vol. 14, No. 1,
  !                 pp. 1-17, March 1988.
  !***ROUTINES CALLED  LSAME, XERBLA
  !***REVISION HISTORY  (YYMMDD)
  !   861022  DATE WRITTEN
  !   910605  Modified to meet SLATEC prologue standards.  Only comment
  !           lines were modified.  (BKS)
  !***END PROLOGUE  CTRSV
  !     .. Scalar Arguments ..
  INTEGER Incx, Lda, N
  CHARACTER :: Diag, Trans, Uplo
  !     .. Array Arguments ..
  COMPLEX A(Lda,*), X(*)
  !     .. Parameters ..
  COMPLEX ZERO
  PARAMETER (ZERO=(0.0E+0,0.0E+0))
  !     .. Local Scalars ..
  COMPLEX temp
  INTEGER i, info, ix, j, jx, kx
  LOGICAL noconj, nounit
  !     .. External Functions ..
  LOGICAL LSAME
  EXTERNAL LSAME
  !     .. External Subroutines ..
  EXTERNAL XERBLA
  !     .. Intrinsic Functions ..
  INTRINSIC CONJG, MAX
  !***FIRST EXECUTABLE STATEMENT  CTRSV
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
  ELSEIF ( Lda<MAX(1,N) ) THEN
    info = 6
  ELSEIF ( Incx==0 ) THEN
    info = 8
  ENDIF
  IF ( info/=0 ) THEN
    CALL XERBLA('CTRSV ',info)
    RETURN
  ENDIF
  !
  !     Quick return if possible.
  !
  IF ( N==0 ) RETURN
  !
  noconj = LSAME(Trans,'T')
  nounit = LSAME(Diag,'N')
  !
  !     Set up the start point in X if the increment is not unity. This
  !     will be  ( N - 1 )*INCX  too small for descending loops.
  !
  IF ( Incx<=0 ) THEN
    kx = 1 - (N-1)*Incx
  ELSEIF ( Incx/=1 ) THEN
    kx = 1
  ENDIF
  !
  !     Start the operations. In this version the elements of A are
  !     accessed sequentially with one pass through A.
  !
  IF ( LSAME(Trans,'N') ) THEN
    !
    !        Form  x := inv( A )*x.
    !
    IF ( LSAME(Uplo,'U') ) THEN
      IF ( Incx==1 ) THEN
        DO j = N, 1, -1
          IF ( X(j)/=ZERO ) THEN
            IF ( nounit ) X(j) = X(j)/A(j,j)
            temp = X(j)
            DO i = j - 1, 1, -1
              X(i) = X(i) - temp*A(i,j)
            ENDDO
          ENDIF
        ENDDO
      ELSE
        jx = kx + (N-1)*Incx
        DO j = N, 1, -1
          IF ( X(jx)/=ZERO ) THEN
            IF ( nounit ) X(jx) = X(jx)/A(j,j)
            temp = X(jx)
            ix = jx
            DO i = j - 1, 1, -1
              ix = ix - Incx
              X(ix) = X(ix) - temp*A(i,j)
            ENDDO
          ENDIF
          jx = jx - Incx
        ENDDO
      ENDIF
    ELSEIF ( Incx==1 ) THEN
      DO j = 1, N
        IF ( X(j)/=ZERO ) THEN
          IF ( nounit ) X(j) = X(j)/A(j,j)
          temp = X(j)
          DO i = j + 1, N
            X(i) = X(i) - temp*A(i,j)
          ENDDO
        ENDIF
      ENDDO
    ELSE
      jx = kx
      DO j = 1, N
        IF ( X(jx)/=ZERO ) THEN
          IF ( nounit ) X(jx) = X(jx)/A(j,j)
          temp = X(jx)
          ix = jx
          DO i = j + 1, N
            ix = ix + Incx
            X(ix) = X(ix) - temp*A(i,j)
          ENDDO
        ENDIF
        jx = jx + Incx
      ENDDO
    ENDIF
    !
    !        Form  x := inv( A' )*x  or  x := inv( conjg( A' ) )*x.
    !
  ELSEIF ( LSAME(Uplo,'U') ) THEN
    IF ( Incx==1 ) THEN
      DO j = 1, N
        temp = X(j)
        IF ( noconj ) THEN
          DO i = 1, j - 1
            temp = temp - A(i,j)*X(i)
          ENDDO
          IF ( nounit ) temp = temp/A(j,j)
        ELSE
          DO i = 1, j - 1
            temp = temp - CONJG(A(i,j))*X(i)
          ENDDO
          IF ( nounit ) temp = temp/CONJG(A(j,j))
        ENDIF
        X(j) = temp
      ENDDO
    ELSE
      jx = kx
      DO j = 1, N
        ix = kx
        temp = X(jx)
        IF ( noconj ) THEN
          DO i = 1, j - 1
            temp = temp - A(i,j)*X(ix)
            ix = ix + Incx
          ENDDO
          IF ( nounit ) temp = temp/A(j,j)
        ELSE
          DO i = 1, j - 1
            temp = temp - CONJG(A(i,j))*X(ix)
            ix = ix + Incx
          ENDDO
          IF ( nounit ) temp = temp/CONJG(A(j,j))
        ENDIF
        X(jx) = temp
        jx = jx + Incx
      ENDDO
    ENDIF
  ELSEIF ( Incx==1 ) THEN
    DO j = N, 1, -1
      temp = X(j)
      IF ( noconj ) THEN
        DO i = N, j + 1, -1
          temp = temp - A(i,j)*X(i)
        ENDDO
        IF ( nounit ) temp = temp/A(j,j)
      ELSE
        DO i = N, j + 1, -1
          temp = temp - CONJG(A(i,j))*X(i)
        ENDDO
        IF ( nounit ) temp = temp/CONJG(A(j,j))
      ENDIF
      X(j) = temp
    ENDDO
  ELSE
    kx = kx + (N-1)*Incx
    jx = kx
    DO j = N, 1, -1
      ix = kx
      temp = X(jx)
      IF ( noconj ) THEN
        DO i = N, j + 1, -1
          temp = temp - A(i,j)*X(ix)
          ix = ix - Incx
        ENDDO
        IF ( nounit ) temp = temp/A(j,j)
      ELSE
        DO i = N, j + 1, -1
          temp = temp - CONJG(A(i,j))*X(ix)
          ix = ix - Incx
        ENDDO
        IF ( nounit ) temp = temp/CONJG(A(j,j))
      ENDIF
      X(jx) = temp
      jx = jx - Incx
    ENDDO
  ENDIF
  !
  !
  !     End of CTRSV .
  !
END SUBROUTINE CTRSV
