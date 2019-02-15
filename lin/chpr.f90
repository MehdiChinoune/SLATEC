!DECK CHPR
SUBROUTINE CHPR(Uplo,N,Alpha,X,Incx,Ap)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  CHPR
  !***PURPOSE  Perform the hermitian rank 1 operation.
  !***LIBRARY   SLATEC (BLAS)
  !***CATEGORY  D1B4
  !***TYPE      COMPLEX (CHPR-C)
  !***KEYWORDS  LEVEL 2 BLAS, LINEAR ALGEBRA
  !***AUTHOR  Dongarra, J. J., (ANL)
  !           Du Croz, J., (NAG)
  !           Hammarling, S., (NAG)
  !           Hanson, R. J., (SNLA)
  !***DESCRIPTION
  !
  !  CHPR    performs the hermitian rank 1 operation
  !
  !     A := alpha*x*conjg( x') + A,
  !
  !  where alpha is a real scalar, x is an n element vector and A is an
  !  n by n hermitian matrix, supplied in packed form.
  !
  !  Parameters
  !  ==========
  !
  !  UPLO   - CHARACTER*1.
  !           On entry, UPLO specifies whether the upper or lower
  !           triangular part of the matrix A is supplied in the packed
  !           array AP as follows:
  !
  !              UPLO = 'U' or 'u'   The upper triangular part of A is
  !                                  supplied in AP.
  !
  !              UPLO = 'L' or 'l'   The lower triangular part of A is
  !                                  supplied in AP.
  !
  !           Unchanged on exit.
  !
  !  N      - INTEGER.
  !           On entry, N specifies the order of the matrix A.
  !           N must be at least zero.
  !           Unchanged on exit.
  !
  !  ALPHA  - REAL            .
  !           On entry, ALPHA specifies the scalar alpha.
  !           Unchanged on exit.
  !
  !  X      - COMPLEX          array of dimension at least
  !           ( 1 + ( n - 1 )*abs( INCX ) ).
  !           Before entry, the incremented array X must contain the n
  !           element vector x.
  !           Unchanged on exit.
  !
  !  INCX   - INTEGER.
  !           On entry, INCX specifies the increment for the elements of
  !           X. INCX must not be zero.
  !           Unchanged on exit.
  !
  !  AP     - COMPLEX          array of DIMENSION at least
  !           ( ( n*( n + 1 ) )/2 ).
  !           Before entry with  UPLO = 'U' or 'u', the array AP must
  !           contain the upper triangular part of the hermitian matrix
  !           packed sequentially, column by column, so that AP( 1 )
  !           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
  !           and a( 2, 2 ) respectively, and so on. On exit, the array
  !           AP is overwritten by the upper triangular part of the
  !           updated matrix.
  !           Before entry with UPLO = 'L' or 'l', the array AP must
  !           contain the lower triangular part of the hermitian matrix
  !           packed sequentially, column by column, so that AP( 1 )
  !           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
  !           and a( 3, 1 ) respectively, and so on. On exit, the array
  !           AP is overwritten by the lower triangular part of the
  !           updated matrix.
  !           Note that the imaginary parts of the diagonal elements need
  !           not be set, they are assumed to be zero, and on exit they
  !           are set to zero.
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
  !***END PROLOGUE  CHPR
  !     .. Scalar Arguments ..
  REAL Alpha
  INTEGER Incx, N
  CHARACTER :: Uplo
  !     .. Array Arguments ..
  COMPLEX Ap(*), X(*)
  !     .. Parameters ..
  COMPLEX ZERO
  PARAMETER (ZERO=(0.0E+0,0.0E+0))
  !     .. Local Scalars ..
  COMPLEX temp
  INTEGER i, info, ix, j, jx, k, kk, kx
  !     .. External Functions ..
  LOGICAL LSAME
  EXTERNAL LSAME
  !     .. External Subroutines ..
  EXTERNAL XERBLA
  !     .. Intrinsic Functions ..
  INTRINSIC CONJG, REAL
  !***FIRST EXECUTABLE STATEMENT  CHPR
  !
  !     Test the input parameters.
  !
  info = 0
  IF ( .NOT.LSAME(Uplo,'U').AND..NOT.LSAME(Uplo,'L') ) THEN
    info = 1
  ELSEIF ( N<0 ) THEN
    info = 2
  ELSEIF ( Incx==0 ) THEN
    info = 5
  ENDIF
  IF ( info/=0 ) THEN
    CALL XERBLA('CHPR  ',info)
    RETURN
  ENDIF
  !
  !     Quick return if possible.
  !
  IF ( (N==0).OR.(Alpha==REAL(ZERO)) ) RETURN
  !
  !     Set the start point in X if the increment is not unity.
  !
  IF ( Incx<=0 ) THEN
    kx = 1 - (N-1)*Incx
  ELSEIF ( Incx/=1 ) THEN
    kx = 1
  ENDIF
  !
  !     Start the operations. In this version the elements of the array AP
  !     are accessed sequentially with one pass through AP.
  !
  kk = 1
  IF ( LSAME(Uplo,'U') ) THEN
    !
    !        Form  A  when upper triangle is stored in AP.
    !
    IF ( Incx==1 ) THEN
      DO j = 1, N
        IF ( X(j)/=ZERO ) THEN
          temp = Alpha*CONJG(X(j))
          k = kk
          DO i = 1, j - 1
            Ap(k) = Ap(k) + X(i)*temp
            k = k + 1
          ENDDO
          Ap(kk+j-1) = REAL(Ap(kk+j-1)) + REAL(X(j)*temp)
        ELSE
          Ap(kk+j-1) = REAL(Ap(kk+j-1))
        ENDIF
        kk = kk + j
      ENDDO
    ELSE
      jx = kx
      DO j = 1, N
        IF ( X(jx)/=ZERO ) THEN
          temp = Alpha*CONJG(X(jx))
          ix = kx
          DO k = kk, kk + j - 2
            Ap(k) = Ap(k) + X(ix)*temp
            ix = ix + Incx
          ENDDO
          Ap(kk+j-1) = REAL(Ap(kk+j-1)) + REAL(X(jx)*temp)
        ELSE
          Ap(kk+j-1) = REAL(Ap(kk+j-1))
        ENDIF
        jx = jx + Incx
        kk = kk + j
      ENDDO
    ENDIF
    !
    !        Form  A  when lower triangle is stored in AP.
    !
  ELSEIF ( Incx==1 ) THEN
    DO j = 1, N
      IF ( X(j)/=ZERO ) THEN
        temp = Alpha*CONJG(X(j))
        Ap(kk) = REAL(Ap(kk)) + REAL(temp*X(j))
        k = kk + 1
        DO i = j + 1, N
          Ap(k) = Ap(k) + X(i)*temp
          k = k + 1
        ENDDO
      ELSE
        Ap(kk) = REAL(Ap(kk))
      ENDIF
      kk = kk + N - j + 1
    ENDDO
  ELSE
    jx = kx
    DO j = 1, N
      IF ( X(jx)/=ZERO ) THEN
        temp = Alpha*CONJG(X(jx))
        Ap(kk) = REAL(Ap(kk)) + REAL(temp*X(jx))
        ix = jx
        DO k = kk + 1, kk + N - j
          ix = ix + Incx
          Ap(k) = Ap(k) + X(ix)*temp
        ENDDO
      ELSE
        Ap(kk) = REAL(Ap(kk))
      ENDIF
      jx = jx + Incx
      kk = kk + N - j + 1
    ENDDO
  ENDIF
  !
  !
  !     End of CHPR  .
  !
END SUBROUTINE CHPR
