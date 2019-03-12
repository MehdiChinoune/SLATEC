!DECK SSYR
SUBROUTINE SSYR(Uplo,N,Alpha,X,Incx,A,Lda)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  SSYR
  !***PURPOSE  Perform symmetric rank 1 update of a real symmetric matrix.
  !***LIBRARY   SLATEC (BLAS)
  !***CATEGORY  D1B4
  !***TYPE      SINGLE PRECISION (SSYR-S)
  !***KEYWORDS  LEVEL 2 BLAS, LINEAR ALGEBRA
  !***AUTHOR  Dongarra, J. J., (ANL)
  !           Du Croz, J., (NAG)
  !           Hammarling, S., (NAG)
  !           Hanson, R. J., (SNLA)
  !***DESCRIPTION
  !
  !  SSYR   performs the symmetric rank 1 operation
  !
  !     A := alpha*x*x' + A,
  !
  !  where alpha is a real scalar, x is an n element vector and A is an
  !  n by n symmetric matrix.
  !
  !  Parameters
  !  ==========
  !
  !  UPLO   - CHARACTER*1.
  !           On entry, UPLO specifies whether the upper or lower
  !           triangular part of the array A is to be referenced as
  !           follows:
  !
  !              UPLO = 'U' or 'u'   Only the upper triangular part of A
  !                                  is to be referenced.
  !
  !              UPLO = 'L' or 'l'   Only the lower triangular part of A
  !                                  is to be referenced.
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
  !  X      - REAL             array of dimension at least
  !           ( 1 + ( n - 1)*abs( INCX)).
  !           Before entry, the incremented array X must contain the n
  !           element vector x.
  !           Unchanged on exit.
  !
  !  INCX   - INTEGER.
  !           On entry, INCX specifies the increment for the elements of
  !           X. INCX must not be zero.
  !           Unchanged on exit.
  !
  !  A      - REAL             array of DIMENSION ( LDA, n ).
  !           Before entry with  UPLO = 'U' or 'u', the leading n by n
  !           upper triangular part of the array A must contain the upper
  !           triangular part of the symmetric matrix and the strictly
  !           lower triangular part of A is not referenced. On exit, the
  !           upper triangular part of the array A is overwritten by the
  !           upper triangular part of the updated matrix.
  !           Before entry with UPLO = 'L' or 'l', the leading n by n
  !           lower triangular part of the array A must contain the lower
  !           triangular part of the symmetric matrix and the strictly
  !           upper triangular part of A is not referenced. On exit, the
  !           lower triangular part of the array A is overwritten by the
  !           lower triangular part of the updated matrix.
  !
  !  LDA    - INTEGER.
  !           On entry, LDA specifies the first dimension of A as declared
  !           in the calling (sub) program. LDA must be at least
  !           max( 1, n ).
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
  !***END PROLOGUE  SSYR
  !     .. Scalar Arguments ..
  REAL Alpha
  INTEGER Incx, Lda, N
  CHARACTER :: Uplo
  !     .. Array Arguments ..
  REAL A(Lda,*), X(*)
  !     .. Parameters ..
  REAL ZERO
  PARAMETER (ZERO=0.0E+0)
  !     .. Local Scalars ..
  REAL temp
  INTEGER i, info, ix, j, jx, kx
  !     .. External Functions ..
  LOGICAL LSAME
  EXTERNAL LSAME
  !     .. External Subroutines ..
  EXTERNAL XERBLA
  !     .. Intrinsic Functions ..
  INTRINSIC MAX
  !***FIRST EXECUTABLE STATEMENT  SSYR
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
  ELSEIF ( Lda<MAX(1,N) ) THEN
    info = 7
  ENDIF
  IF ( info/=0 ) THEN
    CALL XERBLA('SSYR  ',info)
    RETURN
  ENDIF
  !
  !     Quick return if possible.
  !
  IF ( (N==0).OR.(Alpha==ZERO) ) RETURN
  !
  !     Set the start point in X if the increment is not unity.
  !
  IF ( Incx<=0 ) THEN
    kx = 1 - (N-1)*Incx
  ELSEIF ( Incx/=1 ) THEN
    kx = 1
  ENDIF
  !
  !     Start the operations. In this version the elements of A are
  !     accessed sequentially with one pass through the triangular part
  !     of A.
  !
  IF ( LSAME(Uplo,'U') ) THEN
    !
    !        Form  A  when A is stored in upper triangle.
    !
    IF ( Incx==1 ) THEN
      DO j = 1, N
        IF ( X(j)/=ZERO ) THEN
          temp = Alpha*X(j)
          DO i = 1, j
            A(i,j) = A(i,j) + X(i)*temp
          ENDDO
        ENDIF
      ENDDO
    ELSE
      jx = kx
      DO j = 1, N
        IF ( X(jx)/=ZERO ) THEN
          temp = Alpha*X(jx)
          ix = kx
          DO i = 1, j
            A(i,j) = A(i,j) + X(ix)*temp
            ix = ix + Incx
          ENDDO
        ENDIF
        jx = jx + Incx
      ENDDO
    ENDIF
    !
    !        Form  A  when A is stored in lower triangle.
    !
  ELSEIF ( Incx==1 ) THEN
    DO j = 1, N
      IF ( X(j)/=ZERO ) THEN
        temp = Alpha*X(j)
        DO i = j, N
          A(i,j) = A(i,j) + X(i)*temp
        ENDDO
      ENDIF
    ENDDO
  ELSE
    jx = kx
    DO j = 1, N
      IF ( X(jx)/=ZERO ) THEN
        temp = Alpha*X(jx)
        ix = jx
        DO i = j, N
          A(i,j) = A(i,j) + X(ix)*temp
          ix = ix + Incx
        ENDDO
      ENDIF
      jx = jx + Incx
    ENDDO
  ENDIF
  !
  !
  !     End of SSYR  .
  !
END SUBROUTINE SSYR
