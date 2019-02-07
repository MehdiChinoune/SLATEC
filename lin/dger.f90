!*==DGER.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK DGER
SUBROUTINE DGER(M,N,Alpha,X,Incx,Y,Incy,A,Lda)
  IMPLICIT NONE
  !*--DGER5
  !***BEGIN PROLOGUE  DGER
  !***PURPOSE  Perform the rank 1 operation.
  !***LIBRARY   SLATEC (BLAS)
  !***CATEGORY  D1B4
  !***TYPE      DOUBLE PRECISION (DGER-D)
  !***KEYWORDS  LEVEL 2 BLAS, LINEAR ALGEBRA
  !***AUTHOR  Dongarra, J. J., (ANL)
  !           Du Croz, J., (NAG)
  !           Hammarling, S., (NAG)
  !           Hanson, R. J., (SNLA)
  !***DESCRIPTION
  !
  !  DGER   performs the rank 1 operation
  !
  !     A := alpha*x*y' + A,
  !
  !  where alpha is a scalar, x is an m element vector, y is an n element
  !  vector and A is an m by n matrix.
  !
  !  Parameters
  !  ==========
  !
  !  M      - INTEGER.
  !           On entry, M specifies the number of rows of the matrix A.
  !           M must be at least zero.
  !           Unchanged on exit.
  !
  !  N      - INTEGER.
  !           On entry, N specifies the number of columns of the matrix A.
  !           N must be at least zero.
  !           Unchanged on exit.
  !
  !  ALPHA  - DOUBLE PRECISION.
  !           On entry, ALPHA specifies the scalar alpha.
  !           Unchanged on exit.
  !
  !  X      - DOUBLE PRECISION array of dimension at least
  !           ( 1 + ( m - 1)*abs( INCX)).
  !           Before entry, the incremented array X must contain the m
  !           element vector x.
  !           Unchanged on exit.
  !
  !  INCX   - INTEGER.
  !           On entry, INCX specifies the increment for the elements of
  !           X. INCX must not be zero.
  !           Unchanged on exit.
  !
  !  Y      - DOUBLE PRECISION array of dimension at least
  !           ( 1 + ( n - 1 )*abs( INCY ) ).
  !           Before entry, the incremented array Y must contain the n
  !           element vector y.
  !           Unchanged on exit.
  !
  !  INCY   - INTEGER.
  !           On entry, INCY specifies the increment for the elements of
  !           Y. INCY must not be zero.
  !           Unchanged on exit.
  !
  !  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
  !           Before entry, the leading m by n part of the array A must
  !           contain the matrix of coefficients. On exit, A is
  !           overwritten by the updated matrix.
  !
  !  LDA    - INTEGER.
  !           On entry, LDA specifies the first dimension of A as declared
  !           in the calling (sub) program. LDA must be at least
  !           max( 1, m ).
  !           Unchanged on exit.
  !
  !***REFERENCES  Dongarra, J. J., Du Croz, J., Hammarling, S., and
  !                 Hanson, R. J.  An extended set of Fortran basic linear
  !                 algebra subprograms.  ACM TOMS, Vol. 14, No. 1,
  !                 pp. 1-17, March 1988.
  !***ROUTINES CALLED  XERBLA
  !***REVISION HISTORY  (YYMMDD)
  !   861022  DATE WRITTEN
  !   910605  Modified to meet SLATEC prologue standards.  Only comment
  !           lines were modified.  (BKS)
  !***END PROLOGUE  DGER
  !     .. Scalar Arguments ..
  REAL(8) :: Alpha
  INTEGER Incx, Incy, Lda, M, N
  !     .. Array Arguments ..
  REAL(8) :: A(Lda,*), X(*), Y(*)
  !     .. Parameters ..
  REAL(8) :: ZERO
  PARAMETER (ZERO=0.0D+0)
  !     .. Local Scalars ..
  REAL(8) :: temp
  INTEGER i, info, ix, j, jy, kx
  !     .. External Subroutines ..
  EXTERNAL XERBLA
  !     .. Intrinsic Functions ..
  INTRINSIC MAX
  !***FIRST EXECUTABLE STATEMENT  DGER
  !
  !     Test the input parameters.
  !
  info = 0
  IF ( M<0 ) THEN
    info = 1
  ELSEIF ( N<0 ) THEN
    info = 2
  ELSEIF ( Incx==0 ) THEN
    info = 5
  ELSEIF ( Incy==0 ) THEN
    info = 7
  ELSEIF ( Lda<MAX(1,M) ) THEN
    info = 9
  ENDIF
  IF ( info/=0 ) THEN
    CALL XERBLA('DGER  ',info)
    RETURN
  ENDIF
  !
  !     Quick return if possible.
  !
  IF ( (M==0).OR.(N==0).OR.(Alpha==ZERO) ) RETURN
  !
  !     Start the operations. In this version the elements of A are
  !     accessed sequentially with one pass through A.
  !
  IF ( Incy>0 ) THEN
    jy = 1
  ELSE
    jy = 1 - (N-1)*Incy
  ENDIF
  IF ( Incx==1 ) THEN
    DO j = 1, N
      IF ( Y(jy)/=ZERO ) THEN
        temp = Alpha*Y(jy)
        DO i = 1, M
          A(i,j) = A(i,j) + X(i)*temp
        ENDDO
      ENDIF
      jy = jy + Incy
    ENDDO
  ELSE
    IF ( Incx>0 ) THEN
      kx = 1
    ELSE
      kx = 1 - (M-1)*Incx
    ENDIF
    DO j = 1, N
      IF ( Y(jy)/=ZERO ) THEN
        temp = Alpha*Y(jy)
        ix = kx
        DO i = 1, M
          A(i,j) = A(i,j) + X(ix)*temp
          ix = ix + Incx
        ENDDO
      ENDIF
      jy = jy + Incy
    ENDDO
  ENDIF
  !
  !
  !     End of DGER  .
  !
END SUBROUTINE DGER
