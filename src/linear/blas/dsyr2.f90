!** DSYR2
SUBROUTINE DSYR2(Uplo,N,Alpha,X,Incx,Y,Incy,A,Lda)
  !>
  !  Perform the symmetric rank 2 operation.
  !***
  ! **Library:**   SLATEC (BLAS)
  !***
  ! **Category:**  D1B4
  !***
  ! **Type:**      DOUBLE PRECISION (SSYR2-S, DSYR2-D, CSYR2-C)
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
  !  DSYR2  performs the symmetric rank 2 operation
  !
  !     A := alpha*x*y' + alpha*y*x' + A,
  !
  !  where alpha is a scalar, x and y are n element vectors and A is an n
  !  by n symmetric matrix.
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
  !  ALPHA  - DOUBLE PRECISION.
  !           On entry, ALPHA specifies the scalar alpha.
  !           Unchanged on exit.
  !
  !  X      - DOUBLE PRECISION array of dimension at least
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
  REAL(8) :: Alpha
  INTEGER Incx, Incy, Lda, N
  CHARACTER :: Uplo
  !     .. Array Arguments ..
  REAL(8) :: A(Lda,*), X(*), Y(*)
  !     .. Parameters ..
  REAL(8), PARAMETER :: ZERO = 0.0D+0
  !     .. Local Scalars ..
  REAL(8) :: temp1, temp2
  INTEGER i, info, ix, iy, j, jx, jy, kx, ky
  !     .. Intrinsic Functions ..
  INTRINSIC MAX
  !* FIRST EXECUTABLE STATEMENT  DSYR2
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
  ELSEIF ( Incy==0 ) THEN
    info = 7
  ELSEIF ( Lda<MAX(1,N) ) THEN
    info = 9
  END IF
  IF ( info/=0 ) THEN
    CALL XERBLA('DSYR2 ',info)
    RETURN
  END IF
  !
  !     Quick return if possible.
  !
  IF ( (N==0).OR.(Alpha==ZERO) ) RETURN
  !
  !     Set up the start points in X and Y if the increments are not both
  !     unity.
  !
  IF ( (Incx/=1).OR.(Incy/=1) ) THEN
    IF ( Incx>0 ) THEN
      kx = 1
    ELSE
      kx = 1 - (N-1)*Incx
    END IF
    IF ( Incy>0 ) THEN
      ky = 1
    ELSE
      ky = 1 - (N-1)*Incy
    END IF
    jx = kx
    jy = ky
  END IF
  !
  !     Start the operations. In this version the elements of A are
  !     accessed sequentially with one pass through the triangular part
  !     of A.
  !
  IF ( LSAME(Uplo,'U') ) THEN
    !
    !        Form  A  when A is stored in the upper triangle.
    !
    IF ( (Incx==1).AND.(Incy==1) ) THEN
      DO j = 1, N
        IF ( (X(j)/=ZERO).OR.(Y(j)/=ZERO) ) THEN
          temp1 = Alpha*Y(j)
          temp2 = Alpha*X(j)
          DO i = 1, j
            A(i,j) = A(i,j) + X(i)*temp1 + Y(i)*temp2
          END DO
        END IF
      END DO
    ELSE
      DO j = 1, N
        IF ( (X(jx)/=ZERO).OR.(Y(jy)/=ZERO) ) THEN
          temp1 = Alpha*Y(jy)
          temp2 = Alpha*X(jx)
          ix = kx
          iy = ky
          DO i = 1, j
            A(i,j) = A(i,j) + X(ix)*temp1 + Y(iy)*temp2
            ix = ix + Incx
            iy = iy + Incy
          END DO
        END IF
        jx = jx + Incx
        jy = jy + Incy
      END DO
    END IF
    !
    !        Form  A  when A is stored in the lower triangle.
    !
  ELSEIF ( (Incx==1).AND.(Incy==1) ) THEN
    DO j = 1, N
      IF ( (X(j)/=ZERO).OR.(Y(j)/=ZERO) ) THEN
        temp1 = Alpha*Y(j)
        temp2 = Alpha*X(j)
        DO i = j, N
          A(i,j) = A(i,j) + X(i)*temp1 + Y(i)*temp2
        END DO
      END IF
    END DO
  ELSE
    DO j = 1, N
      IF ( (X(jx)/=ZERO).OR.(Y(jy)/=ZERO) ) THEN
        temp1 = Alpha*Y(jy)
        temp2 = Alpha*X(jx)
        ix = jx
        iy = jy
        DO i = j, N
          A(i,j) = A(i,j) + X(ix)*temp1 + Y(iy)*temp2
          ix = ix + Incx
          iy = iy + Incy
        END DO
      END IF
      jx = jx + Incx
      jy = jy + Incy
    END DO
  END IF
  !
  !
  !     End of DSYR2 .
  !
END SUBROUTINE DSYR2
