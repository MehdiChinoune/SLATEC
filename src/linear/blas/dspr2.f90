!** DSPR2
SUBROUTINE DSPR2(Uplo,N,Alpha,X,Incx,Y,Incy,Ap)
  IMPLICIT NONE
  !>
  !***
  !  Perform the symmetric rank 2 operation.
  !***
  ! **Library:**   SLATEC (BLAS)
  !***
  ! **Category:**  D1B4
  !***
  ! **Type:**      DOUBLE PRECISION (SSPR2-S, DSPR2-D, CSPR2-C)
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
  !  DSPR2  performs the symmetric rank 2 operation
  !
  !     A := alpha*x*y' + alpha*y*x' + A,
  !
  !  where alpha is a scalar, x and y are n element vectors and A is an
  !  n by n symmetric matrix, supplied in packed form.
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
  !  AP     - DOUBLE PRECISION array of DIMENSION at least
  !           ( ( n*( n + 1 ) )/2 ).
  !           Before entry with  UPLO = 'U' or 'u', the array AP must
  !           contain the upper triangular part of the symmetric matrix
  !           packed sequentially, column by column, so that AP( 1 )
  !           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
  !           and a( 2, 2 ) respectively, and so on. On exit, the array
  !           AP is overwritten by the upper triangular part of the
  !           updated matrix.
  !           Before entry with UPLO = 'L' or 'l', the array AP must
  !           contain the lower triangular part of the symmetric matrix
  !           packed sequentially, column by column, so that AP( 1 )
  !           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
  !           and a( 3, 1 ) respectively, and so on. On exit, the array
  !           AP is overwritten by the lower triangular part of the
  !           updated matrix.
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

  !     .. Scalar Arguments ..
  REAL(8) :: Alpha
  INTEGER Incx, Incy, N
  CHARACTER :: Uplo
  !     .. Array Arguments ..
  REAL(8) :: Ap(*), X(*), Y(*)
  !     .. Parameters ..
  REAL(8), PARAMETER :: ZERO = 0.0D+0
  !     .. Local Scalars ..
  REAL(8) :: temp1, temp2
  INTEGER i, info, ix, iy, j, jx, jy, k, kk, kx, ky
  !     .. External Functions ..
  LOGICAL, EXTERNAL :: LSAME
  !     .. External Subroutines ..
  EXTERNAL :: XERBLA
  !* FIRST EXECUTABLE STATEMENT  DSPR2
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
  END IF
  IF ( info/=0 ) THEN
    CALL XERBLA('DSPR2 ',info)
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
  !     Start the operations. In this version the elements of the array AP
  !     are accessed sequentially with one pass through AP.
  !
  kk = 1
  IF ( LSAME(Uplo,'U') ) THEN
    !
    !        Form  A  when upper triangle is stored in AP.
    !
    IF ( (Incx==1).AND.(Incy==1) ) THEN
      DO j = 1, N
        IF ( (X(j)/=ZERO).OR.(Y(j)/=ZERO) ) THEN
          temp1 = Alpha*Y(j)
          temp2 = Alpha*X(j)
          k = kk
          DO i = 1, j
            Ap(k) = Ap(k) + X(i)*temp1 + Y(i)*temp2
            k = k + 1
          END DO
        END IF
        kk = kk + j
      END DO
    ELSE
      DO j = 1, N
        IF ( (X(jx)/=ZERO).OR.(Y(jy)/=ZERO) ) THEN
          temp1 = Alpha*Y(jy)
          temp2 = Alpha*X(jx)
          ix = kx
          iy = ky
          DO k = kk, kk + j - 1
            Ap(k) = Ap(k) + X(ix)*temp1 + Y(iy)*temp2
            ix = ix + Incx
            iy = iy + Incy
          END DO
        END IF
        jx = jx + Incx
        jy = jy + Incy
        kk = kk + j
      END DO
    END IF
    !
    !        Form  A  when lower triangle is stored in AP.
    !
  ELSEIF ( (Incx==1).AND.(Incy==1) ) THEN
    DO j = 1, N
      IF ( (X(j)/=ZERO).OR.(Y(j)/=ZERO) ) THEN
        temp1 = Alpha*Y(j)
        temp2 = Alpha*X(j)
        k = kk
        DO i = j, N
          Ap(k) = Ap(k) + X(i)*temp1 + Y(i)*temp2
          k = k + 1
        END DO
      END IF
      kk = kk + N - j + 1
    END DO
  ELSE
    DO j = 1, N
      IF ( (X(jx)/=ZERO).OR.(Y(jy)/=ZERO) ) THEN
        temp1 = Alpha*Y(jy)
        temp2 = Alpha*X(jx)
        ix = jx
        iy = jy
        DO k = kk, kk + N - j
          Ap(k) = Ap(k) + X(ix)*temp1 + Y(iy)*temp2
          ix = ix + Incx
          iy = iy + Incy
        END DO
      END IF
      jx = jx + Incx
      jy = jy + Incy
      kk = kk + N - j + 1
    END DO
  END IF
  !
  !
  !     End of DSPR2 .
  !
END SUBROUTINE DSPR2
