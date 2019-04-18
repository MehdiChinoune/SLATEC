!** DSPMV
SUBROUTINE DSPMV(Uplo,N,Alpha,Ap,X,Incx,Beta,Y,Incy)
  !>
  !  Perform the matrix-vector operation.
  !***
  ! **Library:**   SLATEC (BLAS)
  !***
  ! **Category:**  D1B4
  !***
  ! **Type:**      DOUBLE PRECISION (SSPMV-S, DSPMV-D, CSPMV-C)
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
  !  DSPMV  performs the matrix-vector operation
  !
  !     y := alpha*A*x + beta*y,
  !
  !  where alpha and beta are scalars, x and y are n element vectors and
  !  A is an n by n symmetric matrix, supplied in packed form.
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
  !  AP     - DOUBLE PRECISION array of DIMENSION at least
  !           ( ( n*( n + 1))/2).
  !           Before entry with UPLO = 'U' or 'u', the array AP must
  !           contain the upper triangular part of the symmetric matrix
  !           packed sequentially, column by column, so that AP( 1 )
  !           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
  !           and a( 2, 2 ) respectively, and so on.
  !           Before entry with UPLO = 'L' or 'l', the array AP must
  !           contain the lower triangular part of the symmetric matrix
  !           packed sequentially, column by column, so that AP( 1 )
  !           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
  !           and a( 3, 1 ) respectively, and so on.
  !           Unchanged on exit.
  !
  !  X      - DOUBLE PRECISION array of dimension at least
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
  !  BETA   - DOUBLE PRECISION.
  !           On entry, BETA specifies the scalar beta. When BETA is
  !           supplied as zero then Y need not be set on input.
  !           Unchanged on exit.
  !
  !  Y      - DOUBLE PRECISION array of dimension at least
  !           ( 1 + ( n - 1 )*abs( INCY ) ).
  !           Before entry, the incremented array Y must contain the n
  !           element vector y. On exit, Y is overwritten by the updated
  !           vector y.
  !
  !  INCY   - INTEGER.
  !           On entry, INCY specifies the increment for the elements of
  !           Y. INCY must not be zero.
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
  REAL(8) :: Alpha, Beta
  INTEGER Incx, Incy, N
  CHARACTER :: Uplo
  !     .. Array Arguments ..
  REAL(8) :: Ap(*), X(*), Y(*)
  !     .. Parameters ..
  REAL(8), PARAMETER :: ONE = 1.0D+0, ZERO = 0.0D+0
  !     .. Local Scalars ..
  REAL(8) :: temp1, temp2
  INTEGER i, info, ix, iy, j, jx, jy, k, kk, kx, ky
  !* FIRST EXECUTABLE STATEMENT  DSPMV
  !
  !     Test the input parameters.
  !
  info = 0
  IF ( .NOT.LSAME(Uplo,'U').AND..NOT.LSAME(Uplo,'L') ) THEN
    info = 1
  ELSEIF ( N<0 ) THEN
    info = 2
  ELSEIF ( Incx==0 ) THEN
    info = 6
  ELSEIF ( Incy==0 ) THEN
    info = 9
  END IF
  IF ( info/=0 ) THEN
    CALL XERBLA('DSPMV ',info)
    RETURN
  END IF
  !
  !     Quick return if possible.
  !
  IF ( (N==0).OR.((Alpha==ZERO).AND.(Beta==ONE)) ) RETURN
  !
  !     Set up the start points in  X  and  Y.
  !
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
  !
  !     Start the operations. In this version the elements of the array AP
  !     are accessed sequentially with one pass through AP.
  !
  !     First form  y := beta*y.
  !
  IF ( Beta/=ONE ) THEN
    IF ( Incy/=1 ) THEN
      iy = ky
      IF ( Beta==ZERO ) THEN
        DO i = 1, N
          Y(iy) = ZERO
          iy = iy + Incy
        END DO
      ELSE
        DO i = 1, N
          Y(iy) = Beta*Y(iy)
          iy = iy + Incy
        END DO
      END IF
    ELSEIF ( Beta==ZERO ) THEN
      DO i = 1, N
        Y(i) = ZERO
      END DO
    ELSE
      DO i = 1, N
        Y(i) = Beta*Y(i)
      END DO
    END IF
  END IF
  IF ( Alpha==ZERO ) RETURN
  kk = 1
  IF ( LSAME(Uplo,'U') ) THEN
    !
    !        Form  y  when AP contains the upper triangle.
    !
    IF ( (Incx==1).AND.(Incy==1) ) THEN
      DO j = 1, N
        temp1 = Alpha*X(j)
        temp2 = ZERO
        k = kk
        DO i = 1, j - 1
          Y(i) = Y(i) + temp1*Ap(k)
          temp2 = temp2 + Ap(k)*X(i)
          k = k + 1
        END DO
        Y(j) = Y(j) + temp1*Ap(kk+j-1) + Alpha*temp2
        kk = kk + j
      END DO
    ELSE
      jx = kx
      jy = ky
      DO j = 1, N
        temp1 = Alpha*X(jx)
        temp2 = ZERO
        ix = kx
        iy = ky
        DO k = kk, kk + j - 2
          Y(iy) = Y(iy) + temp1*Ap(k)
          temp2 = temp2 + Ap(k)*X(ix)
          ix = ix + Incx
          iy = iy + Incy
        END DO
        Y(jy) = Y(jy) + temp1*Ap(kk+j-1) + Alpha*temp2
        jx = jx + Incx
        jy = jy + Incy
        kk = kk + j
      END DO
    END IF
    !
    !        Form  y  when AP contains the lower triangle.
    !
  ELSEIF ( (Incx==1).AND.(Incy==1) ) THEN
    DO j = 1, N
      temp1 = Alpha*X(j)
      temp2 = ZERO
      Y(j) = Y(j) + temp1*Ap(kk)
      k = kk + 1
      DO i = j + 1, N
        Y(i) = Y(i) + temp1*Ap(k)
        temp2 = temp2 + Ap(k)*X(i)
        k = k + 1
      END DO
      Y(j) = Y(j) + Alpha*temp2
      kk = kk + (N-j+1)
    END DO
  ELSE
    jx = kx
    jy = ky
    DO j = 1, N
      temp1 = Alpha*X(jx)
      temp2 = ZERO
      Y(jy) = Y(jy) + temp1*Ap(kk)
      ix = jx
      iy = jy
      DO k = kk + 1, kk + N - j
        ix = ix + Incx
        iy = iy + Incy
        Y(iy) = Y(iy) + temp1*Ap(k)
        temp2 = temp2 + Ap(k)*X(ix)
      END DO
      Y(jy) = Y(jy) + Alpha*temp2
      jx = jx + Incx
      jy = jy + Incy
      kk = kk + (N-j+1)
    END DO
  END IF
  !
  !
  !     End of DSPMV .
  !
END SUBROUTINE DSPMV
