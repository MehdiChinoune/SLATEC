!** CHPR2
SUBROUTINE CHPR2(Uplo,N,Alpha,X,Incx,Y,Incy,Ap)
  !>
  !***
  !  Perform the hermitian rank 2 operation.
  !***
  ! **Library:**   SLATEC (BLAS)
  !***
  ! **Category:**  D1B4
  !***
  ! **Type:**      COMPLEX (SHPR2-S, DHPR2-D, CHPR2-C)
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
  !  CHPR2  performs the hermitian rank 2 operation
  !
  !     A := alpha*x*conjg( y') + conjg( alpha)*y*conjg( x') + A,
  !
  !  where alpha is a scalar, x and y are n element vectors and A is an
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
  !  ALPHA  - COMPLEX         .
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
  !  Y      - COMPLEX          array of dimension at least
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
  COMPLEX Alpha
  INTEGER Incx, Incy, N
  CHARACTER :: Uplo
  !     .. Array Arguments ..
  COMPLEX Ap(*), X(*), Y(*)
  !     .. Parameters ..
  COMPLEX, PARAMETER :: ZERO = (0.0E+0,0.0E+0)
  !     .. Local Scalars ..
  COMPLEX temp1, temp2
  INTEGER i, info, ix, iy, j, jx, jy, k, kk, kx, ky
  !     .. Intrinsic Functions ..
  INTRINSIC CONJG, REAL
  !* FIRST EXECUTABLE STATEMENT  CHPR2
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
    CALL XERBLA('CHPR2 ',info)
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
          temp1 = Alpha*CONJG(Y(j))
          temp2 = CONJG(Alpha*X(j))
          k = kk
          DO i = 1, j - 1
            Ap(k) = Ap(k) + X(i)*temp1 + Y(i)*temp2
            k = k + 1
          END DO
          Ap(kk+j-1) = REAL(Ap(kk+j-1)) + REAL(X(j)*temp1+Y(j)*temp2)
        ELSE
          Ap(kk+j-1) = REAL(Ap(kk+j-1))
        END IF
        kk = kk + j
      END DO
    ELSE
      DO j = 1, N
        IF ( (X(jx)/=ZERO).OR.(Y(jy)/=ZERO) ) THEN
          temp1 = Alpha*CONJG(Y(jy))
          temp2 = CONJG(Alpha*X(jx))
          ix = kx
          iy = ky
          DO k = kk, kk + j - 2
            Ap(k) = Ap(k) + X(ix)*temp1 + Y(iy)*temp2
            ix = ix + Incx
            iy = iy + Incy
          END DO
          Ap(kk+j-1) = REAL(Ap(kk+j-1)) + REAL(X(jx)*temp1+Y(jy)*temp2)
        ELSE
          Ap(kk+j-1) = REAL(Ap(kk+j-1))
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
        temp1 = Alpha*CONJG(Y(j))
        temp2 = CONJG(Alpha*X(j))
        Ap(kk) = REAL(Ap(kk)) + REAL(X(j)*temp1+Y(j)*temp2)
        k = kk + 1
        DO i = j + 1, N
          Ap(k) = Ap(k) + X(i)*temp1 + Y(i)*temp2
          k = k + 1
        END DO
      ELSE
        Ap(kk) = REAL(Ap(kk))
      END IF
      kk = kk + N - j + 1
    END DO
  ELSE
    DO j = 1, N
      IF ( (X(jx)/=ZERO).OR.(Y(jy)/=ZERO) ) THEN
        temp1 = Alpha*CONJG(Y(jy))
        temp2 = CONJG(Alpha*X(jx))
        Ap(kk) = REAL(Ap(kk)) + REAL(X(jx)*temp1+Y(jy)*temp2)
        ix = jx
        iy = jy
        DO k = kk + 1, kk + N - j
          ix = ix + Incx
          iy = iy + Incy
          Ap(k) = Ap(k) + X(ix)*temp1 + Y(iy)*temp2
        END DO
      ELSE
        Ap(kk) = REAL(Ap(kk))
      END IF
      jx = jx + Incx
      jy = jy + Incy
      kk = kk + N - j + 1
    END DO
  END IF
  !
  !
  !     End of CHPR2 .
  !
END SUBROUTINE CHPR2
