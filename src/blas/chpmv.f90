!** CHPMV
SUBROUTINE CHPMV(Uplo,N,Alpha,Ap,X,Incx,Beta,Y,Incy)
  IMPLICIT NONE
  !>
  !***
  !  Perform the matrix-vector operation.
  !***
  ! **Library:**   SLATEC (BLAS)
  !***
  ! **Category:**  D1B4
  !***
  ! **Type:**      COMPLEX (SHPMV-S, DHPMV-D, CHPMV-C)
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
  !  CHPMV  performs the matrix-vector operation
  !
  !     y := alpha*A*x + beta*y,
  !
  !  where alpha and beta are scalars, x and y are n element vectors and
  !  A is an n by n hermitian matrix, supplied in packed form.
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
  !  AP     - COMPLEX          array of DIMENSION at least
  !           ( ( n*( n + 1))/2).
  !           Before entry with UPLO = 'U' or 'u', the array AP must
  !           contain the upper triangular part of the hermitian matrix
  !           packed sequentially, column by column, so that AP( 1 )
  !           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
  !           and a( 2, 2 ) respectively, and so on.
  !           Before entry with UPLO = 'L' or 'l', the array AP must
  !           contain the lower triangular part of the hermitian matrix
  !           packed sequentially, column by column, so that AP( 1 )
  !           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
  !           and a( 3, 1 ) respectively, and so on.
  !           Note that the imaginary parts of the diagonal elements need
  !           not be set and are assumed to be zero.
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
  !  BETA   - COMPLEX         .
  !           On entry, BETA specifies the scalar beta. When BETA is
  !           supplied as zero then Y need not be set on input.
  !           Unchanged on exit.
  !
  !  Y      - COMPLEX          array of dimension at least
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

  !     .. Scalar Arguments ..
  COMPLEX Alpha, Beta
  INTEGER Incx, Incy, N
  CHARACTER :: Uplo
  !     .. Array Arguments ..
  COMPLEX Ap(*), X(*), Y(*)
  !     .. Parameters ..
  COMPLEX, PARAMETER :: ONE = (1.0E+0,0.0E+0)
  COMPLEX, PARAMETER :: ZERO = (0.0E+0,0.0E+0)
  !     .. Local Scalars ..
  COMPLEX temp1, temp2
  INTEGER i, info, ix, iy, j, jx, jy, k, kk, kx, ky
  !     .. External Functions ..
  LOGICAL LSAME
  EXTERNAL LSAME
  !     .. External Subroutines ..
  EXTERNAL XERBLA
  !     .. Intrinsic Functions ..
  INTRINSIC CONJG, REAL
  !* FIRST EXECUTABLE STATEMENT  CHPMV
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
  ENDIF
  IF ( info/=0 ) THEN
    CALL XERBLA('CHPMV ',info)
    RETURN
  ENDIF
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
  ENDIF
  IF ( Incy>0 ) THEN
    ky = 1
  ELSE
    ky = 1 - (N-1)*Incy
  ENDIF
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
        ENDDO
      ELSE
        DO i = 1, N
          Y(iy) = Beta*Y(iy)
          iy = iy + Incy
        ENDDO
      ENDIF
    ELSEIF ( Beta==ZERO ) THEN
      DO i = 1, N
        Y(i) = ZERO
      ENDDO
    ELSE
      DO i = 1, N
        Y(i) = Beta*Y(i)
      ENDDO
    ENDIF
  ENDIF
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
          temp2 = temp2 + CONJG(Ap(k))*X(i)
          k = k + 1
        ENDDO
        Y(j) = Y(j) + temp1*REAL(Ap(kk+j-1)) + Alpha*temp2
        kk = kk + j
      ENDDO
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
          temp2 = temp2 + CONJG(Ap(k))*X(ix)
          ix = ix + Incx
          iy = iy + Incy
        ENDDO
        Y(jy) = Y(jy) + temp1*REAL(Ap(kk+j-1)) + Alpha*temp2
        jx = jx + Incx
        jy = jy + Incy
        kk = kk + j
      ENDDO
    ENDIF
    !
    !        Form  y  when AP contains the lower triangle.
    !
  ELSEIF ( (Incx==1).AND.(Incy==1) ) THEN
    DO j = 1, N
      temp1 = Alpha*X(j)
      temp2 = ZERO
      Y(j) = Y(j) + temp1*REAL(Ap(kk))
      k = kk + 1
      DO i = j + 1, N
        Y(i) = Y(i) + temp1*Ap(k)
        temp2 = temp2 + CONJG(Ap(k))*X(i)
        k = k + 1
      ENDDO
      Y(j) = Y(j) + Alpha*temp2
      kk = kk + (N-j+1)
    ENDDO
  ELSE
    jx = kx
    jy = ky
    DO j = 1, N
      temp1 = Alpha*X(jx)
      temp2 = ZERO
      Y(jy) = Y(jy) + temp1*REAL(Ap(kk))
      ix = jx
      iy = jy
      DO k = kk + 1, kk + N - j
        ix = ix + Incx
        iy = iy + Incy
        Y(iy) = Y(iy) + temp1*Ap(k)
        temp2 = temp2 + CONJG(Ap(k))*X(ix)
      ENDDO
      Y(jy) = Y(jy) + Alpha*temp2
      jx = jx + Incx
      jy = jy + Incy
      kk = kk + (N-j+1)
    ENDDO
  ENDIF
  !
  !
  !     End of CHPMV .
  !
END SUBROUTINE CHPMV
