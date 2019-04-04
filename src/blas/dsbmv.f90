!** DSBMV
SUBROUTINE DSBMV(Uplo,N,K,Alpha,A,Lda,X,Incx,Beta,Y,Incy)
  IMPLICIT NONE
  !>
  !***
  !  Perform the matrix-vector operation.
  !***
  ! **Library:**   SLATEC (BLAS)
  !***
  ! **Category:**  D1B4
  !***
  ! **Type:**      DOUBLE PRECISION (SSBMV-S, DSBMV-D, CSBMV-C)
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
  !  DSBMV  performs the matrix-vector  operation
  !
  !     y := alpha*A*x + beta*y,
  !
  !  where alpha and beta are scalars, x and y are n element vectors and
  !  A is an n by n symmetric band matrix, with k super-diagonals.
  !
  !  Parameters
  !  ==========
  !
  !  UPLO   - CHARACTER*1.
  !           On entry, UPLO specifies whether the upper or lower
  !           triangular part of the band matrix A is being supplied as
  !           follows:
  !
  !              UPLO = 'U' or 'u'   The upper triangular part of A is
  !                                  being supplied.
  !
  !              UPLO = 'L' or 'l'   The lower triangular part of A is
  !                                  being supplied.
  !
  !           Unchanged on exit.
  !
  !  N      - INTEGER.
  !           On entry, N specifies the order of the matrix A.
  !           N must be at least zero.
  !           Unchanged on exit.
  !
  !  K      - INTEGER.
  !           On entry, K specifies the number of super-diagonals of the
  !           matrix A. K must satisfy  0 .le. K.
  !           Unchanged on exit.
  !
  !  ALPHA  - DOUBLE PRECISION.
  !           On entry, ALPHA specifies the scalar alpha.
  !           Unchanged on exit.
  !
  !  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
  !           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
  !           by n part of the array A must contain the upper triangular
  !           band part of the symmetric matrix, supplied column by
  !           column, with the leading diagonal of the matrix in row
  !           ( k + 1 ) of the array, the first super-diagonal starting at
  !           position 2 in row k, and so on. The top left k by k triangle
  !           of the array A is not referenced.
  !           The following program segment will transfer the upper
  !           triangular part of a symmetric band matrix from conventional
  !           full matrix storage to band storage:
  !
  !                 DO 20, J = 1, N
  !                    M = K + 1 - J
  !                    DO 10, I = MAX( 1, J - K ), J
  !                       A( M + I, J ) = matrix( I, J )
  !              10    CONTINUE
  !              20 CONTINUE
  !
  !           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
  !           by n part of the array A must contain the lower triangular
  !           band part of the symmetric matrix, supplied column by
  !           column, with the leading diagonal of the matrix in row 1 of
  !           the array, the first sub-diagonal starting at position 1 in
  !           row 2, and so on. The bottom right k by k triangle of the
  !           array A is not referenced.
  !           The following program segment will transfer the lower
  !           triangular part of a symmetric band matrix from conventional
  !           full matrix storage to band storage:
  !
  !                 DO 20, J = 1, N
  !                    M = 1 - J
  !                    DO 10, I = J, MIN( N, J + K )
  !                       A( M + I, J ) = matrix( I, J )
  !              10    CONTINUE
  !              20 CONTINUE
  !
  !           Unchanged on exit.
  !
  !  LDA    - INTEGER.
  !           On entry, LDA specifies the first dimension of A as declared
  !           in the calling (sub) program. LDA must be at least
  !           ( k + 1 ).
  !           Unchanged on exit.
  !
  !  X      - DOUBLE PRECISION array of DIMENSION at least
  !           ( 1 + ( n - 1 )*abs( INCX ) ).
  !           Before entry, the incremented array X must contain the
  !           vector x.
  !           Unchanged on exit.
  !
  !  INCX   - INTEGER.
  !           On entry, INCX specifies the increment for the elements of
  !           X. INCX must not be zero.
  !           Unchanged on exit.
  !
  !  BETA   - DOUBLE PRECISION.
  !           On entry, BETA specifies the scalar beta.
  !           Unchanged on exit.
  !
  !  Y      - DOUBLE PRECISION array of DIMENSION at least
  !           ( 1 + ( n - 1 )*abs( INCY ) ).
  !           Before entry, the incremented array Y must contain the
  !           vector y. On exit, Y is overwritten by the updated vector y.
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
  REAL(8) :: Alpha, Beta
  INTEGER Incx, Incy, K, Lda, N
  CHARACTER :: Uplo
  !     .. Array Arguments ..
  REAL(8) :: A(Lda,*), X(*), Y(*)
  !     .. Parameters ..
  REAL(8), PARAMETER :: ONE = 1.0D+0, ZERO = 0.0D+0
  !     .. Local Scalars ..
  REAL(8) :: temp1, temp2
  INTEGER i, info, ix, iy, j, jx, jy, kplus1, kx, ky, l
  !     .. External Functions ..
  LOGICAL, EXTERNAL :: LSAME
  !     .. External Subroutines ..
  EXTERNAL :: XERBLA
  !     .. Intrinsic Functions ..
  INTRINSIC MAX, MIN
  !* FIRST EXECUTABLE STATEMENT  DSBMV
  !
  !     Test the input parameters.
  !
  info = 0
  IF ( .NOT.LSAME(Uplo,'U').AND..NOT.LSAME(Uplo,'L') ) THEN
    info = 1
  ELSEIF ( N<0 ) THEN
    info = 2
  ELSEIF ( K<0 ) THEN
    info = 3
  ELSEIF ( Lda<(K+1) ) THEN
    info = 6
  ELSEIF ( Incx==0 ) THEN
    info = 8
  ELSEIF ( Incy==0 ) THEN
    info = 11
  END IF
  IF ( info/=0 ) THEN
    CALL XERBLA('DSBMV ',info)
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
  !     Start the operations. In this version the elements of the array A
  !     are accessed sequentially with one pass through A.
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
  IF ( LSAME(Uplo,'U') ) THEN
    !
    !        Form  y  when upper triangle of A is stored.
    !
    kplus1 = K + 1
    IF ( (Incx==1).AND.(Incy==1) ) THEN
      DO j = 1, N
        temp1 = Alpha*X(j)
        temp2 = ZERO
        l = kplus1 - j
        DO i = MAX(1,j-K), j - 1
          Y(i) = Y(i) + temp1*A(l+i,j)
          temp2 = temp2 + A(l+i,j)*X(i)
        END DO
        Y(j) = Y(j) + temp1*A(kplus1,j) + Alpha*temp2
      END DO
    ELSE
      jx = kx
      jy = ky
      DO j = 1, N
        temp1 = Alpha*X(jx)
        temp2 = ZERO
        ix = kx
        iy = ky
        l = kplus1 - j
        DO i = MAX(1,j-K), j - 1
          Y(iy) = Y(iy) + temp1*A(l+i,j)
          temp2 = temp2 + A(l+i,j)*X(ix)
          ix = ix + Incx
          iy = iy + Incy
        END DO
        Y(jy) = Y(jy) + temp1*A(kplus1,j) + Alpha*temp2
        jx = jx + Incx
        jy = jy + Incy
        IF ( j>K ) THEN
          kx = kx + Incx
          ky = ky + Incy
        END IF
      END DO
    END IF
    !
    !        Form  y  when lower triangle of A is stored.
    !
  ELSEIF ( (Incx==1).AND.(Incy==1) ) THEN
    DO j = 1, N
      temp1 = Alpha*X(j)
      temp2 = ZERO
      Y(j) = Y(j) + temp1*A(1,j)
      l = 1 - j
      DO i = j + 1, MIN(N,j+K)
        Y(i) = Y(i) + temp1*A(l+i,j)
        temp2 = temp2 + A(l+i,j)*X(i)
      END DO
      Y(j) = Y(j) + Alpha*temp2
    END DO
  ELSE
    jx = kx
    jy = ky
    DO j = 1, N
      temp1 = Alpha*X(jx)
      temp2 = ZERO
      Y(jy) = Y(jy) + temp1*A(1,j)
      l = 1 - j
      ix = jx
      iy = jy
      DO i = j + 1, MIN(N,j+K)
        ix = ix + Incx
        iy = iy + Incy
        Y(iy) = Y(iy) + temp1*A(l+i,j)
        temp2 = temp2 + A(l+i,j)*X(ix)
      END DO
      Y(jy) = Y(jy) + Alpha*temp2
      jx = jx + Incx
      jy = jy + Incy
    END DO
  END IF
  !
  !
  !     End of DSBMV .
  !
END SUBROUTINE DSBMV
