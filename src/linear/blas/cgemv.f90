!** CGEMV
SUBROUTINE CGEMV(Trans,M,N,Alpha,A,Lda,X,Incx,Beta,Y,Incy)
  !>
  !  Multiply a complex vector by a complex general matrix.
  !***
  ! **Library:**   SLATEC (BLAS)
  !***
  ! **Category:**  D1B4
  !***
  ! **Type:**      COMPLEX (SGEMV-S, DGEMV-D, CGEMV-C)
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
  !  CGEMV  performs one of the matrix-vector operations
  !
  !     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,   or
  !
  !     y := alpha*conjg( A' )*x + beta*y,
  !
  !  where alpha and beta are scalars, x and y are vectors and A is an
  !  m by n matrix.
  !
  !  Parameters
  !  ==========
  !
  !  TRANS  - CHARACTER*1.
  !           On entry, TRANS specifies the operation to be performed as
  !           follows:
  !
  !              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
  !
  !              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
  !
  !              TRANS = 'C' or 'c'   y := alpha*conjg( A' )*x + beta*y.
  !
  !           Unchanged on exit.
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
  !  ALPHA  - COMPLEX         .
  !           On entry, ALPHA specifies the scalar alpha.
  !           Unchanged on exit.
  !
  !  A      - COMPLEX          array of DIMENSION ( LDA, n ).
  !           Before entry, the leading m by n part of the array A must
  !           contain the matrix of coefficients.
  !           Unchanged on exit.
  !
  !  LDA    - INTEGER.
  !           On entry, LDA specifies the first dimension of A as declared
  !           in the calling (sub) program. LDA must be at least
  !           max( 1, m ).
  !           Unchanged on exit.
  !
  !  X      - COMPLEX          array of DIMENSION at least
  !           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
  !           and at least
  !           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
  !           Before entry, the incremented array X must contain the
  !           vector x.
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
  !  Y      - COMPLEX          array of DIMENSION at least
  !           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
  !           and at least
  !           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
  !           Before entry with BETA non-zero, the incremented array Y
  !           must contain the vector y. On exit, Y is overwritten by the
  !           updated vector y.
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
  COMPLEX(SP) Alpha, Beta
  INTEGER Incx, Incy, Lda, M, N
  CHARACTER :: Trans
  !     .. Array Arguments ..
  COMPLEX(SP) A(Lda,*), X(*), Y(*)
  !     .. Parameters ..
  COMPLEX(SP), PARAMETER :: ONE = (1.0E+0,0.0E+0)
  COMPLEX(SP), PARAMETER :: ZERO = (0.0E+0,0.0E+0)
  !     .. Local Scalars ..
  COMPLEX(SP) temp
  INTEGER i, info, ix, iy, j, jx, jy, kx, ky, lenx, leny
  LOGICAL noconj
  !     .. Intrinsic Functions ..
  INTRINSIC CONJG, MAX
  !* FIRST EXECUTABLE STATEMENT  CGEMV
  !
  !     Test the input parameters.
  !
  info = 0
  IF ( .NOT.LSAME(Trans,'N').AND..NOT.LSAME(Trans,'T').AND.&
      .NOT.LSAME(Trans,'C') ) THEN
    info = 1
  ELSEIF ( M<0 ) THEN
    info = 2
  ELSEIF ( N<0 ) THEN
    info = 3
  ELSEIF ( Lda<MAX(1,M) ) THEN
    info = 6
  ELSEIF ( Incx==0 ) THEN
    info = 8
  ELSEIF ( Incy==0 ) THEN
    info = 11
  END IF
  IF ( info/=0 ) THEN
    CALL XERBLA('CGEMV ',info)
    RETURN
  END IF
  !
  !     Quick return if possible.
  !
  IF ( (M==0).OR.(N==0).OR.((Alpha==ZERO).AND.(Beta==ONE)) ) RETURN
  !
  noconj = LSAME(Trans,'T')
  !
  !     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
  !     up the start points in  X  and  Y.
  !
  IF ( LSAME(Trans,'N') ) THEN
    lenx = N
    leny = M
  ELSE
    lenx = M
    leny = N
  END IF
  IF ( Incx>0 ) THEN
    kx = 1
  ELSE
    kx = 1 - (lenx-1)*Incx
  END IF
  IF ( Incy>0 ) THEN
    ky = 1
  ELSE
    ky = 1 - (leny-1)*Incy
  END IF
  !
  !     Start the operations. In this version the elements of A are
  !     accessed sequentially with one pass through A.
  !
  !     First form  y := beta*y.
  !
  IF ( Beta/=ONE ) THEN
    IF ( Incy/=1 ) THEN
      iy = ky
      IF ( Beta==ZERO ) THEN
        DO i = 1, leny
          Y(iy) = ZERO
          iy = iy + Incy
        END DO
      ELSE
        DO i = 1, leny
          Y(iy) = Beta*Y(iy)
          iy = iy + Incy
        END DO
      END IF
    ELSEIF ( Beta==ZERO ) THEN
      DO i = 1, leny
        Y(i) = ZERO
      END DO
    ELSE
      DO i = 1, leny
        Y(i) = Beta*Y(i)
      END DO
    END IF
  END IF
  IF ( Alpha==ZERO ) RETURN
  IF ( LSAME(Trans,'N') ) THEN
    !
    !        Form  y := alpha*A*x + y.
    !
    jx = kx
    IF ( Incy==1 ) THEN
      DO j = 1, N
        IF ( X(jx)/=ZERO ) THEN
          temp = Alpha*X(jx)
          DO i = 1, M
            Y(i) = Y(i) + temp*A(i,j)
          END DO
        END IF
        jx = jx + Incx
      END DO
    ELSE
      DO j = 1, N
        IF ( X(jx)/=ZERO ) THEN
          temp = Alpha*X(jx)
          iy = ky
          DO i = 1, M
            Y(iy) = Y(iy) + temp*A(i,j)
            iy = iy + Incy
          END DO
        END IF
        jx = jx + Incx
      END DO
    END IF
  ELSE
    !
    !        Form  y := alpha*A'*x + y  or  y := alpha*conjg( A' )*x + y.
    !
    jy = ky
    IF ( Incx==1 ) THEN
      DO j = 1, N
        temp = ZERO
        IF ( noconj ) THEN
          DO i = 1, M
            temp = temp + A(i,j)*X(i)
          END DO
        ELSE
          DO i = 1, M
            temp = temp + CONJG(A(i,j))*X(i)
          END DO
        END IF
        Y(jy) = Y(jy) + Alpha*temp
        jy = jy + Incy
      END DO
    ELSE
      DO j = 1, N
        temp = ZERO
        ix = kx
        IF ( noconj ) THEN
          DO i = 1, M
            temp = temp + A(i,j)*X(ix)
            ix = ix + Incx
          END DO
        ELSE
          DO i = 1, M
            temp = temp + CONJG(A(i,j))*X(ix)
            ix = ix + Incx
          END DO
        END IF
        Y(jy) = Y(jy) + Alpha*temp
        jy = jy + Incy
      END DO
    END IF
  END IF
  !
  !
  !     End of CGEMV .
  !
END SUBROUTINE CGEMV
