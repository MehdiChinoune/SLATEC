!*==DGBMV.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK DGBMV
SUBROUTINE DGBMV(Trans,M,N,Kl,Ku,Alpha,A,Lda,X,Incx,Beta,Y,Incy)
  IMPLICIT NONE
  !*--DGBMV5
  !***BEGIN PROLOGUE  DGBMV
  !***PURPOSE  Perform one of the matrix-vector operations.
  !***LIBRARY   SLATEC (BLAS)
  !***CATEGORY  D1B4
  !***TYPE      DOUBLE PRECISION (SGBMV-S, DGBMV-D, CGBMV-C)
  !***KEYWORDS  LEVEL 2 BLAS, LINEAR ALGEBRA
  !***AUTHOR  Dongarra, J. J., (ANL)
  !           Du Croz, J., (NAG)
  !           Hammarling, S., (NAG)
  !           Hanson, R. J., (SNLA)
  !***DESCRIPTION
  !
  !  DGBMV  performs one of the matrix-vector operations
  !
  !     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
  !
  !  where alpha and beta are scalars, x and y are vectors and A is an
  !  m by n band matrix, with kl sub-diagonals and ku super-diagonals.
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
  !              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
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
  !  KL     - INTEGER.
  !           On entry, KL specifies the number of sub-diagonals of the
  !           matrix A. KL must satisfy  0 .le. KL.
  !           Unchanged on exit.
  !
  !  KU     - INTEGER.
  !           On entry, KU specifies the number of super-diagonals of the
  !           matrix A. KU must satisfy  0 .le. KU.
  !           Unchanged on exit.
  !
  !  ALPHA  - DOUBLE PRECISION.
  !           On entry, ALPHA specifies the scalar alpha.
  !           Unchanged on exit.
  !
  !  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
  !           Before entry, the leading ( kl + ku + 1 ) by n part of the
  !           array A must contain the matrix of coefficients, supplied
  !           column by column, with the leading diagonal of the matrix in
  !           row ( ku + 1 ) of the array, the first super-diagonal
  !           starting at position 2 in row ku, the first sub-diagonal
  !           starting at position 1 in row ( ku + 2 ), and so on.
  !           Elements in the array A that do not correspond to elements
  !           in the band matrix (such as the top left ku by ku triangle)
  !           are not referenced.
  !           The following program segment will transfer a band matrix
  !           from conventional full matrix storage to band storage:
  !
  !                 DO 20, J = 1, N
  !                    K = KU + 1 - J
  !                    DO 10, I = MAX( 1, J - KU ), MIN( M, J + KL )
  !                       A( K + I, J ) = matrix( I, J )
  !              10    CONTINUE
  !              20 CONTINUE
  !
  !           Unchanged on exit.
  !
  !  LDA    - INTEGER.
  !           On entry, LDA specifies the first dimension of A as declared
  !           in the calling (sub) program. LDA must be at least
  !           ( kl + ku + 1 ).
  !           Unchanged on exit.
  !
  !  X      - DOUBLE PRECISION array of DIMENSION at least
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
  !  BETA   - DOUBLE PRECISION.
  !           On entry, BETA specifies the scalar beta. When BETA is
  !           supplied as zero then Y need not be set on input.
  !           Unchanged on exit.
  !
  !  Y      - DOUBLE PRECISION array of DIMENSION at least
  !           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
  !           and at least
  !           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
  !           Before entry, the incremented array Y must contain the
  !           vector y. On exit, Y is overwritten by the updated vector y.
  !
  !  INCY   - INTEGER.
  !           On entry, INCY specifies the increment for the elements of
  !           Y. INCY must not be zero.
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
  !***END PROLOGUE  DGBMV
  !     .. Scalar Arguments ..
  REAL(8) :: Alpha, Beta
  INTEGER Incx, Incy, Kl, Ku, Lda, M, N
  CHARACTER :: Trans
  !     .. Array Arguments ..
  REAL(8) :: A(Lda,*), X(*), Y(*)
  REAL(8) :: ONE, ZERO
  PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
  !     .. Local Scalars ..
  REAL(8) :: temp
  INTEGER i, info, ix, iy, j, jx, jy, k, kup1, kx, ky, lenx ,&
    leny
  !     .. External Functions ..
  LOGICAL LSAME
  EXTERNAL LSAME
  !     .. External Subroutines ..
  EXTERNAL XERBLA
  !     .. Intrinsic Functions ..
  INTRINSIC MAX, MIN
  !***FIRST EXECUTABLE STATEMENT  DGBMV
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
  ELSEIF ( Kl<0 ) THEN
    info = 4
  ELSEIF ( Ku<0 ) THEN
    info = 5
  ELSEIF ( Lda<(Kl+Ku+1) ) THEN
    info = 8
  ELSEIF ( Incx==0 ) THEN
    info = 10
  ELSEIF ( Incy==0 ) THEN
    info = 13
  ENDIF
  IF ( info/=0 ) THEN
    CALL XERBLA('DGBMV ',info)
    RETURN
  ENDIF
  !
  !     Quick return if possible.
  !
  IF ( (M==0).OR.(N==0).OR.((Alpha==ZERO).AND.(Beta==ONE)) ) RETURN
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
  ENDIF
  IF ( Incx>0 ) THEN
    kx = 1
  ELSE
    kx = 1 - (lenx-1)*Incx
  ENDIF
  IF ( Incy>0 ) THEN
    ky = 1
  ELSE
    ky = 1 - (leny-1)*Incy
  ENDIF
  !
  !     Start the operations. In this version the elements of A are
  !     accessed sequentially with one pass through the band part of A.
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
        ENDDO
      ELSE
        DO i = 1, leny
          Y(iy) = Beta*Y(iy)
          iy = iy + Incy
        ENDDO
      ENDIF
    ELSEIF ( Beta==ZERO ) THEN
      DO i = 1, leny
        Y(i) = ZERO
      ENDDO
    ELSE
      DO i = 1, leny
        Y(i) = Beta*Y(i)
      ENDDO
    ENDIF
  ENDIF
  IF ( Alpha==ZERO ) RETURN
  kup1 = Ku + 1
  IF ( LSAME(Trans,'N') ) THEN
    !
    !        Form  y := alpha*A*x + y.
    !
    jx = kx
    IF ( Incy==1 ) THEN
      DO j = 1, N
        IF ( X(jx)/=ZERO ) THEN
          temp = Alpha*X(jx)
          k = kup1 - j
          DO i = MAX(1,j-Ku), MIN(M,j+Kl)
            Y(i) = Y(i) + temp*A(k+i,j)
          ENDDO
        ENDIF
        jx = jx + Incx
      ENDDO
    ELSE
      DO j = 1, N
        IF ( X(jx)/=ZERO ) THEN
          temp = Alpha*X(jx)
          iy = ky
          k = kup1 - j
          DO i = MAX(1,j-Ku), MIN(M,j+Kl)
            Y(iy) = Y(iy) + temp*A(k+i,j)
            iy = iy + Incy
          ENDDO
        ENDIF
        jx = jx + Incx
        IF ( j>Ku ) ky = ky + Incy
      ENDDO
    ENDIF
  ELSE
    !
    !        Form  y := alpha*A'*x + y.
    !
    jy = ky
    IF ( Incx==1 ) THEN
      DO j = 1, N
        temp = ZERO
        k = kup1 - j
        DO i = MAX(1,j-Ku), MIN(M,j+Kl)
          temp = temp + A(k+i,j)*X(i)
        ENDDO
        Y(jy) = Y(jy) + Alpha*temp
        jy = jy + Incy
      ENDDO
    ELSE
      DO j = 1, N
        temp = ZERO
        ix = kx
        k = kup1 - j
        DO i = MAX(1,j-Ku), MIN(M,j+Kl)
          temp = temp + A(k+i,j)*X(ix)
          ix = ix + Incx
        ENDDO
        Y(jy) = Y(jy) + Alpha*temp
        jy = jy + Incy
        IF ( j>Ku ) kx = kx + Incx
      ENDDO
    ENDIF
  ENDIF
  !
  !
  !     End of DGBMV .
  !
END SUBROUTINE DGBMV
