!*==STRSM.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK STRSM
SUBROUTINE STRSM(Side,Uplo,Transa,Diag,M,N,Alpha,A,Lda,B,Ldb)
  IMPLICIT NONE
  !*--STRSM5
  !***BEGIN PROLOGUE  STRSM
  !***PURPOSE  Solve a real triangular system of equations with multiple
  !            right-hand sides.
  !***LIBRARY   SLATEC (BLAS)
  !***CATEGORY  D1B6
  !***TYPE      SINGLE PRECISION (STRSM-S, DTRSM-D, CTRSM-C)
  !***KEYWORDS  LEVEL 3 BLAS, LINEAR ALGEBRA
  !***AUTHOR  Dongarra, J., (ANL)
  !           Duff, I., (AERE)
  !           Du Croz, J., (NAG)
  !           Hammarling, S. (NAG)
  !***DESCRIPTION
  !
  !  STRSM  solves one of the matrix equations
  !
  !     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
  !
  !  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
  !  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
  !
  !     op( A ) = A   or   op( A ) = A'.
  !
  !  The matrix X is overwritten on B.
  !
  !  Parameters
  !  ==========
  !
  !  SIDE   - CHARACTER*1.
  !           On entry, SIDE specifies whether op( A ) appears on the left
  !           or right of X as follows:
  !
  !              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
  !
  !              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
  !
  !           Unchanged on exit.
  !
  !  UPLO   - CHARACTER*1.
  !           On entry, UPLO specifies whether the matrix A is an upper or
  !           lower triangular matrix as follows:
  !
  !              UPLO = 'U' or 'u'   A is an upper triangular matrix.
  !
  !              UPLO = 'L' or 'l'   A is a lower triangular matrix.
  !
  !           Unchanged on exit.
  !
  !  TRANSA - CHARACTER*1.
  !           On entry, TRANSA specifies the form of op( A ) to be used in
  !           the matrix multiplication as follows:
  !
  !              TRANSA = 'N' or 'n'   op( A ) = A.
  !
  !              TRANSA = 'T' or 't'   op( A ) = A'.
  !
  !              TRANSA = 'C' or 'c'   op( A ) = A'.
  !
  !           Unchanged on exit.
  !
  !  DIAG   - CHARACTER*1.
  !           On entry, DIAG specifies whether or not A is unit triangular
  !           as follows:
  !
  !              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
  !
  !              DIAG = 'N' or 'n'   A is not assumed to be unit
  !                                  triangular.
  !
  !           Unchanged on exit.
  !
  !  M      - INTEGER.
  !           On entry, M specifies the number of rows of B. M must be at
  !           least zero.
  !           Unchanged on exit.
  !
  !  N      - INTEGER.
  !           On entry, N specifies the number of columns of B.  N must be
  !           at least zero.
  !           Unchanged on exit.
  !
  !  ALPHA  - REAL            .
  !           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
  !           zero then  A is not referenced and  B need not be set before
  !           entry.
  !           Unchanged on exit.
  !
  !  A      - REAL             array of DIMENSION ( LDA, k ), where k is m
  !           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
  !           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
  !           upper triangular part of the array  A must contain the upper
  !           triangular matrix  and the strictly lower triangular part of
  !           A is not referenced.
  !           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
  !           lower triangular part of the array  A must contain the lower
  !           triangular matrix  and the strictly upper triangular part of
  !           A is not referenced.
  !           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
  !           A  are not referenced either,  but are assumed to be  unity.
  !           Unchanged on exit.
  !
  !  LDA    - INTEGER.
  !           On entry, LDA specifies the first dimension of A as declared
  !           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
  !           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
  !           then LDA must be at least max( 1, n ).
  !           Unchanged on exit.
  !
  !  B      - REAL             array of DIMENSION ( LDB, n ).
  !           Before entry,  the leading  m by n part of the array  B must
  !           contain  the  right-hand  side  matrix  B,  and  on exit  is
  !           overwritten by the solution matrix  X.
  !
  !  LDB    - INTEGER.
  !           On entry, LDB specifies the first dimension of B as declared
  !           in  the  calling  (sub)  program.   LDB  must  be  at  least
  !           max( 1, m ).
  !           Unchanged on exit.
  !
  !***REFERENCES  Dongarra, J., Du Croz, J., Duff, I., and Hammarling, S.
  !                 A set of level 3 basic linear algebra subprograms.
  !                 ACM TOMS, Vol. 16, No. 1, pp. 1-17, March 1990.
  !***ROUTINES CALLED  LSAME, XERBLA
  !***REVISION HISTORY  (YYMMDD)
  !   890208  DATE WRITTEN
  !   910605  Modified to meet SLATEC prologue standards.  Only comment
  !           lines were modified.  (BKS)
  !***END PROLOGUE  STRSM
  !     .. Scalar Arguments ..
  CHARACTER :: Side, Uplo, Transa, Diag
  INTEGER M, N, Lda, Ldb
  REAL Alpha
  !     .. Array Arguments ..
  REAL A(Lda,*), B(Ldb,*)
  !
  !     .. External Functions ..
  LOGICAL LSAME
  EXTERNAL LSAME
  !     .. External Subroutines ..
  EXTERNAL XERBLA
  !     .. Intrinsic Functions ..
  INTRINSIC MAX
  !     .. Local Scalars ..
  LOGICAL lside, nounit, upper
  INTEGER i, info, j, k, nrowa
  REAL temp
  !     .. Parameters ..
  REAL ONE, ZERO
  PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
  !***FIRST EXECUTABLE STATEMENT  STRSM
  !
  !     Test the input parameters.
  !
  lside = LSAME(Side,'L')
  IF ( lside ) THEN
    nrowa = M
  ELSE
    nrowa = N
  ENDIF
  nounit = LSAME(Diag,'N')
  upper = LSAME(Uplo,'U')
  !
  info = 0
  IF ( (.NOT.lside).AND.(.NOT.LSAME(Side,'R')) ) THEN
    info = 1
  ELSEIF ( (.NOT.upper).AND.(.NOT.LSAME(Uplo,'L')) ) THEN
    info = 2
  ELSEIF ( (.NOT.LSAME(Transa,'N')).AND.(.NOT.LSAME(Transa,'T')).AND.&
      (.NOT.LSAME(Transa,'C')) ) THEN
    info = 3
  ELSEIF ( (.NOT.LSAME(Diag,'U')).AND.(.NOT.LSAME(Diag,'N')) ) THEN
    info = 4
  ELSEIF ( M<0 ) THEN
    info = 5
  ELSEIF ( N<0 ) THEN
    info = 6
  ELSEIF ( Lda<MAX(1,nrowa) ) THEN
    info = 9
  ELSEIF ( Ldb<MAX(1,M) ) THEN
    info = 11
  ENDIF
  IF ( info/=0 ) THEN
    CALL XERBLA('STRSM ',info)
    RETURN
  ENDIF
  !
  !     Quick return if possible.
  !
  IF ( N==0 ) RETURN
  !
  !     And when  alpha.eq.zero.
  !
  IF ( Alpha==ZERO ) THEN
    DO j = 1, N
      DO i = 1, M
        B(i,j) = ZERO
      ENDDO
    ENDDO
    RETURN
  ENDIF
  !
  !     Start the operations.
  !
  IF ( lside ) THEN
    IF ( LSAME(Transa,'N') ) THEN
      !
      !           Form  B := alpha*inv( A )*B.
      !
      IF ( upper ) THEN
        DO j = 1, N
          IF ( Alpha/=ONE ) THEN
            DO i = 1, M
              B(i,j) = Alpha*B(i,j)
            ENDDO
          ENDIF
          DO k = M, 1, -1
            IF ( B(k,j)/=ZERO ) THEN
              IF ( nounit ) B(k,j) = B(k,j)/A(k,k)
              DO i = 1, k - 1
                B(i,j) = B(i,j) - B(k,j)*A(i,k)
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ELSE
        DO j = 1, N
          IF ( Alpha/=ONE ) THEN
            DO i = 1, M
              B(i,j) = Alpha*B(i,j)
            ENDDO
          ENDIF
          DO k = 1, M
            IF ( B(k,j)/=ZERO ) THEN
              IF ( nounit ) B(k,j) = B(k,j)/A(k,k)
              DO i = k + 1, M
                B(i,j) = B(i,j) - B(k,j)*A(i,k)
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDIF
      !
      !           Form  B := alpha*inv( A' )*B.
      !
    ELSEIF ( upper ) THEN
      DO j = 1, N
        DO i = 1, M
          temp = Alpha*B(i,j)
          DO k = 1, i - 1
            temp = temp - A(k,i)*B(k,j)
          ENDDO
          IF ( nounit ) temp = temp/A(i,i)
          B(i,j) = temp
        ENDDO
      ENDDO
    ELSE
      DO j = 1, N
        DO i = M, 1, -1
          temp = Alpha*B(i,j)
          DO k = i + 1, M
            temp = temp - A(k,i)*B(k,j)
          ENDDO
          IF ( nounit ) temp = temp/A(i,i)
          B(i,j) = temp
        ENDDO
      ENDDO
    ENDIF
  ELSEIF ( LSAME(Transa,'N') ) THEN
    !
    !           Form  B := alpha*B*inv( A ).
    !
    IF ( upper ) THEN
      DO j = 1, N
        IF ( Alpha/=ONE ) THEN
          DO i = 1, M
            B(i,j) = Alpha*B(i,j)
          ENDDO
        ENDIF
        DO k = 1, j - 1
          IF ( A(k,j)/=ZERO ) THEN
            DO i = 1, M
              B(i,j) = B(i,j) - A(k,j)*B(i,k)
            ENDDO
          ENDIF
        ENDDO
        IF ( nounit ) THEN
          temp = ONE/A(j,j)
          DO i = 1, M
            B(i,j) = temp*B(i,j)
          ENDDO
        ENDIF
      ENDDO
    ELSE
      DO j = N, 1, -1
        IF ( Alpha/=ONE ) THEN
          DO i = 1, M
            B(i,j) = Alpha*B(i,j)
          ENDDO
        ENDIF
        DO k = j + 1, N
          IF ( A(k,j)/=ZERO ) THEN
            DO i = 1, M
              B(i,j) = B(i,j) - A(k,j)*B(i,k)
            ENDDO
          ENDIF
        ENDDO
        IF ( nounit ) THEN
          temp = ONE/A(j,j)
          DO i = 1, M
            B(i,j) = temp*B(i,j)
          ENDDO
        ENDIF
      ENDDO
    ENDIF
    !
    !           Form  B := alpha*B*inv( A' ).
    !
  ELSEIF ( upper ) THEN
    DO k = N, 1, -1
      IF ( nounit ) THEN
        temp = ONE/A(k,k)
        DO i = 1, M
          B(i,k) = temp*B(i,k)
        ENDDO
      ENDIF
      DO j = 1, k - 1
        IF ( A(j,k)/=ZERO ) THEN
          temp = A(j,k)
          DO i = 1, M
            B(i,j) = B(i,j) - temp*B(i,k)
          ENDDO
        ENDIF
      ENDDO
      IF ( Alpha/=ONE ) THEN
        DO i = 1, M
          B(i,k) = Alpha*B(i,k)
        ENDDO
      ENDIF
    ENDDO
  ELSE
    DO k = 1, N
      IF ( nounit ) THEN
        temp = ONE/A(k,k)
        DO i = 1, M
          B(i,k) = temp*B(i,k)
        ENDDO
      ENDIF
      DO j = k + 1, N
        IF ( A(j,k)/=ZERO ) THEN
          temp = A(j,k)
          DO i = 1, M
            B(i,j) = B(i,j) - temp*B(i,k)
          ENDDO
        ENDIF
      ENDDO
      IF ( Alpha/=ONE ) THEN
        DO i = 1, M
          B(i,k) = Alpha*B(i,k)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
  !
  !
  !     End of STRSM .
  !
END SUBROUTINE STRSM
