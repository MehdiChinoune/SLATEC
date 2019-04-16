!** CTRSM
SUBROUTINE CTRSM(Side,Uplo,Transa,Diag,M,N,Alpha,A,Lda,B,Ldb)
  !>
  !***
  !  Solve a complex triangular system of equations with
  !            multiple right-hand sides.
  !***
  ! **Library:**   SLATEC (BLAS)
  !***
  ! **Category:**  D1B6
  !***
  ! **Type:**      COMPLEX (STRSM-S, DTRSM-D, CTRSM-C)
  !***
  ! **Keywords:**  LEVEL 3 BLAS, LINEAR ALGEBRA
  !***
  ! **Author:**  Dongarra, J., (ANL)
  !           Duff, I., (AERE)
  !           Du Croz, J., (NAG)
  !           Hammarling, S. (NAG)
  !***
  ! **Description:**
  !
  !  CTRSM  solves one of the matrix equations
  !
  !     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
  !
  !  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
  !  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
  !
  !     op( A ) = A   or   op( A ) = A'   or   op( A ) = conjg( A' ).
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
  !              TRANSA = 'C' or 'c'   op( A ) = conjg( A' ).
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
  !  ALPHA  - COMPLEX         .
  !           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
  !           zero then  A is not referenced and  B need not be set before
  !           entry.
  !           Unchanged on exit.
  !
  !  A      - COMPLEX          array of DIMENSION ( LDA, k ), where k is m
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
  !  B      - COMPLEX          array of DIMENSION ( LDB, n ).
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
  !***
  ! **References:**  Dongarra, J., Du Croz, J., Duff, I., and Hammarling, S.
  !                 A set of level 3 basic linear algebra subprograms.
  !                 ACM TOMS, Vol. 16, No. 1, pp. 1-17, March 1990.
  !***
  ! **Routines called:**  LSAME, XERBLA

  !* REVISION HISTORY  (YYMMDD)
  !   890208  DATE WRITTEN
  !   910605  Modified to meet SLATEC prologue standards.  Only comment
  !           lines were modified.  (BKS)

  !     .. Scalar Arguments ..
  CHARACTER :: Side, Uplo, Transa, Diag
  INTEGER M, N, Lda, Ldb
  COMPLEX Alpha
  !     .. Array Arguments ..
  COMPLEX A(Lda,*), B(Ldb,*)
  !     .. Intrinsic Functions ..
  INTRINSIC CONJG, MAX
  !     .. Local Scalars ..
  LOGICAL lside, noconj, nounit, upper
  INTEGER i, info, j, k, nrowa
  COMPLEX temp
  !     .. Parameters ..
  COMPLEX, PARAMETER :: ONE = (1.0E+0,0.0E+0)
  COMPLEX, PARAMETER :: ZERO = (0.0E+0,0.0E+0)
  !* FIRST EXECUTABLE STATEMENT  CTRSM
  !
  !     Test the input parameters.
  !
  lside = LSAME(Side,'L')
  IF ( lside ) THEN
    nrowa = M
  ELSE
    nrowa = N
  END IF
  noconj = LSAME(Transa,'T')
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
  END IF
  IF ( info/=0 ) THEN
    CALL XERBLA('CTRSM ',info)
    RETURN
  END IF
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
      END DO
    END DO
    RETURN
  END IF
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
            END DO
          END IF
          DO k = M, 1, -1
            IF ( B(k,j)/=ZERO ) THEN
              IF ( nounit ) B(k,j) = B(k,j)/A(k,k)
              DO i = 1, k - 1
                B(i,j) = B(i,j) - B(k,j)*A(i,k)
              END DO
            END IF
          END DO
        END DO
      ELSE
        DO j = 1, N
          IF ( Alpha/=ONE ) THEN
            DO i = 1, M
              B(i,j) = Alpha*B(i,j)
            END DO
          END IF
          DO k = 1, M
            IF ( B(k,j)/=ZERO ) THEN
              IF ( nounit ) B(k,j) = B(k,j)/A(k,k)
              DO i = k + 1, M
                B(i,j) = B(i,j) - B(k,j)*A(i,k)
              END DO
            END IF
          END DO
        END DO
      END IF
      !
      !           Form  B := alpha*inv( A' )*B
      !           or    B := alpha*inv( conjg( A' ) )*B.
      !
    ELSEIF ( upper ) THEN
      DO j = 1, N
        DO i = 1, M
          temp = Alpha*B(i,j)
          IF ( noconj ) THEN
            DO k = 1, i - 1
              temp = temp - A(k,i)*B(k,j)
            END DO
            IF ( nounit ) temp = temp/A(i,i)
          ELSE
            DO k = 1, i - 1
              temp = temp - CONJG(A(k,i))*B(k,j)
            END DO
            IF ( nounit ) temp = temp/CONJG(A(i,i))
          END IF
          B(i,j) = temp
        END DO
      END DO
    ELSE
      DO j = 1, N
        DO i = M, 1, -1
          temp = Alpha*B(i,j)
          IF ( noconj ) THEN
            DO k = i + 1, M
              temp = temp - A(k,i)*B(k,j)
            END DO
            IF ( nounit ) temp = temp/A(i,i)
          ELSE
            DO k = i + 1, M
              temp = temp - CONJG(A(k,i))*B(k,j)
            END DO
            IF ( nounit ) temp = temp/CONJG(A(i,i))
          END IF
          B(i,j) = temp
        END DO
      END DO
    END IF
  ELSEIF ( LSAME(Transa,'N') ) THEN
    !
    !           Form  B := alpha*B*inv( A ).
    !
    IF ( upper ) THEN
      DO j = 1, N
        IF ( Alpha/=ONE ) THEN
          DO i = 1, M
            B(i,j) = Alpha*B(i,j)
          END DO
        END IF
        DO k = 1, j - 1
          IF ( A(k,j)/=ZERO ) THEN
            DO i = 1, M
              B(i,j) = B(i,j) - A(k,j)*B(i,k)
            END DO
          END IF
        END DO
        IF ( nounit ) THEN
          temp = ONE/A(j,j)
          DO i = 1, M
            B(i,j) = temp*B(i,j)
          END DO
        END IF
      END DO
    ELSE
      DO j = N, 1, -1
        IF ( Alpha/=ONE ) THEN
          DO i = 1, M
            B(i,j) = Alpha*B(i,j)
          END DO
        END IF
        DO k = j + 1, N
          IF ( A(k,j)/=ZERO ) THEN
            DO i = 1, M
              B(i,j) = B(i,j) - A(k,j)*B(i,k)
            END DO
          END IF
        END DO
        IF ( nounit ) THEN
          temp = ONE/A(j,j)
          DO i = 1, M
            B(i,j) = temp*B(i,j)
          END DO
        END IF
      END DO
    END IF
    !
    !           Form  B := alpha*B*inv( A' )
    !           or    B := alpha*B*inv( conjg( A' ) ).
    !
  ELSEIF ( upper ) THEN
    DO k = N, 1, -1
      IF ( nounit ) THEN
        IF ( noconj ) THEN
          temp = ONE/A(k,k)
        ELSE
          temp = ONE/CONJG(A(k,k))
        END IF
        DO i = 1, M
          B(i,k) = temp*B(i,k)
        END DO
      END IF
      DO j = 1, k - 1
        IF ( A(j,k)/=ZERO ) THEN
          IF ( noconj ) THEN
            temp = A(j,k)
          ELSE
            temp = CONJG(A(j,k))
          END IF
          DO i = 1, M
            B(i,j) = B(i,j) - temp*B(i,k)
          END DO
        END IF
      END DO
      IF ( Alpha/=ONE ) THEN
        DO i = 1, M
          B(i,k) = Alpha*B(i,k)
        END DO
      END IF
    END DO
  ELSE
    DO k = 1, N
      IF ( nounit ) THEN
        IF ( noconj ) THEN
          temp = ONE/A(k,k)
        ELSE
          temp = ONE/CONJG(A(k,k))
        END IF
        DO i = 1, M
          B(i,k) = temp*B(i,k)
        END DO
      END IF
      DO j = k + 1, N
        IF ( A(j,k)/=ZERO ) THEN
          IF ( noconj ) THEN
            temp = A(j,k)
          ELSE
            temp = CONJG(A(j,k))
          END IF
          DO i = 1, M
            B(i,j) = B(i,j) - temp*B(i,k)
          END DO
        END IF
      END DO
      IF ( Alpha/=ONE ) THEN
        DO i = 1, M
          B(i,k) = Alpha*B(i,k)
        END DO
      END IF
    END DO
  END IF
  !
  !
  !     End of CTRSM .
  !
END SUBROUTINE CTRSM
