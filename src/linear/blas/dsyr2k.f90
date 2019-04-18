!** DSYR2K
SUBROUTINE DSYR2K(Uplo,Trans,N,K,Alpha,A,Lda,B,Ldb,Beta,C,Ldc)
  !>
  !  Perform one of the symmetric rank 2k operations.
  !***
  ! **Library:**   SLATEC (BLAS)
  !***
  ! **Category:**  D1B6
  !***
  ! **Type:**      DOUBLE PRECISION (SSYR2-S, DSYR2-D, CSYR2-C, DSYR2K-D)
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
  !  DSYR2K  performs one of the symmetric rank 2k operations
  !
  !     C := alpha*A*B' + alpha*B*A' + beta*C,
  !
  !  or
  !
  !     C := alpha*A'*B + alpha*B'*A + beta*C,
  !
  !  where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
  !  and  A and B  are  n by k  matrices  in the  first  case  and  k by n
  !  matrices in the second case.
  !
  !  Parameters
  !  ==========
  !
  !  UPLO   - CHARACTER*1.
  !           On  entry,   UPLO  specifies  whether  the  upper  or  lower
  !           triangular  part  of the  array  C  is to be  referenced  as
  !           follows:
  !
  !              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
  !                                  is to be referenced.
  !
  !              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
  !                                  is to be referenced.
  !
  !           Unchanged on exit.
  !
  !  TRANS  - CHARACTER*1.
  !           On entry,  TRANS  specifies the operation to be performed as
  !           follows:
  !
  !              TRANS = 'N' or 'n'   C := alpha*A*B' + alpha*B*A' +
  !                                        beta*C.
  !
  !              TRANS = 'T' or 't'   C := alpha*A'*B + alpha*B'*A +
  !                                        beta*C.
  !
  !              TRANS = 'C' or 'c'   C := alpha*A'*B + alpha*B'*A +
  !                                        beta*C.
  !
  !           Unchanged on exit.
  !
  !  N      - INTEGER.
  !           On entry,  N specifies the order of the matrix C.  N must be
  !           at least zero.
  !           Unchanged on exit.
  !
  !  K      - INTEGER.
  !           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
  !           of  columns  of the  matrices  A and B,  and on  entry  with
  !           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number
  !           of rows of the matrices  A and B.  K must be at least  zero.
  !           Unchanged on exit.
  !
  !  ALPHA  - DOUBLE PRECISION.
  !           On entry, ALPHA specifies the scalar alpha.
  !           Unchanged on exit.
  !
  !  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
  !           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
  !           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
  !           part of the array  A  must contain the matrix  A,  otherwise
  !           the leading  k by n  part of the array  A  must contain  the
  !           matrix A.
  !           Unchanged on exit.
  !
  !  LDA    - INTEGER.
  !           On entry, LDA specifies the first dimension of A as declared
  !           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
  !           then  LDA must be at least  max( 1, n ), otherwise  LDA must
  !           be at least  max( 1, k ).
  !           Unchanged on exit.
  !
  !  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
  !           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
  !           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
  !           part of the array  B  must contain the matrix  B,  otherwise
  !           the leading  k by n  part of the array  B  must contain  the
  !           matrix B.
  !           Unchanged on exit.
  !
  !  LDB    - INTEGER.
  !           On entry, LDB specifies the first dimension of B as declared
  !           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
  !           then  LDB must be at least  max( 1, n ), otherwise  LDB must
  !           be at least  max( 1, k ).
  !           Unchanged on exit.
  !
  !  BETA   - DOUBLE PRECISION.
  !           On entry, BETA specifies the scalar beta.
  !           Unchanged on exit.
  !
  !  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
  !           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
  !           upper triangular part of the array C must contain the upper
  !           triangular part  of the  symmetric matrix  and the strictly
  !           lower triangular part of C is not referenced.  On exit, the
  !           upper triangular part of the array  C is overwritten by the
  !           upper triangular part of the updated matrix.
  !           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
  !           lower triangular part of the array C must contain the lower
  !           triangular part  of the  symmetric matrix  and the strictly
  !           upper triangular part of C is not referenced.  On exit, the
  !           lower triangular part of the array  C is overwritten by the
  !           lower triangular part of the updated matrix.
  !
  !  LDC    - INTEGER.
  !           On entry, LDC specifies the first dimension of C as declared
  !           in  the  calling  (sub)  program.   LDC  must  be  at  least
  !           max( 1, n ).
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
  USE service, ONLY : XERBLA
  !     .. Scalar Arguments ..
  CHARACTER :: Uplo, Trans
  INTEGER N, K, Lda, Ldb, Ldc
  REAL(8) :: Alpha, Beta
  !     .. Array Arguments ..
  REAL(8) :: A(Lda,*), B(Ldb,*), C(Ldc,*)
  !
  !     .. Intrinsic Functions ..
  INTRINSIC MAX
  !     .. Local Scalars ..
  LOGICAL upper
  INTEGER i, info, j, l, nrowa
  REAL(8) :: temp1, temp2
  !     .. Parameters ..
  REAL(8), PARAMETER :: ONE = 1.0D+0, ZERO = 0.0D+0
  !* FIRST EXECUTABLE STATEMENT  DSYR2K
  !
  !     Test the input parameters.
  !
  IF ( LSAME(Trans,'N') ) THEN
    nrowa = N
  ELSE
    nrowa = K
  END IF
  upper = LSAME(Uplo,'U')
  !
  info = 0
  IF ( (.NOT.upper).AND.(.NOT.LSAME(Uplo,'L')) ) THEN
    info = 1
  ELSEIF ( (.NOT.LSAME(Trans,'N')).AND.(.NOT.LSAME(Trans,'T')).AND.&
      (.NOT.LSAME(Trans,'C')) ) THEN
    info = 2
  ELSEIF ( N<0 ) THEN
    info = 3
  ELSEIF ( K<0 ) THEN
    info = 4
  ELSEIF ( Lda<MAX(1,nrowa) ) THEN
    info = 7
  ELSEIF ( Ldb<MAX(1,nrowa) ) THEN
    info = 9
  ELSEIF ( Ldc<MAX(1,N) ) THEN
    info = 12
  END IF
  IF ( info/=0 ) THEN
    CALL XERBLA('DSYR2K',info)
    RETURN
  END IF
  !
  !     Quick return if possible.
  !
  IF ( (N==0).OR.(((Alpha==ZERO).OR.(K==0)).AND.(Beta==ONE)) ) RETURN
  !
  !     And when  alpha.eq.zero.
  !
  IF ( Alpha==ZERO ) THEN
    IF ( upper ) THEN
      IF ( Beta==ZERO ) THEN
        DO j = 1, N
          DO i = 1, j
            C(i,j) = ZERO
          END DO
        END DO
      ELSE
        DO j = 1, N
          DO i = 1, j
            C(i,j) = Beta*C(i,j)
          END DO
        END DO
      END IF
    ELSEIF ( Beta==ZERO ) THEN
      DO j = 1, N
        DO i = j, N
          C(i,j) = ZERO
        END DO
      END DO
    ELSE
      DO j = 1, N
        DO i = j, N
          C(i,j) = Beta*C(i,j)
        END DO
      END DO
    END IF
    RETURN
  END IF
  !
  !     Start the operations.
  !
  IF ( LSAME(Trans,'N') ) THEN
    !
    !        Form  C := alpha*A*B' + alpha*B*A' + C.
    !
    IF ( upper ) THEN
      DO j = 1, N
        IF ( Beta==ZERO ) THEN
          DO i = 1, j
            C(i,j) = ZERO
          END DO
        ELSEIF ( Beta/=ONE ) THEN
          DO i = 1, j
            C(i,j) = Beta*C(i,j)
          END DO
        END IF
        DO l = 1, K
          IF ( (A(j,l)/=ZERO).OR.(B(j,l)/=ZERO) ) THEN
            temp1 = Alpha*B(j,l)
            temp2 = Alpha*A(j,l)
            DO i = 1, j
              C(i,j) = C(i,j) + A(i,l)*temp1 + B(i,l)*temp2
            END DO
          END IF
        END DO
      END DO
    ELSE
      DO j = 1, N
        IF ( Beta==ZERO ) THEN
          DO i = j, N
            C(i,j) = ZERO
          END DO
        ELSEIF ( Beta/=ONE ) THEN
          DO i = j, N
            C(i,j) = Beta*C(i,j)
          END DO
        END IF
        DO l = 1, K
          IF ( (A(j,l)/=ZERO).OR.(B(j,l)/=ZERO) ) THEN
            temp1 = Alpha*B(j,l)
            temp2 = Alpha*A(j,l)
            DO i = j, N
              C(i,j) = C(i,j) + A(i,l)*temp1 + B(i,l)*temp2
            END DO
          END IF
        END DO
      END DO
    END IF
    !
    !        Form  C := alpha*A'*B + alpha*B'*A + C.
    !
  ELSEIF ( upper ) THEN
    DO j = 1, N
      DO i = 1, j
        temp1 = ZERO
        temp2 = ZERO
        DO l = 1, K
          temp1 = temp1 + A(l,i)*B(l,j)
          temp2 = temp2 + B(l,i)*A(l,j)
        END DO
        IF ( Beta==ZERO ) THEN
          C(i,j) = Alpha*temp1 + Alpha*temp2
        ELSE
          C(i,j) = Beta*C(i,j) + Alpha*temp1 + Alpha*temp2
        END IF
      END DO
    END DO
  ELSE
    DO j = 1, N
      DO i = j, N
        temp1 = ZERO
        temp2 = ZERO
        DO l = 1, K
          temp1 = temp1 + A(l,i)*B(l,j)
          temp2 = temp2 + B(l,i)*A(l,j)
        END DO
        IF ( Beta==ZERO ) THEN
          C(i,j) = Alpha*temp1 + Alpha*temp2
        ELSE
          C(i,j) = Beta*C(i,j) + Alpha*temp1 + Alpha*temp2
        END IF
      END DO
    END DO
  END IF
  !
  !
  !     End of DSYR2K.
  !
END SUBROUTINE DSYR2K
