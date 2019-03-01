!DECK CHER2K
SUBROUTINE CHER2K(Uplo,Trans,N,K,Alpha,A,Lda,B,Ldb,Beta,C,Ldc)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  CHER2K
  !***PURPOSE  Perform Hermitian rank 2k update of a complex.
  !***LIBRARY   SLATEC (BLAS)
  !***CATEGORY  D1B6
  !***TYPE      COMPLEX (SHER2-S, DHER2-D, CHER2-C, CHER2K-C)
  !***KEYWORDS  LEVEL 3 BLAS, LINEAR ALGEBRA
  !***AUTHOR  Dongarra, J., (ANL)
  !           Duff, I., (AERE)
  !           Du Croz, J., (NAG)
  !           Hammarling, S. (NAG)
  !***DESCRIPTION
  !
  !  CHER2K  performs one of the hermitian rank 2k operations
  !
  !     C := alpha*A*conjg( B' ) + conjg( alpha )*B*conjg( A' ) + beta*C,
  !
  !  or
  !
  !     C := alpha*conjg( A' )*B + conjg( alpha )*conjg( B' )*A + beta*C,
  !
  !  where  alpha and beta  are scalars with  beta  real,  C is an  n by n
  !  hermitian matrix and  A and B  are  n by k matrices in the first case
  !  and  k by n  matrices in the second case.
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
  !              TRANS = 'N' or 'n'    C := alpha*A*conjg( B' )          +
  !                                         conjg( alpha )*B*conjg( A' ) +
  !                                         beta*C.
  !
  !              TRANS = 'C' or 'c'    C := alpha*conjg( A' )*B          +
  !                                         conjg( alpha )*conjg( B' )*A +
  !                                         beta*C.
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
  !           TRANS = 'C' or 'c',  K  specifies  the number of rows of the
  !           matrices  A and B.  K must be at least zero.
  !           Unchanged on exit.
  !
  !  ALPHA  - COMPLEX         .
  !           On entry, ALPHA specifies the scalar alpha.
  !           Unchanged on exit.
  !
  !  A      - COMPLEX          array of DIMENSION ( LDA, ka ), where ka is
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
  !  B      - COMPLEX          array of DIMENSION ( LDB, kb ), where kb is
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
  !  BETA   - REAL            .
  !           On entry, BETA specifies the scalar beta.
  !           Unchanged on exit.
  !
  !  C      - COMPLEX          array of DIMENSION ( LDC, n ).
  !           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
  !           upper triangular part of the array C must contain the upper
  !           triangular part  of the  hermitian matrix  and the strictly
  !           lower triangular part of C is not referenced.  On exit, the
  !           upper triangular part of the array  C is overwritten by the
  !           upper triangular part of the updated matrix.
  !           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
  !           lower triangular part of the array C must contain the lower
  !           triangular part  of the  hermitian matrix  and the strictly
  !           upper triangular part of C is not referenced.  On exit, the
  !           lower triangular part of the array  C is overwritten by the
  !           lower triangular part of the updated matrix.
  !           Note that the imaginary parts of the diagonal elements need
  !           not be set,  they are assumed to be zero,  and on exit they
  !           are set to zero.
  !
  !  LDC    - INTEGER.
  !           On entry, LDC specifies the first dimension of C as declared
  !           in  the  calling  (sub)  program.   LDC  must  be  at  least
  !           max( 1, n ).
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
  !***END PROLOGUE  CHER2K
  !     .. Scalar Arguments ..
  CHARACTER :: Uplo, Trans
  INTEGER N, K, Lda, Ldb, Ldc
  REAL Beta
  COMPLEX Alpha
  !     .. Array Arguments ..
  COMPLEX A(Lda,*), B(Ldb,*), C(Ldc,*)
  !     .. External Functions ..
  LOGICAL LSAME
  EXTERNAL LSAME
  !     .. External Subroutines ..
  EXTERNAL XERBLA
  !     .. Intrinsic Functions ..
  INTRINSIC CONJG, MAX, REAL
  !     .. Local Scalars ..
  LOGICAL upper
  INTEGER i, info, j, l, nrowa
  COMPLEX temp1, temp2
  !     .. Parameters ..
  REAL ONE
  PARAMETER (ONE=1.0E+0)
  COMPLEX ZERO
  PARAMETER (ZERO=(0.0E+0,0.0E+0))
  !***FIRST EXECUTABLE STATEMENT  CHER2K
  !
  !     Test the input parameters.
  !
  IF ( LSAME(Trans,'N') ) THEN
    nrowa = N
  ELSE
    nrowa = K
  ENDIF
  upper = LSAME(Uplo,'U')
  !
  info = 0
  IF ( (.NOT.upper).AND.(.NOT.LSAME(Uplo,'L')) ) THEN
    info = 1
  ELSEIF ( (.NOT.LSAME(Trans,'N')).AND.(.NOT.LSAME(Trans,'C')) ) THEN
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
  ENDIF
  IF ( info/=0 ) THEN
    CALL XERBLA('CHER2K',info)
    RETURN
  ENDIF
  !
  !     Quick return if possible.
  !
  IF ( (N==0).OR.(((Alpha==ZERO).OR.(K==0)).AND.(Beta==ONE)) ) RETURN
  !
  !     And when  alpha.eq.zero.
  !
  IF ( Alpha==ZERO ) THEN
    IF ( upper ) THEN
      IF ( Beta==REAL(ZERO) ) THEN
        DO j = 1, N
          DO i = 1, j
            C(i,j) = ZERO
          ENDDO
        ENDDO
      ELSE
        DO j = 1, N
          DO i = 1, j - 1
            C(i,j) = Beta*C(i,j)
          ENDDO
          C(j,j) = Beta*REAL(C(j,j))
        ENDDO
      ENDIF
    ELSEIF ( Beta==REAL(ZERO) ) THEN
      DO j = 1, N
        DO i = j, N
          C(i,j) = ZERO
        ENDDO
      ENDDO
    ELSE
      DO j = 1, N
        C(j,j) = Beta*REAL(C(j,j))
        DO i = j + 1, N
          C(i,j) = Beta*C(i,j)
        ENDDO
      ENDDO
    ENDIF
    RETURN
  ENDIF
  !
  !     Start the operations.
  !
  IF ( LSAME(Trans,'N') ) THEN
    !
    !        Form  C := alpha*A*conjg( B' ) + conjg( alpha )*B*conjg( A' ) +
    !                   C.
    !
    IF ( upper ) THEN
      DO j = 1, N
        IF ( Beta==REAL(ZERO) ) THEN
          DO i = 1, j
            C(i,j) = ZERO
          ENDDO
        ELSEIF ( Beta/=ONE ) THEN
          DO i = 1, j - 1
            C(i,j) = Beta*C(i,j)
          ENDDO
          C(j,j) = Beta*REAL(C(j,j))
        ENDIF
        DO l = 1, K
          IF ( (A(j,l)/=ZERO).OR.(B(j,l)/=ZERO) ) THEN
            temp1 = Alpha*CONJG(B(j,l))
            temp2 = CONJG(Alpha*A(j,l))
            DO i = 1, j - 1
              C(i,j) = C(i,j) + A(i,l)*temp1 + B(i,l)*temp2
            ENDDO
            C(j,j) = REAL(C(j,j)) + REAL(A(j,l)*temp1+B(j,l)*temp2)
          ENDIF
        ENDDO
      ENDDO
    ELSE
      DO j = 1, N
        IF ( Beta==REAL(ZERO) ) THEN
          DO i = j, N
            C(i,j) = ZERO
          ENDDO
        ELSEIF ( Beta/=ONE ) THEN
          DO i = j + 1, N
            C(i,j) = Beta*C(i,j)
          ENDDO
          C(j,j) = Beta*REAL(C(j,j))
        ENDIF
        DO l = 1, K
          IF ( (A(j,l)/=ZERO).OR.(B(j,l)/=ZERO) ) THEN
            temp1 = Alpha*CONJG(B(j,l))
            temp2 = CONJG(Alpha*A(j,l))
            DO i = j + 1, N
              C(i,j) = C(i,j) + A(i,l)*temp1 + B(i,l)*temp2
            ENDDO
            C(j,j) = REAL(C(j,j)) + REAL(A(j,l)*temp1+B(j,l)*temp2)
          ENDIF
        ENDDO
      ENDDO
    ENDIF
    !
    !        Form  C := alpha*conjg( A' )*B + conjg( alpha )*conjg( B' )*A +
    !                   C.
    !
  ELSEIF ( upper ) THEN
    DO j = 1, N
      DO i = 1, j
        temp1 = ZERO
        temp2 = ZERO
        DO l = 1, K
          temp1 = temp1 + CONJG(A(l,i))*B(l,j)
          temp2 = temp2 + CONJG(B(l,i))*A(l,j)
        ENDDO
        IF ( i==j ) THEN
          IF ( Beta==REAL(ZERO) ) THEN
            C(j,j) = REAL(Alpha*temp1+CONJG(Alpha)*temp2)
          ELSE
            C(j,j) = Beta*REAL(C(j,j))&
              + REAL(Alpha*temp1+CONJG(Alpha)*temp2)
          ENDIF
        ELSEIF ( Beta==REAL(ZERO) ) THEN
          C(i,j) = Alpha*temp1 + CONJG(Alpha)*temp2
        ELSE
          C(i,j) = Beta*C(i,j) + Alpha*temp1 + CONJG(Alpha)*temp2
        ENDIF
      ENDDO
    ENDDO
  ELSE
    DO j = 1, N
      DO i = j, N
        temp1 = ZERO
        temp2 = ZERO
        DO l = 1, K
          temp1 = temp1 + CONJG(A(l,i))*B(l,j)
          temp2 = temp2 + CONJG(B(l,i))*A(l,j)
        ENDDO
        IF ( i==j ) THEN
          IF ( Beta==REAL(ZERO) ) THEN
            C(j,j) = REAL(Alpha*temp1+CONJG(Alpha)*temp2)
          ELSE
            C(j,j) = Beta*REAL(C(j,j))&
              + REAL(Alpha*temp1+CONJG(Alpha)*temp2)
          ENDIF
        ELSEIF ( Beta==REAL(ZERO) ) THEN
          C(i,j) = Alpha*temp1 + CONJG(Alpha)*temp2
        ELSE
          C(i,j) = Beta*C(i,j) + Alpha*temp1 + CONJG(Alpha)*temp2
        ENDIF
      ENDDO
    ENDDO
  ENDIF
  !
  !
  !     End of CHER2K.
  !
END SUBROUTINE CHER2K
