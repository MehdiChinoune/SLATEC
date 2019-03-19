!** CSYMM
SUBROUTINE CSYMM(Side,Uplo,M,N,Alpha,A,Lda,B,Ldb,Beta,C,Ldc)
  IMPLICIT NONE
  !>
  !***
  !  Multiply a complex general matrix by a complex symmetric
  !            matrix.
  !***
  ! **Library:**   SLATEC (BLAS)
  !***
  ! **Category:**  D1B6
  !***
  ! **Type:**      COMPLEX (SSYMM-S, DSYMM-D, CSYMM-C)
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
  !  CSYMM  performs one of the matrix-matrix operations
  !
  !     C := alpha*A*B + beta*C,
  !
  !  or
  !
  !     C := alpha*B*A + beta*C,
  !
  !  where  alpha and beta are scalars, A is a symmetric matrix and  B and
  !  C are m by n matrices.
  !
  !  Parameters
  !  ==========
  !
  !  SIDE   - CHARACTER*1.
  !           On entry,  SIDE  specifies whether  the  symmetric matrix  A
  !           appears on the  left or right  in the  operation as follows:
  !
  !              SIDE = 'L' or 'l'   C := alpha*A*B + beta*C,
  !
  !              SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,
  !
  !           Unchanged on exit.
  !
  !  UPLO   - CHARACTER*1.
  !           On  entry,   UPLO  specifies  whether  the  upper  or  lower
  !           triangular  part  of  the  symmetric  matrix   A  is  to  be
  !           referenced as follows:
  !
  !              UPLO = 'U' or 'u'   Only the upper triangular part of the
  !                                  symmetric matrix is to be referenced.
  !
  !              UPLO = 'L' or 'l'   Only the lower triangular part of the
  !                                  symmetric matrix is to be referenced.
  !
  !           Unchanged on exit.
  !
  !  M      - INTEGER.
  !           On entry,  M  specifies the number of rows of the matrix  C.
  !           M  must be at least zero.
  !           Unchanged on exit.
  !
  !  N      - INTEGER.
  !           On entry, N specifies the number of columns of the matrix C.
  !           N  must be at least zero.
  !           Unchanged on exit.
  !
  !  ALPHA  - COMPLEX         .
  !           On entry, ALPHA specifies the scalar alpha.
  !           Unchanged on exit.
  !
  !  A      - COMPLEX          array of DIMENSION ( LDA, ka ), where ka is
  !           m  when  SIDE = 'L' or 'l'  and is n  otherwise.
  !           Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of
  !           the array  A  must contain the  symmetric matrix,  such that
  !           when  UPLO = 'U' or 'u', the leading m by m upper triangular
  !           part of the array  A  must contain the upper triangular part
  !           of the  symmetric matrix and the  strictly  lower triangular
  !           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
  !           the leading  m by m  lower triangular part  of the  array  A
  !           must  contain  the  lower triangular part  of the  symmetric
  !           matrix and the  strictly upper triangular part of  A  is not
  !           referenced.
  !           Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of
  !           the array  A  must contain the  symmetric matrix,  such that
  !           when  UPLO = 'U' or 'u', the leading n by n upper triangular
  !           part of the array  A  must contain the upper triangular part
  !           of the  symmetric matrix and the  strictly  lower triangular
  !           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
  !           the leading  n by n  lower triangular part  of the  array  A
  !           must  contain  the  lower triangular part  of the  symmetric
  !           matrix and the  strictly upper triangular part of  A  is not
  !           referenced.
  !           Unchanged on exit.
  !
  !  LDA    - INTEGER.
  !           On entry, LDA specifies the first dimension of A as declared
  !           in the  calling (sub) program. When  SIDE = 'L' or 'l'  then
  !           LDA must be at least  max( 1, m ), otherwise  LDA must be at
  !           least max( 1, n ).
  !           Unchanged on exit.
  !
  !  B      - COMPLEX          array of DIMENSION ( LDB, n ).
  !           Before entry, the leading  m by n part of the array  B  must
  !           contain the matrix B.
  !           Unchanged on exit.
  !
  !  LDB    - INTEGER.
  !           On entry, LDB specifies the first dimension of B as declared
  !           in  the  calling  (sub)  program.   LDB  must  be  at  least
  !           max( 1, m ).
  !           Unchanged on exit.
  !
  !  BETA   - COMPLEX         .
  !           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
  !           supplied as zero then C need not be set on input.
  !           Unchanged on exit.
  !
  !  C      - COMPLEX          array of DIMENSION ( LDC, n ).
  !           Before entry, the leading  m by n  part of the array  C must
  !           contain the matrix  C,  except when  beta  is zero, in which
  !           case C need not be set on entry.
  !           On exit, the array  C  is overwritten by the  m by n updated
  !           matrix.
  !
  !  LDC    - INTEGER.
  !           On entry, LDC specifies the first dimension of C as declared
  !           in  the  calling  (sub)  program.   LDC  must  be  at  least
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
  CHARACTER :: Side, Uplo
  INTEGER M, N, Lda, Ldb, Ldc
  COMPLEX Alpha, Beta
  !     .. Array Arguments ..
  COMPLEX A(Lda,*), B(Ldb,*), C(Ldc,*)
  !     .. External Functions ..
  LOGICAL LSAME
  EXTERNAL LSAME
  !     .. External Subroutines ..
  EXTERNAL XERBLA
  !     .. Intrinsic Functions ..
  INTRINSIC MAX
  !     .. Local Scalars ..
  LOGICAL upper
  INTEGER i, info, j, k, nrowa
  COMPLEX temp1, temp2
  !     .. Parameters ..
  COMPLEX ONE
  PARAMETER (ONE=(1.0E+0,0.0E+0))
  COMPLEX ZERO
  PARAMETER (ZERO=(0.0E+0,0.0E+0))
  !* FIRST EXECUTABLE STATEMENT  CSYMM
  !
  !     Set NROWA as the number of rows of A.
  !
  IF ( LSAME(Side,'L') ) THEN
    nrowa = M
  ELSE
    nrowa = N
  ENDIF
  upper = LSAME(Uplo,'U')
  !
  !     Test the input parameters.
  !
  info = 0
  IF ( (.NOT.LSAME(Side,'L')).AND.(.NOT.LSAME(Side,'R')) ) THEN
    info = 1
  ELSEIF ( (.NOT.upper).AND.(.NOT.LSAME(Uplo,'L')) ) THEN
    info = 2
  ELSEIF ( M<0 ) THEN
    info = 3
  ELSEIF ( N<0 ) THEN
    info = 4
  ELSEIF ( Lda<MAX(1,nrowa) ) THEN
    info = 7
  ELSEIF ( Ldb<MAX(1,M) ) THEN
    info = 9
  ELSEIF ( Ldc<MAX(1,M) ) THEN
    info = 12
  ENDIF
  IF ( info/=0 ) THEN
    CALL XERBLA('CSYMM ',info)
    RETURN
  ENDIF
  !
  !     Quick return if possible.
  !
  IF ( (M==0).OR.(N==0).OR.((Alpha==ZERO).AND.(Beta==ONE)) ) RETURN
  !
  !     And when  alpha.eq.zero.
  !
  IF ( Alpha==ZERO ) THEN
    IF ( Beta==ZERO ) THEN
      DO j = 1, N
        DO i = 1, M
          C(i,j) = ZERO
        ENDDO
      ENDDO
    ELSE
      DO j = 1, N
        DO i = 1, M
          C(i,j) = Beta*C(i,j)
        ENDDO
      ENDDO
    ENDIF
    RETURN
  ENDIF
  !
  !     Start the operations.
  !
  IF ( .NOT.(LSAME(Side,'L')) ) THEN
    !
    !        Form  C := alpha*B*A + beta*C.
    !
    DO j = 1, N
      temp1 = Alpha*A(j,j)
      IF ( Beta==ZERO ) THEN
        DO i = 1, M
          C(i,j) = temp1*B(i,j)
        ENDDO
      ELSE
        DO i = 1, M
          C(i,j) = Beta*C(i,j) + temp1*B(i,j)
        ENDDO
      ENDIF
      DO k = 1, j - 1
        IF ( upper ) THEN
          temp1 = Alpha*A(k,j)
        ELSE
          temp1 = Alpha*A(j,k)
        ENDIF
        DO i = 1, M
          C(i,j) = C(i,j) + temp1*B(i,k)
        ENDDO
      ENDDO
      DO k = j + 1, N
        IF ( upper ) THEN
          temp1 = Alpha*A(j,k)
        ELSE
          temp1 = Alpha*A(k,j)
        ENDIF
        DO i = 1, M
          C(i,j) = C(i,j) + temp1*B(i,k)
        ENDDO
      ENDDO
    ENDDO
    !
    !        Form  C := alpha*A*B + beta*C.
    !
  ELSEIF ( upper ) THEN
    DO j = 1, N
      DO i = 1, M
        temp1 = Alpha*B(i,j)
        temp2 = ZERO
        DO k = 1, i - 1
          C(k,j) = C(k,j) + temp1*A(k,i)
          temp2 = temp2 + B(k,j)*A(k,i)
        ENDDO
        IF ( Beta==ZERO ) THEN
          C(i,j) = temp1*A(i,i) + Alpha*temp2
        ELSE
          C(i,j) = Beta*C(i,j) + temp1*A(i,i) + Alpha*temp2
        ENDIF
      ENDDO
    ENDDO
  ELSE
    DO j = 1, N
      DO i = M, 1, -1
        temp1 = Alpha*B(i,j)
        temp2 = ZERO
        DO k = i + 1, M
          C(k,j) = C(k,j) + temp1*A(k,i)
          temp2 = temp2 + B(k,j)*A(k,i)
        ENDDO
        IF ( Beta==ZERO ) THEN
          C(i,j) = temp1*A(i,i) + Alpha*temp2
        ELSE
          C(i,j) = Beta*C(i,j) + temp1*A(i,i) + Alpha*temp2
        ENDIF
      ENDDO
    ENDDO
  ENDIF
  !
  !
  !     End of CSYMM .
  !
END SUBROUTINE CSYMM
