!*==CGEMM.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK CGEMM
SUBROUTINE CGEMM(Transa,Transb,M,N,K,Alpha,A,Lda,B,Ldb,Beta,C,Ldc)
  IMPLICIT NONE
  !*--CGEMM5
  !***BEGIN PROLOGUE  CGEMM
  !***PURPOSE  Multiply a complex general matrix by a complex general
  !            matrix.
  !***LIBRARY   SLATEC (BLAS)
  !***CATEGORY  D1B6
  !***TYPE      COMPLEX (SGEMM-S, DGEMM-D, CGEMM-C)
  !***KEYWORDS  LEVEL 3 BLAS, LINEAR ALGEBRA
  !***AUTHOR  Dongarra, J., (ANL)
  !           Duff, I., (AERE)
  !           Du Croz, J., (NAG)
  !           Hammarling, S. (NAG)
  !***DESCRIPTION
  !
  !  CGEMM  performs one of the matrix-matrix operations
  !
  !     C := alpha*op( A )*op( B ) + beta*C,
  !
  !  where  op( X ) is one of
  !
  !     op( X ) = X   or   op( X ) = X'   or   op( X ) = conjg( X' ),
  !
  !  alpha and beta are scalars, and A, B and C are matrices, with op( A )
  !  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
  !
  !  Parameters
  !  ==========
  !
  !  TRANSA - CHARACTER*1.
  !           On entry, TRANSA specifies the form of op( A ) to be used in
  !           the matrix multiplication as follows:
  !
  !              TRANSA = 'N' or 'n',  op( A ) = A.
  !
  !              TRANSA = 'T' or 't',  op( A ) = A'.
  !
  !              TRANSA = 'C' or 'c',  op( A ) = conjg( A' ).
  !
  !           Unchanged on exit.
  !
  !  TRANSB - CHARACTER*1.
  !           On entry, TRANSB specifies the form of op( B ) to be used in
  !           the matrix multiplication as follows:
  !
  !              TRANSB = 'N' or 'n',  op( B ) = B.
  !
  !              TRANSB = 'T' or 't',  op( B ) = B'.
  !
  !              TRANSB = 'C' or 'c',  op( B ) = conjg( B' ).
  !
  !           Unchanged on exit.
  !
  !  M      - INTEGER.
  !           On entry,  M  specifies  the number  of rows  of the  matrix
  !           op( A )  and of the  matrix  C.  M  must  be at least  zero.
  !           Unchanged on exit.
  !
  !  N      - INTEGER.
  !           On entry,  N  specifies the number  of columns of the matrix
  !           op( B ) and the number of columns of the matrix C. N must be
  !           at least zero.
  !           Unchanged on exit.
  !
  !  K      - INTEGER.
  !           On entry,  K  specifies  the number of columns of the matrix
  !           op( A ) and the number of rows of the matrix op( B ). K must
  !           be at least  zero.
  !           Unchanged on exit.
  !
  !  ALPHA  - COMPLEX         .
  !           On entry, ALPHA specifies the scalar alpha.
  !           Unchanged on exit.
  !
  !  A      - COMPLEX          array of DIMENSION ( LDA, ka ), where ka is
  !           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
  !           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
  !           part of the array  A  must contain the matrix  A,  otherwise
  !           the leading  k by m  part of the array  A  must contain  the
  !           matrix A.
  !           Unchanged on exit.
  !
  !  LDA    - INTEGER.
  !           On entry, LDA specifies the first dimension of A as declared
  !           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
  !           LDA must be at least  max( 1, m ), otherwise  LDA must be at
  !           least  max( 1, k ).
  !           Unchanged on exit.
  !
  !  B      - COMPLEX          array of DIMENSION ( LDB, kb ), where kb is
  !           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
  !           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
  !           part of the array  B  must contain the matrix  B,  otherwise
  !           the leading  n by k  part of the array  B  must contain  the
  !           matrix B.
  !           Unchanged on exit.
  !
  !  LDB    - INTEGER.
  !           On entry, LDB specifies the first dimension of B as declared
  !           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
  !           LDB must be at least  max( 1, k ), otherwise  LDB must be at
  !           least  max( 1, n ).
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
  !           On exit, the array  C  is overwritten by the  m by n  matrix
  !           ( alpha*op( A )*op( B ) + beta*C ).
  !
  !  LDC    - INTEGER.
  !           On entry, LDC specifies the first dimension of C as declared
  !           in  the  calling  (sub)  program.   LDC  must  be  at  least
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
  !***END PROLOGUE  CGEMM
  !     .. Scalar Arguments ..
  CHARACTER :: Transa , Transb
  INTEGER M , N , K , Lda , Ldb , Ldc
  COMPLEX Alpha , Beta
  !     .. Array Arguments ..
  COMPLEX A(Lda,*) , B(Ldb,*) , C(Ldc,*)
  !     .. External Functions ..
  LOGICAL LSAME
  EXTERNAL LSAME
  !     .. External Subroutines ..
  EXTERNAL XERBLA
  !     .. Intrinsic Functions ..
  INTRINSIC CONJG , MAX
  !     .. Local Scalars ..
  LOGICAL conja , conjb , nota , notb
  INTEGER i , info , j , l , ncola , nrowa , nrowb
  COMPLEX temp
  !     .. Parameters ..
  COMPLEX ONE
  PARAMETER (ONE=(1.0E+0,0.0E+0))
  COMPLEX ZERO
  PARAMETER (ZERO=(0.0E+0,0.0E+0))
  !***FIRST EXECUTABLE STATEMENT  CGEMM
  !
  !     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
  !     conjugated or transposed, set  CONJA and CONJB  as true if  A  and
  !     B  respectively are to be  transposed but  not conjugated  and set
  !     NROWA, NCOLA and  NROWB  as the number of rows and  columns  of  A
  !     and the number of rows of  B  respectively.
  !
  nota = LSAME(Transa,'N')
  notb = LSAME(Transb,'N')
  conja = LSAME(Transa,'C')
  conjb = LSAME(Transb,'C')
  IF ( nota ) THEN
    nrowa = M
    ncola = K
  ELSE
    nrowa = K
    ncola = M
  ENDIF
  IF ( notb ) THEN
    nrowb = K
  ELSE
    nrowb = N
  ENDIF
  !
  !     Test the input parameters.
  !
  info = 0
  IF ( (.NOT.nota).AND.(.NOT.conja).AND.(.NOT.LSAME(Transa,'T')) ) THEN
    info = 1
  ELSEIF ( (.NOT.notb).AND.(.NOT.conjb).AND.(.NOT.LSAME(Transb,'T')) ) THEN
    info = 2
  ELSEIF ( M<0 ) THEN
    info = 3
  ELSEIF ( N<0 ) THEN
    info = 4
  ELSEIF ( K<0 ) THEN
    info = 5
  ELSEIF ( Lda<MAX(1,nrowa) ) THEN
    info = 8
  ELSEIF ( Ldb<MAX(1,nrowb) ) THEN
    info = 10
  ELSEIF ( Ldc<MAX(1,M) ) THEN
    info = 13
  ENDIF
  IF ( info/=0 ) THEN
    CALL XERBLA('CGEMM ',info)
    RETURN
  ENDIF
  !
  !     Quick return if possible.
  !
  IF ( (M==0).OR.(N==0).OR.(((Alpha==ZERO).OR.(K==0)).AND.(Beta==ONE)) )&
    RETURN
  !
  !     And when  alpha.eq.zero.
  !
  IF ( Alpha==ZERO ) THEN
    IF ( Beta==ZERO ) THEN
      DO j = 1 , N
        DO i = 1 , M
          C(i,j) = ZERO
        ENDDO
      ENDDO
    ELSE
      DO j = 1 , N
        DO i = 1 , M
          C(i,j) = Beta*C(i,j)
        ENDDO
      ENDDO
    ENDIF
    RETURN
  ENDIF
  !
  !     Start the operations.
  !
  IF ( notb ) THEN
    IF ( nota ) THEN
      !
      !           Form  C := alpha*A*B + beta*C.
      !
      DO j = 1 , N
        IF ( Beta==ZERO ) THEN
          DO i = 1 , M
            C(i,j) = ZERO
          ENDDO
        ELSEIF ( Beta/=ONE ) THEN
          DO i = 1 , M
            C(i,j) = Beta*C(i,j)
          ENDDO
        ENDIF
        DO l = 1 , K
          IF ( B(l,j)/=ZERO ) THEN
            temp = Alpha*B(l,j)
            DO i = 1 , M
              C(i,j) = C(i,j) + temp*A(i,l)
            ENDDO
          ENDIF
        ENDDO
      ENDDO
    ELSEIF ( conja ) THEN
      !
      !           Form  C := alpha*conjg( A' )*B + beta*C.
      !
      DO j = 1 , N
        DO i = 1 , M
          temp = ZERO
          DO l = 1 , K
            temp = temp + CONJG(A(l,i))*B(l,j)
          ENDDO
          IF ( Beta==ZERO ) THEN
            C(i,j) = Alpha*temp
          ELSE
            C(i,j) = Alpha*temp + Beta*C(i,j)
          ENDIF
        ENDDO
      ENDDO
    ELSE
      !
      !           Form  C := alpha*A'*B + beta*C
      !
      DO j = 1 , N
        DO i = 1 , M
          temp = ZERO
          DO l = 1 , K
            temp = temp + A(l,i)*B(l,j)
          ENDDO
          IF ( Beta==ZERO ) THEN
            C(i,j) = Alpha*temp
          ELSE
            C(i,j) = Alpha*temp + Beta*C(i,j)
          ENDIF
        ENDDO
      ENDDO
    ENDIF
  ELSEIF ( nota ) THEN
    IF ( conjb ) THEN
      !
      !           Form  C := alpha*A*conjg( B' ) + beta*C.
      !
      DO j = 1 , N
        IF ( Beta==ZERO ) THEN
          DO i = 1 , M
            C(i,j) = ZERO
          ENDDO
        ELSEIF ( Beta/=ONE ) THEN
          DO i = 1 , M
            C(i,j) = Beta*C(i,j)
          ENDDO
        ENDIF
        DO l = 1 , K
          IF ( B(j,l)/=ZERO ) THEN
            temp = Alpha*CONJG(B(j,l))
            DO i = 1 , M
              C(i,j) = C(i,j) + temp*A(i,l)
            ENDDO
          ENDIF
        ENDDO
      ENDDO
    ELSE
      !
      !           Form  C := alpha*A*B'          + beta*C
      !
      DO j = 1 , N
        IF ( Beta==ZERO ) THEN
          DO i = 1 , M
            C(i,j) = ZERO
          ENDDO
        ELSEIF ( Beta/=ONE ) THEN
          DO i = 1 , M
            C(i,j) = Beta*C(i,j)
          ENDDO
        ENDIF
        DO l = 1 , K
          IF ( B(j,l)/=ZERO ) THEN
            temp = Alpha*B(j,l)
            DO i = 1 , M
              C(i,j) = C(i,j) + temp*A(i,l)
            ENDDO
          ENDIF
        ENDDO
      ENDDO
    ENDIF
  ELSEIF ( conja ) THEN
    IF ( conjb ) THEN
      !
      !           Form  C := alpha*conjg( A' )*conjg( B' ) + beta*C.
      !
      DO j = 1 , N
        DO i = 1 , M
          temp = ZERO
          DO l = 1 , K
            temp = temp + CONJG(A(l,i))*CONJG(B(j,l))
          ENDDO
          IF ( Beta==ZERO ) THEN
            C(i,j) = Alpha*temp
          ELSE
            C(i,j) = Alpha*temp + Beta*C(i,j)
          ENDIF
        ENDDO
      ENDDO
    ELSE
      !
      !           Form  C := alpha*conjg( A' )*B' + beta*C
      !
      DO j = 1 , N
        DO i = 1 , M
          temp = ZERO
          DO l = 1 , K
            temp = temp + CONJG(A(l,i))*B(j,l)
          ENDDO
          IF ( Beta==ZERO ) THEN
            C(i,j) = Alpha*temp
          ELSE
            C(i,j) = Alpha*temp + Beta*C(i,j)
          ENDIF
        ENDDO
      ENDDO
    ENDIF
  ELSEIF ( conjb ) THEN
    !
    !           Form  C := alpha*A'*conjg( B' ) + beta*C
    !
    DO j = 1 , N
      DO i = 1 , M
        temp = ZERO
        DO l = 1 , K
          temp = temp + A(l,i)*CONJG(B(j,l))
        ENDDO
        IF ( Beta==ZERO ) THEN
          C(i,j) = Alpha*temp
        ELSE
          C(i,j) = Alpha*temp + Beta*C(i,j)
        ENDIF
      ENDDO
    ENDDO
  ELSE
    !
    !           Form  C := alpha*A'*B' + beta*C
    !
    DO j = 1 , N
      DO i = 1 , M
        temp = ZERO
        DO l = 1 , K
          temp = temp + A(l,i)*B(j,l)
        ENDDO
        IF ( Beta==ZERO ) THEN
          C(i,j) = Alpha*temp
        ELSE
          C(i,j) = Alpha*temp + Beta*C(i,j)
        ENDIF
      ENDDO
    ENDDO
  ENDIF
  !
  !
  !     End of CGEMM .
  !
END SUBROUTINE CGEMM
