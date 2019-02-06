!*==CHERK.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK CHERK
      SUBROUTINE CHERK(Uplo,Trans,N,K,Alpha,A,Lda,Beta,C,Ldc)
      IMPLICIT NONE
!*--CHERK5
!***BEGIN PROLOGUE  CHERK
!***PURPOSE  Perform Hermitian rank k update of a complex Hermitian
!            matrix.
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1B6
!***TYPE      COMPLEX (SHERK-S, DHERK-D, CHERK-C)
!***KEYWORDS  LEVEL 3 BLAS, LINEAR ALGEBRA
!***AUTHOR  Dongarra, J., (ANL)
!           Duff, I., (AERE)
!           Du Croz, J., (NAG)
!           Hammarling, S. (NAG)
!***DESCRIPTION
!
!  CHERK  performs one of the hermitian rank k operations
!
!     C := alpha*A*conjg( A' ) + beta*C,
!
!  or
!
!     C := alpha*conjg( A' )*A + beta*C,
!
!  where  alpha and beta  are  real scalars,  C is an  n by n  hermitian
!  matrix and  A  is an  n by k  matrix in the  first case and a  k by n
!  matrix in the second case.
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
!              TRANS = 'N' or 'n'   C := alpha*A*conjg( A' ) + beta*C.
!
!              TRANS = 'C' or 'c'   C := alpha*conjg( A' )*A + beta*C.
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
!           of  columns   of  the   matrix   A,   and  on   entry   with
!           TRANS = 'C' or 'c',  K  specifies  the number of rows of the
!           matrix A.  K must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - REAL            .
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
!***END PROLOGUE  CHERK
!     .. Scalar Arguments ..
      CHARACTER*1 Uplo , Trans
      INTEGER N , K , Lda , Ldc
      REAL Alpha , Beta
!     .. Array Arguments ..
      COMPLEX A(Lda,*) , C(Ldc,*)
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     .. External Subroutines ..
      EXTERNAL XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC CMPLX , CONJG , MAX , REAL
!     .. Local Scalars ..
      LOGICAL upper
      INTEGER i , info , j , l , nrowa
      REAL rtemp
      COMPLEX temp
!     .. Parameters ..
      REAL ONE , ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
!***FIRST EXECUTABLE STATEMENT  CHERK
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
      ELSEIF ( Ldc<MAX(1,N) ) THEN
        info = 10
      ENDIF
      IF ( info/=0 ) THEN
        CALL XERBLA('CHERK ',info)
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
          IF ( Beta==ZERO ) THEN
            DO j = 1 , N
              DO i = 1 , j
                C(i,j) = ZERO
              ENDDO
            ENDDO
          ELSE
            DO j = 1 , N
              DO i = 1 , j - 1
                C(i,j) = Beta*C(i,j)
              ENDDO
              C(j,j) = Beta*REAL(C(j,j))
            ENDDO
          ENDIF
        ELSEIF ( Beta==ZERO ) THEN
          DO j = 1 , N
            DO i = j , N
              C(i,j) = ZERO
            ENDDO
          ENDDO
        ELSE
          DO j = 1 , N
            C(j,j) = Beta*REAL(C(j,j))
            DO i = j + 1 , N
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
!        Form  C := alpha*A*conjg( A' ) + beta*C.
!
        IF ( upper ) THEN
          DO j = 1 , N
            IF ( Beta==ZERO ) THEN
              DO i = 1 , j
                C(i,j) = ZERO
              ENDDO
            ELSEIF ( Beta/=ONE ) THEN
              DO i = 1 , j - 1
                C(i,j) = Beta*C(i,j)
              ENDDO
              C(j,j) = Beta*REAL(C(j,j))
            ENDIF
            DO l = 1 , K
              IF ( A(j,l)/=CMPLX(ZERO) ) THEN
                temp = Alpha*CONJG(A(j,l))
                DO i = 1 , j - 1
                  C(i,j) = C(i,j) + temp*A(i,l)
                ENDDO
                C(j,j) = REAL(C(j,j)) + REAL(temp*A(i,l))
              ENDIF
            ENDDO
          ENDDO
        ELSE
          DO j = 1 , N
            IF ( Beta==ZERO ) THEN
              DO i = j , N
                C(i,j) = ZERO
              ENDDO
            ELSEIF ( Beta/=ONE ) THEN
              C(j,j) = Beta*REAL(C(j,j))
              DO i = j + 1 , N
                C(i,j) = Beta*C(i,j)
              ENDDO
            ENDIF
            DO l = 1 , K
              IF ( A(j,l)/=CMPLX(ZERO) ) THEN
                temp = Alpha*CONJG(A(j,l))
                C(j,j) = REAL(C(j,j)) + REAL(temp*A(j,l))
                DO i = j + 1 , N
                  C(i,j) = C(i,j) + temp*A(i,l)
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDIF
!
!        Form  C := alpha*conjg( A' )*A + beta*C.
!
      ELSEIF ( upper ) THEN
        DO j = 1 , N
          DO i = 1 , j - 1
            temp = ZERO
            DO l = 1 , K
              temp = temp + CONJG(A(l,i))*A(l,j)
            ENDDO
            IF ( Beta==ZERO ) THEN
              C(i,j) = Alpha*temp
            ELSE
              C(i,j) = Alpha*temp + Beta*C(i,j)
            ENDIF
          ENDDO
          rtemp = ZERO
          DO l = 1 , K
            rtemp = rtemp + CONJG(A(l,j))*A(l,j)
          ENDDO
          IF ( Beta==ZERO ) THEN
            C(j,j) = Alpha*rtemp
          ELSE
            C(j,j) = Alpha*rtemp + Beta*REAL(C(j,j))
          ENDIF
        ENDDO
      ELSE
        DO j = 1 , N
          rtemp = ZERO
          DO l = 1 , K
            rtemp = rtemp + CONJG(A(l,j))*A(l,j)
          ENDDO
          IF ( Beta==ZERO ) THEN
            C(j,j) = Alpha*rtemp
          ELSE
            C(j,j) = Alpha*rtemp + Beta*REAL(C(j,j))
          ENDIF
          DO i = j + 1 , N
            temp = ZERO
            DO l = 1 , K
              temp = temp + CONJG(A(l,i))*A(l,j)
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
!     End of CHERK .
!
      END SUBROUTINE CHERK
