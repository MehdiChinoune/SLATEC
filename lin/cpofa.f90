!DECK CPOFA
SUBROUTINE CPOFA(A,Lda,N,Info)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  CPOFA
  !***PURPOSE  Factor a complex Hermitian positive definite matrix.
  !***LIBRARY   SLATEC (LINPACK)
  !***CATEGORY  D2D1B
  !***TYPE      COMPLEX (SPOFA-S, DPOFA-D, CPOFA-C)
  !***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION,
  !             POSITIVE DEFINITE
  !***AUTHOR  Moler, C. B., (U. of New Mexico)
  !***DESCRIPTION
  !
  !     CPOFA factors a complex Hermitian positive definite matrix.
  !
  !     CPOFA is usually called by CPOCO, but it can be called
  !     directly with a saving in time if  RCOND  is not needed.
  !     (Time for CPOCO) = (1 + 18/N)*(Time for CPOFA) .
  !
  !     On Entry
  !
  !        A       COMPLEX(LDA, N)
  !                the Hermitian matrix to be factored.  Only the
  !                diagonal and upper triangle are used.
  !
  !        LDA     INTEGER
  !                the leading dimension of the array  A .
  !
  !        N       INTEGER
  !                the order of the matrix  A .
  !
  !     On Return
  !
  !        A       an upper triangular matrix  R  so that  A =
  !                CTRANS(R)*R where  CTRANS(R)  is the conjugate
  !                transpose.  The strict lower triangle is unaltered.
  !                If  INFO .NE. 0, the factorization is not complete.
  !
  !        INFO    INTEGER
  !                = 0  for normal return.
  !                = K  signals an error condition.  The leading minor
  !                     of order  K  is not positive definite.
  !
  !***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***ROUTINES CALLED  CDOTC
  !***REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  CPOFA
  INTEGER Lda, N, Info
  COMPLEX A(Lda,*)
  !
  COMPLEX CDOTC, t
  REAL s
  INTEGER j, jm1, k
  !***FIRST EXECUTABLE STATEMENT  CPOFA
  DO j = 1, N
    Info = j
    s = 0.0E0
    jm1 = j - 1
    IF ( jm1>=1 ) THEN
      DO k = 1, jm1
        t = A(k,j) - CDOTC(k-1,A(1,k),1,A(1,j),1)
        t = t/A(k,k)
        A(k,j) = t
        s = s + REAL(t*CONJG(t))
      ENDDO
    ENDIF
    s = REAL(A(j,j)) - s
    IF ( s<=0.0E0.OR.AIMAG(A(j,j))/=0.0E0 ) RETURN
    A(j,j) = CMPLX(SQRT(s),0.0E0)
  ENDDO
  Info = 0
  RETURN
END SUBROUTINE CPOFA
