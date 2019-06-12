!** CPOFA
SUBROUTINE CPOFA(A,Lda,N,Info)
  !>
  !  Factor a complex Hermitian positive definite matrix.
  !***
  ! **Library:**   SLATEC (LINPACK)
  !***
  ! **Category:**  D2D1B
  !***
  ! **Type:**      COMPLEX (SPOFA-S, DPOFA-D, CPOFA-C)
  !***
  ! **Keywords:**  LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION,
  !             POSITIVE DEFINITE
  !***
  ! **Author:**  Moler, C. B., (U. of New Mexico)
  !***
  ! **Description:**
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
  !***
  ! **References:**  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***
  ! **Routines called:**  CDOTC

  !* REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTEGER :: Lda, N, Info
  COMPLEX(SP) :: A(Lda,N)
  !
  COMPLEX(SP) :: t
  REAL(SP) :: s
  INTEGER :: j, jm1, k
  !* FIRST EXECUTABLE STATEMENT  CPOFA
  DO j = 1, N
    Info = j
    s = 0.0E0
    jm1 = j - 1
    IF ( jm1>=1 ) THEN
      DO k = 1, jm1
        t = A(k,j) - DOT_PRODUCT(A(1:k-1,k),A(1:k-1,j))
        t = t/A(k,k)
        A(k,j) = t
        s = s + REAL(t*CONJG(t))
      END DO
    END IF
    s = REAL(A(j,j)) - s
    IF ( s<=0.0E0.OR.AIMAG(A(j,j))/=0.0E0 ) RETURN
    A(j,j) = CMPLX(SQRT(s),0.0E0)
  END DO
  Info = 0
  RETURN
END SUBROUTINE CPOFA
