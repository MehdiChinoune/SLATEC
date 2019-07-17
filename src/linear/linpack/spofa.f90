!** SPOFA
PURE SUBROUTINE SPOFA(A,Lda,N,Info)
  !> Factor a real symmetric positive definite matrix.
  !***
  ! **Library:**   SLATEC (LINPACK)
  !***
  ! **Category:**  D2B1B
  !***
  ! **Type:**      SINGLE PRECISION (SPOFA-S, DPOFA-D, CPOFA-C)
  !***
  ! **Keywords:**  LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION, POSITIVE DEFINITE
  !***
  ! **Author:**  Moler, C. B., (U. of New Mexico)
  !***
  ! **Description:**
  !
  !     SPOFA factors a real symmetric positive definite matrix.
  !
  !     SPOFA is usually called by SPOCO, but it can be called
  !     directly with a saving in time if  RCOND  is not needed.
  !     (Time for SPOCO) = (1 + 18/N)*(Time for SPOFA) .
  !
  !     On Entry
  !
  !        A       REAL(LDA, N)
  !                the symmetric matrix to be factored.  Only the
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
  !        A       an upper triangular matrix  R  so that  A = TRANS(R)*R
  !                where  TRANS(R)  is the transpose.
  !                The strict lower triangle is unaltered.
  !                If  INFO /= 0, the factorization is not complete.
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
  ! **Routines called:**  SDOT

  !* REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTIONsection.  (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTEGER, INTENT(IN) :: Lda, N
  INTEGER, INTENT(OUT) :: Info
  REAL(SP), INTENT(INOUT) :: A(Lda,N)
  !
  REAL(SP) :: t
  REAL(SP) :: s
  INTEGER :: j, jm1, k
  !* FIRST EXECUTABLE STATEMENT  SPOFA
  DO j = 1, N
    Info = j
    s = 0._SP
    jm1 = j - 1
    IF( jm1>=1 ) THEN
      DO k = 1, jm1
        t = A(k,j) - DOT_PRODUCT(A(1:k-1,k),A(1:k-1,j))
        t = t/A(k,k)
        A(k,j) = t
        s = s + t*t
      END DO
    END IF
    s = A(j,j) - s
    IF( s<=0._SP ) RETURN
    A(j,j) = SQRT(s)
  END DO
  Info = 0

  RETURN
END SUBROUTINE SPOFA