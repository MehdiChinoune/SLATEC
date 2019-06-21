!** DPOFA
SUBROUTINE DPOFA(A,Lda,N,Info)
  !> Factor a real symmetric positive definite matrix.
  !***
  ! **Library:**   SLATEC (LINPACK)
  !***
  ! **Category:**  D2B1B
  !***
  ! **Type:**      DOUBLE PRECISION (SPOFA-S, DPOFA-D, CPOFA-C)
  !***
  ! **Keywords:**  LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION,
  !             POSITIVE DEFINITE
  !***
  ! **Author:**  Moler, C. B., (U. of New Mexico)
  !***
  ! **Description:**
  !
  !     DPOFA factors a double precision symmetric positive definite
  !     matrix.
  !
  !     DPOFA is usually called by DPOCO, but it can be called
  !     directly with a saving in time if  RCOND  is not needed.
  !     (time for DPOCO) = (1 + 18/N)*(time for DPOFA) .
  !
  !     On Entry
  !
  !        A       DOUBLE PRECISION(LDA, N)
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
  ! **Routines called:**  DDOT

  !* REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTEGER :: Lda, N, Info
  REAL(DP) :: A(Lda,N)
  !
  REAL(DP) :: t
  REAL(DP) :: s
  INTEGER :: j, jm1, k
  !* FIRST EXECUTABLE STATEMENT  DPOFA
  DO j = 1, N
    Info = j
    s = 0._DP
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
    IF( s<=0._DP ) RETURN
    A(j,j) = SQRT(s)
  END DO
  Info = 0
  RETURN
END SUBROUTINE DPOFA
