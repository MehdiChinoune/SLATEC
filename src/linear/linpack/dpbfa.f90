!** DPBFA
SUBROUTINE DPBFA(Abd,Lda,N,M,Info)
  !>
  !***
  !  Factor a real symmetric positive definite matrix stored in
  !            in band form.
  !***
  ! **Library:**   SLATEC (LINPACK)
  !***
  ! **Category:**  D2B2
  !***
  ! **Type:**      DOUBLE PRECISION (SPBFA-S, DPBFA-D, CPBFA-C)
  !***
  ! **Keywords:**  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION,
  !             POSITIVE DEFINITE
  !***
  ! **Author:**  Moler, C. B., (U. of New Mexico)
  !***
  ! **Description:**
  !
  !     DPBFA factors a double precision symmetric positive definite
  !     matrix stored in band form.
  !
  !     DPBFA is usually called by DPBCO, but it can be called
  !     directly with a saving in time if  RCOND  is not needed.
  !
  !     On Entry
  !
  !        ABD     DOUBLE PRECISION(LDA, N)
  !                the matrix to be factored.  The columns of the upper
  !                triangle are stored in the columns of ABD and the
  !                diagonals of the upper triangle are stored in the
  !                rows of ABD .  See the comments below for details.
  !
  !        LDA     INTEGER
  !                the leading dimension of the array  ABD .
  !                LDA must be .GE. M + 1 .
  !
  !        N       INTEGER
  !                the order of the matrix  A .
  !
  !        M       INTEGER
  !                the number of diagonals above the main diagonal.
  !                0 .LE. M .LT. N .
  !
  !     On Return
  !
  !        ABD     an upper triangular matrix  R, stored in band
  !                form, so that  A = TRANS(R)*R .
  !
  !        INFO    INTEGER
  !                = 0  for normal return.
  !                = K  if the leading minor of order  K  is not
  !                     positive definite.
  !
  !     Band Storage
  !
  !           If  A  is a symmetric positive definite band matrix,
  !           the following program segment will set up the input.
  !
  !                   M = (band width above diagonal)
  !                   DO 20 J = 1, N
  !                      I1 = MAX(1, J-M)
  !                      DO 10 I = I1, J
  !                         K = I-J+M+1
  !                         ABD(K,J) = A(I,J)
  !                10    CONTINUE
  !                20 CONTINUE
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

  INTEGER Lda, N, M, Info
  REAL(8) :: Abd(Lda,*)
  !
  REAL(8) :: t
  REAL(8) :: s
  INTEGER ik, j, jk, k, mu
  !* FIRST EXECUTABLE STATEMENT  DPBFA
  DO j = 1, N
    Info = j
    s = 0.0D0
    ik = M + 1
    jk = MAX(j-M,1)
    mu = MAX(M+2-j,1)
    IF ( M>=mu ) THEN
      DO k = mu, M
        t = Abd(k,j) - DDOT(k-mu,Abd(ik,jk),1,Abd(mu,j),1)
        t = t/Abd(M+1,jk)
        Abd(k,j) = t
        s = s + t*t
        ik = ik - 1
        jk = jk + 1
      END DO
    END IF
    s = Abd(M+1,j) - s
    IF ( s<=0.0D0 ) RETURN
    Abd(M+1,j) = SQRT(s)
  END DO
  Info = 0
  RETURN
END SUBROUTINE DPBFA
