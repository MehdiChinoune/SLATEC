!*==DPBFA.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK DPBFA
      SUBROUTINE DPBFA(Abd,Lda,N,M,Info)
      IMPLICIT NONE
!*--DPBFA5
!***BEGIN PROLOGUE  DPBFA
!***PURPOSE  Factor a real symmetric positive definite matrix stored in
!            in band form.
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2B2
!***TYPE      DOUBLE PRECISION (SPBFA-S, DPBFA-D, CPBFA-C)
!***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION,
!             POSITIVE DEFINITE
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
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
!        ABD     an upper triangular matrix  R , stored in band
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
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  DDOT
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DPBFA
      INTEGER Lda , N , M , Info
      DOUBLE PRECISION Abd(Lda,*)
!
      DOUBLE PRECISION DDOT , t
      DOUBLE PRECISION s
      INTEGER ik , j , jk , k , mu
!***FIRST EXECUTABLE STATEMENT  DPBFA
      DO j = 1 , N
        Info = j
        s = 0.0D0
        ik = M + 1
        jk = MAX(j-M,1)
        mu = MAX(M+2-j,1)
        IF ( M>=mu ) THEN
          DO k = mu , M
            t = Abd(k,j) - DDOT(k-mu,Abd(ik,jk),1,Abd(mu,j),1)
            t = t/Abd(M+1,jk)
            Abd(k,j) = t
            s = s + t*t
            ik = ik - 1
            jk = jk + 1
          ENDDO
        ENDIF
        s = Abd(M+1,j) - s
        IF ( s<=0.0D0 ) GOTO 99999
        Abd(M+1,j) = SQRT(s)
      ENDDO
      Info = 0
99999 END SUBROUTINE DPBFA
