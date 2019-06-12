!** SGBFA
SUBROUTINE SGBFA(Abd,Lda,N,Ml,Mu,Ipvt,Info)
  !>
  !  Factor a band matrix using Gaussian elimination.
  !***
  ! **Library:**   SLATEC (LINPACK)
  !***
  ! **Category:**  D2A2
  !***
  ! **Type:**      SINGLE PRECISION (SGBFA-S, DGBFA-D, CGBFA-C)
  !***
  ! **Keywords:**  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION
  !***
  ! **Author:**  Moler, C. B., (U. of New Mexico)
  !***
  ! **Description:**
  !
  !     SGBFA factors a real band matrix by elimination.
  !
  !     SGBFA is usually called by SBGCO, but it can be called
  !     directly with a saving in time if  RCOND  is not needed.
  !
  !     On Entry
  !
  !        ABD     REAL(LDA, N)
  !                contains the matrix in band storage.  The columns
  !                of the matrix are stored in the columns of  ABD  and
  !                the diagonals of the matrix are stored in rows
  !                ML+1 through 2*ML+MU+1 of  ABD .
  !                See the comments below for details.
  !
  !        LDA     INTEGER
  !                the leading dimension of the array  ABD .
  !                LDA must be .GE. 2*ML + MU + 1 .
  !
  !        N       INTEGER
  !                the order of the original matrix.
  !
  !        ML      INTEGER
  !                number of diagonals below the main diagonal.
  !                0 .LE. ML .LT. N .
  !
  !        MU      INTEGER
  !                number of diagonals above the main diagonal.
  !                0 .LE. MU .LT. N .
  !                More efficient if  ML .LE. MU .
  !     On Return
  !
  !        ABD     an upper triangular matrix in band storage and
  !                the multipliers which were used to obtain it.
  !                The factorization can be written  A = L*U, where
  !                L  is a product of permutation and unit lower
  !                triangular matrices and  U  is upper triangular.
  !
  !        IPVT    INTEGER(N)
  !                an integer vector of pivot indices.
  !
  !        INFO    INTEGER
  !                = 0  normal value.
  !                = K  if  U(K,K) .EQ. 0.0 .  This is not an error
  !                     condition for this subroutine, but it does
  !                     indicate that SGBSL will divide by zero if
  !                     called.  Use  RCOND  in SBGCO for a reliable
  !                     indication of singularity.
  !
  !     Band Storage
  !
  !           If  A  is a band matrix, the following program segment
  !           will set up the input.
  !
  !                   ML = (band width below the diagonal)
  !                   MU = (band width above the diagonal)
  !                   M = ML + MU + 1
  !                   DO 20 J = 1, N
  !                      I1 = MAX(1, J-MU)
  !                      I2 = MIN(N, J+ML)
  !                      DO 10 I = I1, I2
  !                         K = I - J + M
  !                         ABD(K,J) = A(I,J)
  !                10    CONTINUE
  !                20 CONTINUE
  !
  !           This uses rows  ML+1  through  2*ML+MU+1  of  ABD .
  !           In addition, the first  ML  rows in  ABD  are used for
  !           elements generated during the triangularization.
  !           The total number of rows needed in  ABD  is  2*ML+MU+1 .
  !           The  ML+MU by ML+MU  upper left triangle and the
  !           ML by ML  lower right triangle are not referenced.
  !
  !***
  ! **References:**  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***
  ! **Routines called:**  ISAMAX, SAXPY, SSCAL

  !* REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  USE blas, ONLY : SAXPY

  INTEGER :: Lda, N, Ml, Mu, Ipvt(N), Info
  REAL(SP) :: Abd(Lda,N)
  !
  REAL(SP) :: t
  INTEGER :: i, i0, j, ju, jz, j0, j1, k, kp1, l, lm, m, mm, nm1
  !
  !* FIRST EXECUTABLE STATEMENT  SGBFA
  m = Ml + Mu + 1
  Info = 0
  !
  !     ZERO INITIAL FILL-IN COLUMNS
  !
  j0 = Mu + 2
  j1 = MIN(N,m) - 1
  IF ( j1>=j0 ) THEN
    DO jz = j0, j1
      i0 = m + 1 - jz
      DO i = i0, Ml
        Abd(i,jz) = 0.0E0
      END DO
    END DO
  END IF
  jz = j1
  ju = 0
  !
  !     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
  !
  nm1 = N - 1
  IF ( nm1>=1 ) THEN
    DO k = 1, nm1
      kp1 = k + 1
      !
      !        ZERO NEXT FILL-IN COLUMN
      !
      jz = jz + 1
      IF ( jz<=N ) THEN
        IF ( Ml>=1 ) THEN
          DO i = 1, Ml
            Abd(i,jz) = 0.0E0
          END DO
        END IF
      END IF
      !
      !        FIND L = PIVOT INDEX
      !
      lm = MIN(Ml,N-k)
      l = MAXLOC( ABS(Abd(m:m+lm,k)),1) + m - 1
      Ipvt(k) = l + k - m
      !
      !        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
      !
      IF ( Abd(l,k)==0.0E0 ) THEN
        Info = k
      ELSE
        !
        !           INTERCHANGE IF NECESSARY
        !
        IF ( l/=m ) THEN
          t = Abd(l,k)
          Abd(l,k) = Abd(m,k)
          Abd(m,k) = t
        END IF
        !
        !           COMPUTE MULTIPLIERS
        !
        t = -1.0E0/Abd(m,k)
        Abd(m+1:m+lm,k) = t*Abd(m+1:m+lm,k)
        !
        !           ROW ELIMINATION WITH COLUMN INDEXING
        !
        ju = MIN(MAX(ju,Mu+Ipvt(k)),N)
        mm = m
        IF ( ju>=kp1 ) THEN
          DO j = kp1, ju
            l = l - 1
            mm = mm - 1
            t = Abd(l,j)
            IF ( l/=mm ) THEN
              Abd(l,j) = Abd(mm,j)
              Abd(mm,j) = t
            END IF
            CALL SAXPY(lm,t,Abd(m+1:m+lm,k),1,Abd(mm+1:mm+lm,j),1)
          END DO
        END IF
      END IF
    END DO
  END IF
  Ipvt(N) = N
  IF ( Abd(m,N)==0.0E0 ) Info = N
END SUBROUTINE SGBFA
