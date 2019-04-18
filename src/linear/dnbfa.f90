!** DNBFA
SUBROUTINE DNBFA(Abe,Lda,N,Ml,Mu,Ipvt,Info)
  !>
  !  Factor a band matrix by elimination.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  D2A2
  !***
  ! **Type:**      DOUBLE PRECISION (SNBFA-S, DNBFA-D, CNBFA-C)
  !***
  ! **Keywords:**  BANDED, LINEAR EQUATIONS, MATRIX FACTORIZATION,
  !             NONSYMMETRIC
  !***
  ! **Author:**  Voorhees, E. A., (LANL)
  !***
  ! **Description:**
  !
  !     DNBFA factors a double precision band matrix by elimination.
  !
  !     DNBFA is usually called by DNBCO, but it can be called
  !     directly with a saving in time if RCOND is not needed.
  !
  !     On Entry
  !
  !        ABE     DOUBLE PRECISION(LDA, NC)
  !                contains the matrix in band storage.  The rows
  !                of the original matrix are stored in the rows
  !                of ABE and the diagonals of the original matrix
  !                are stored in columns 1 through ML+MU+1 of ABE.
  !                NC must be .GE. 2*ML+MU+1 .
  !                See the comments below for details.
  !
  !        LDA     INTEGER
  !                the leading dimension of the array ABE.
  !                LDA must be .GE.  N .
  !
  !        N       INTEGER
  !                the order of the original matrix.
  !
  !        ML      INTEGER
  !                number of diagonals below the main diagonal.
  !                0 .LE. ML .LT.  N .
  !
  !        MU      INTEGER
  !                number of diagonals above the main diagonal.
  !                0 .LE. MU .LT.  N .
  !                More efficient if ML .LE. MU .
  !
  !     On Return
  !
  !        ABE     an upper triangular matrix in band storage
  !                and the multipliers which were used to obtain it.
  !                The factorization can be written  A = L*U  where
  !                L is a product of permutation and unit lower
  !                triangular matrices and  U  is upper triangular.
  !
  !        IPVT    INTEGER(N)
  !                an integer vector of pivot indices.
  !
  !        INFO    INTEGER
  !                =0  normal value
  !                =K  if  U(K,K) .EQ. 0.0 .  This is not an error
  !                condition for this subroutine, but it does
  !                indicate that DNBSL will divide by zero if
  !                called.  Use RCOND in DNBCO for a reliable
  !                indication of singularity.
  !
  !     Band Storage
  !
  !           If  A  is a band matrix, the following program segment
  !           will set up the input.
  !
  !                   ML = (band width below the diagonal)
  !                   MU = (band width above the diagonal)
  !                   DO 20 I = 1, N
  !                      J1 = MAX(1, I-ML)
  !                      J2 = MIN(N, I+MU)
  !                      DO 10 J = J1, J2
  !                         K = J - I + ML + 1
  !                         ABE(I,K) = A(I,J)
  !                10    CONTINUE
  !                20 CONTINUE
  !
  !           This uses columns  1  through  ML+MU+1  of ABE .
  !           Furthermore,  ML  additional columns are needed in
  !           ABE  starting with column  ML+MU+2  for elements
  !           generated during the triangularization.  The total
  !           number of columns needed in  ABE  is  2*ML+MU+1 .
  !
  !     Example:  If the original matrix is
  !
  !           11 12 13  0  0  0
  !           21 22 23 24  0  0
  !            0 32 33 34 35  0
  !            0  0 43 44 45 46
  !            0  0  0 54 55 56
  !            0  0  0  0 65 66
  !
  !      then  N = 6, ML = 1, MU = 2, LDA .GE. 5  and ABE should contain
  !
  !            * 11 12 13  +    , * = not used
  !           21 22 23 24  +    , + = used for pivoting
  !           32 33 34 35  +
  !           43 44 45 46  +
  !           54 55 56  *  +
  !           65 66  *  *  +
  !
  !***
  ! **References:**  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***
  ! **Routines called:**  DAXPY, DSCAL, DSWAP, IDAMAX

  !* REVISION HISTORY  (YYMMDD)
  !   800728  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  INTEGER Lda, N, Ml, Mu, Ipvt(*), Info
  REAL(8) :: Abe(Lda,*)
  !
  INTEGER ml1, mb, m, n1, ldb, i, j, k, l, lm, lm1, lm2, mp
  REAL(8) :: t
  !* FIRST EXECUTABLE STATEMENT  DNBFA
  ml1 = Ml + 1
  mb = Ml + Mu
  m = Ml + Mu + 1
  n1 = N - 1
  ldb = Lda - 1
  Info = 0
  !
  !     SET FILL-IN COLUMNS TO ZERO
  !
  IF ( N>1 ) THEN
    IF ( Ml>0 ) THEN
      DO j = 1, Ml
        DO i = 1, N
          Abe(i,m+j) = 0.0D0
        END DO
      END DO
    END IF
    !
    !     GAUSSIAN ELIMINATION WITH PARTIAL ELIMINATION
    !
    DO k = 1, n1
      lm = MIN(N-k,Ml)
      lm1 = lm + 1
      lm2 = ml1 - lm
      !
      !     SEARCH FOR PIVOT INDEX
      !
      l = -IDAMAX(lm1,Abe(lm+k,lm2),ldb) + lm1 + k
      Ipvt(k) = l
      mp = MIN(mb,N-k)
      !
      !     SWAP ROWS IF NECESSARY
      !
      IF ( l/=k ) CALL DSWAP(mp+1,Abe(k,ml1),Lda,Abe(l,ml1+k-l),Lda)
      !
      !     SKIP COLUMN REDUCTION IF PIVOT IS ZERO
      !
      IF ( Abe(k,ml1)==0.0D0 ) THEN
        Info = k
      ELSE
        !
        !     COMPUTE MULTIPLIERS
        !
        t = -1.0/Abe(k,ml1)
        CALL DSCAL(lm,t,Abe(lm+k,lm2),ldb)
        !
        !     ROW ELIMINATION WITH COLUMN INDEXING
        !
        DO j = 1, mp
          CALL DAXPY(lm,Abe(k,ml1+j),Abe(lm+k,lm2),ldb,Abe(lm+k,lm2+j),ldb)
        END DO
      END IF
    END DO
  END IF
  Ipvt(N) = N
  IF ( Abe(N,ml1)==0.0D0 ) Info = N
END SUBROUTINE DNBFA
