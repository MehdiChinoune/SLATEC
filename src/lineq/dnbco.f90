!** DNBCO
SUBROUTINE DNBCO(Abe,Lda,N,Ml,Mu,Ipvt,Rcond,Z)
  IMPLICIT NONE
  !>
  !***
  !  Factor a band matrix using Gaussian elimination and
  !            estimate the condition number.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  D2A2
  !***
  ! **Type:**      DOUBLE PRECISION (SNBCO-S, DNBCO-D, CNBCO-C)
  !***
  ! **Keywords:**  BANDED, LINEAR EQUATIONS, MATRIX FACTORIZATION,
  !             NONSYMMETRIC
  !***
  ! **Author:**  Voorhees, E. A., (LANL)
  !***
  ! **Description:**
  !
  !     DNBCO factors a double precision band matrix by Gaussian
  !     elimination and estimates the condition of the matrix.
  !
  !     If RCOND is not needed, DNBFA is slightly faster.
  !     To solve  A*X = B, follow DNBCO by DNBSL.
  !     To compute  INVERSE(A)*C, follow DNBCO by DNBSL.
  !     To compute  DETERMINANT(A), follow DNBCO by DNBDI.
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
  !        RCOND   DOUBLE PRECISION
  !                an estimate of the reciprocal condition of  A .
  !                For the system  A*X = B, relative perturbations
  !                in  A  and  B  of size  EPSILON  may cause
  !                relative perturbations in  X  of size  EPSILON/RCOND .
  !                If  RCOND  is so small that the logical expression
  !                         1.0 + RCOND .EQ. 1.0
  !                is true, then  A  may be singular to working
  !                precision.  In particular,  RCOND  is zero  if
  !                exact singularity is detected or the estimate
  !                underflows.
  !
  !        Z       DOUBLE PRECISION(N)
  !                a work vector whose contents are usually unimportant.
  !                If  A  is close to a singular matrix, then  Z  is
  !                an approximate null vector in the sense that
  !                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
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
  ! **Routines called:**  DASUM, DAXPY, DDOT, DNBFA, DSCAL

  !* REVISION HISTORY  (YYMMDD)
  !   800728  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  
  INTEGER Lda, N, Ml, Mu, Ipvt(*)
  REAL(8) :: Abe(Lda,*), Z(*)
  REAL(8) :: Rcond
  !
  REAL(8) :: DDOT, ek, t, wk, wkm
  REAL(8) :: anorm, s, DASUM, sm, ynorm
  INTEGER i, info, j, ju, k, kb, kp1, l, ldb, lm, lz, m, ml1, mm, nl, nu
  !* FIRST EXECUTABLE STATEMENT  DNBCO
  ml1 = Ml + 1
  ldb = Lda - 1
  anorm = 0.0D0
  DO j = 1, N
    nu = MIN(Mu,j-1)
    nl = MIN(Ml,N-j)
    l = 1 + nu + nl
    anorm = MAX(anorm,DASUM(l,Abe(j+nl,ml1-nl),ldb))
  ENDDO
  !
  !     FACTOR
  !
  CALL DNBFA(Abe,Lda,N,Ml,Mu,Ipvt,info)
  !
  !     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
  !     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND TRANS(A)*Y = E .
  !     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
  !     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF  W WHERE
  !     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
  !     OVERFLOW.
  !
  !     SOLVE TRANS(U)*W = E
  !
  ek = 1.0D0
  DO j = 1, N
    Z(j) = 0.0D0
  ENDDO
  m = Ml + Mu + 1
  ju = 0
  DO k = 1, N
    IF ( Z(k)/=0.0D0 ) ek = SIGN(ek,-Z(k))
    IF ( ABS(ek-Z(k))>ABS(Abe(k,ml1)) ) THEN
      s = ABS(Abe(k,ml1))/ABS(ek-Z(k))
      CALL DSCAL(N,s,Z,1)
      ek = s*ek
    ENDIF
    wk = ek - Z(k)
    wkm = -ek - Z(k)
    s = ABS(wk)
    sm = ABS(wkm)
    IF ( Abe(k,ml1)==0.0D0 ) THEN
      wk = 1.0D0
      wkm = 1.0D0
    ELSE
      wk = wk/Abe(k,ml1)
      wkm = wkm/Abe(k,ml1)
    ENDIF
    kp1 = k + 1
    ju = MIN(MAX(ju,Mu+Ipvt(k)),N)
    mm = ml1
    IF ( kp1<=ju ) THEN
      DO i = kp1, ju
        mm = mm + 1
        sm = sm + ABS(Z(i)+wkm*Abe(k,mm))
        Z(i) = Z(i) + wk*Abe(k,mm)
        s = s + ABS(Z(i))
      ENDDO
      IF ( s<sm ) THEN
        t = wkm - wk
        wk = wkm
        mm = ml1
        DO i = kp1, ju
          mm = mm + 1
          Z(i) = Z(i) + t*Abe(k,mm)
        ENDDO
      ENDIF
    ENDIF
    Z(k) = wk
  ENDDO
  s = 1.0D0/DASUM(N,Z,1)
  CALL DSCAL(N,s,Z,1)
  !
  !     SOLVE TRANS(L)*Y = W
  !
  DO kb = 1, N
    k = N + 1 - kb
    nl = MIN(Ml,N-k)
    IF ( k<N ) Z(k) = Z(k) + DDOT(nl,Abe(k+nl,ml1-nl),-ldb,Z(k+1),1)
    IF ( ABS(Z(k))>1.0D0 ) THEN
      s = 1.0D0/ABS(Z(k))
      CALL DSCAL(N,s,Z,1)
    ENDIF
    l = Ipvt(k)
    t = Z(l)
    Z(l) = Z(k)
    Z(k) = t
  ENDDO
  s = 1.0D0/DASUM(N,Z,1)
  CALL DSCAL(N,s,Z,1)
  !
  ynorm = 1.0D0
  !
  !     SOLVE L*V = Y
  !
  DO k = 1, N
    l = Ipvt(k)
    t = Z(l)
    Z(l) = Z(k)
    Z(k) = t
    nl = MIN(Ml,N-k)
    IF ( k<N ) CALL DAXPY(nl,t,Abe(k+nl,ml1-nl),-ldb,Z(k+1),1)
    IF ( ABS(Z(k))>1.0D0 ) THEN
      s = 1.0D0/ABS(Z(k))
      CALL DSCAL(N,s,Z,1)
      ynorm = s*ynorm
    ENDIF
  ENDDO
  s = 1.0D0/DASUM(N,Z,1)
  CALL DSCAL(N,s,Z,1)
  ynorm = s*ynorm
  !
  !     SOLVE  U*Z = V
  !
  DO kb = 1, N
    k = N + 1 - kb
    IF ( ABS(Z(k))>ABS(Abe(k,ml1)) ) THEN
      s = ABS(Abe(k,ml1))/ABS(Z(k))
      CALL DSCAL(N,s,Z,1)
      ynorm = s*ynorm
    ENDIF
    IF ( Abe(k,ml1)/=0.0D0 ) Z(k) = Z(k)/Abe(k,ml1)
    IF ( Abe(k,ml1)==0.0D0 ) Z(k) = 1.0D0
    lm = MIN(k,m) - 1
    lz = k - lm
    t = -Z(k)
    IF(k<=1) CYCLE
    CALL DAXPY(lm,t,Abe(k-1,Ml+2),-ldb,Z(lz),1)
  ENDDO
  !     MAKE ZNORM = 1.0D0
  s = 1.0D0/DASUM(N,Z,1)
  CALL DSCAL(N,s,Z,1)
  ynorm = s*ynorm
  !
  IF ( anorm/=0.0D0 ) Rcond = ynorm/anorm
  IF ( anorm==0.0D0 ) Rcond = 0.0D0
END SUBROUTINE DNBCO
