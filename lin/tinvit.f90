!DECK TINVIT
SUBROUTINE TINVIT(Nm,N,D,E,E2,M,W,Ind,Z,Ierr,Rv1,Rv2,Rv3,Rv4,Rv6)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  TINVIT
  !***PURPOSE  Compute the eigenvectors of symmetric tridiagonal matrix
  !            corresponding to specified eigenvalues, using inverse
  !            iteration.
  !***LIBRARY   SLATEC (EISPACK)
  !***CATEGORY  D4C3
  !***TYPE      SINGLE PRECISION (TINVIT-S)
  !***KEYWORDS  EIGENVECTORS, EISPACK
  !***AUTHOR  Smith, B. T., et al.
  !***DESCRIPTION
  !
  !     This subroutine is a translation of the inverse iteration tech-
  !     nique in the ALGOL procedure TRISTURM by Peters and Wilkinson.
  !     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971).
  !
  !     This subroutine finds those eigenvectors of a TRIDIAGONAL
  !     SYMMETRIC matrix corresponding to specified eigenvalues,
  !     using inverse iteration.
  !
  !     On Input
  !
  !        NM must be set to the row dimension of the two-dimensional
  !          array parameter, Z, as declared in the calling program
  !          dimension statement.  NM is an INTEGER variable.
  !
  !        N is the order of the matrix.  N is an INTEGER variable.
  !          N must be less than or equal to NM.
  !
  !        D contains the diagonal elements of the symmetric tridiagonal
  !          matrix.  D is a one-dimensional REAL array, dimensioned D(N).
  !
  !        E contains the subdiagonal elements of the symmetric
  !          tridiagonal matrix in its last N-1 positions.  E(1) is
  !          arbitrary.  E is a one-dimensional REAL array, dimensioned
  !          E(N).
  !
  !        E2 contains the squares of the corresponding elements of E,
  !          with zeros corresponding to negligible elements of E.
  !          E(I) is considered negligible if it is not larger than
  !          the product of the relative machine precision and the sum
  !          of the magnitudes of D(I) and D(I-1).  E2(1) must contain
  !          0.0e0 if the eigenvalues are in ascending order, or 2.0e0
  !          if the eigenvalues are in descending order.  If  BISECT,
  !          TRIDIB, or  IMTQLV  has been used to find the eigenvalues,
  !          their output E2 array is exactly what is expected here.
  !          E2 is a one-dimensional REAL array, dimensioned E2(N).
  !
  !        M is the number of specified eigenvalues for which eigenvectors
  !          are to be determined.  M is an INTEGER variable.
  !
  !        W contains the M eigenvalues in ascending or descending order.
  !          W is a one-dimensional REAL array, dimensioned W(M).
  !
  !        IND contains in its first M positions the submatrix indices
  !          associated with the corresponding eigenvalues in W --
  !          1 for eigenvalues belonging to the first submatrix from
  !          the top, 2 for those belonging to the second submatrix, etc.
  !          If  BISECT  or  TRIDIB  has been used to determine the
  !          eigenvalues, their output IND array is suitable for input
  !          to TINVIT.  IND is a one-dimensional INTEGER array,
  !          dimensioned IND(M).
  !
  !     On Output
  !
  !       ** All input arrays are unaltered.**
  !
  !        Z contains the associated set of orthonormal eigenvectors.
  !          Any vector which fails to converge is set to zero.
  !          Z is a two-dimensional REAL array, dimensioned Z(NM,M).
  !
  !        IERR is an INTEGER flag set to
  !          Zero       for normal return,
  !          -J         if the eigenvector corresponding to the J-th
  !                     eigenvalue fails to converge in 5 iterations.
  !
  !        RV1, RV2 and RV3 are one-dimensional REAL arrays used for
  !          temporary storage.  They are used to store the main diagonal
  !          and the two adjacent diagonals of the triangular matrix
  !          produced in the inverse iteration process.  RV1, RV2 and
  !          RV3 are dimensioned RV1(N), RV2(N) and RV3(N).
  !
  !        RV4 and RV6 are one-dimensional REAL arrays used for temporary
  !          storage.  RV4 holds the multipliers of the Gaussian
  !          elimination process.  RV6 holds the approximate eigenvectors
  !          in this process.  RV4 and RV6 are dimensioned RV4(N) and
  !          RV6(N).
  !
  !     Questions and comments should be directed to B. S. Garbow,
  !     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
  !     ------------------------------------------------------------------
  !
  !***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
  !                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
  !                 system Routines - EISPACK Guide, Springer-Verlag,
  !                 1976.
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   760101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  TINVIT
  !
  INTEGER i, j, M, N, p, q, r, s, ii, ip, jj, Nm, its, tag, &
    Ierr, group
  INTEGER Ind(*)
  REAL D(*), E(*), E2(*), W(*), Z(Nm,*)
  REAL Rv1(*), Rv2(*), Rv3(*), Rv4(*), Rv6(*)
  REAL u, v, uk, xu, x0, x1, eps2, eps3, eps4, norm, order
  !
  !***FIRST EXECUTABLE STATEMENT  TINVIT
  Ierr = 0
  IF ( M==0 ) GOTO 99999
  tag = 0
  order = 1.0E0 - E2(1)
  q = 0
  !     .......... ESTABLISH AND PROCESS NEXT SUBMATRIX ..........
  100  p = q + 1
  !
  DO q = p, N
    IF ( q==N ) EXIT
    IF ( E2(q+1)==0.0E0 ) EXIT
  ENDDO
  !     .......... FIND VECTORS BY INVERSE ITERATION ..........
  tag = tag + 1
  s = 0
  !
  DO r = 1, M
    IF ( Ind(r)/=tag ) CYCLE
    its = 1
    x1 = W(r)
    IF ( s==0 ) THEN
      !     .......... CHECK FOR ISOLATED ROOT ..........
      xu = 1.0E0
      IF ( p/=q ) THEN
        norm = ABS(D(p))
        ip = p + 1
        !
        DO i = ip, q
          norm = MAX(norm,ABS(D(i))+ABS(E(i)))
        ENDDO
        !     .......... EPS2 IS THE CRITERION FOR GROUPING,
        !                EPS3 REPLACES ZERO PIVOTS AND EQUAL
        !                ROOTS ARE MODIFIED BY EPS3,
        !                EPS4 IS TAKEN VERY SMALL TO AVOID OVERFLOW ..........
        eps2 = 1.0E-3*norm
        eps3 = norm
        DO
          eps3 = 0.5E0*eps3
          IF ( norm+eps3<=norm ) THEN
            uk = SQRT(REAL(q-p+5))
            eps3 = uk*eps3
            eps4 = uk*eps3
            uk = eps4/uk
            s = p
            group = 0
            EXIT
          ENDIF
        ENDDO
      ELSE
        Rv6(p) = 1.0E0
        GOTO 150
      ENDIF
      !     .......... LOOK FOR CLOSE OR COINCIDENT ROOTS ..........
    ELSEIF ( ABS(x1-x0)>=eps2 ) THEN
      group = 0
    ELSE
      group = group + 1
      IF ( order*(x1-x0)<=0.0E0 ) x1 = x0 + order*eps3
    ENDIF
    !     .......... ELIMINATION WITH INTERCHANGES AND
    !                INITIALIZATION OF VECTOR ..........
    v = 0.0E0
    !
    DO i = p, q
      Rv6(i) = uk
      IF ( i/=p ) THEN
        IF ( ABS(E(i))<ABS(u) ) THEN
          xu = E(i)/u
          Rv4(i) = xu
          Rv1(i-1) = u
          Rv2(i-1) = v
          Rv3(i-1) = 0.0E0
        ELSE
          !     .......... WARNING -- A DIVIDE CHECK MAY OCCUR HERE IF
          !                E2 ARRAY HAS NOT BEEN SPECIFIED CORRECTLY ..........
          xu = u/E(i)
          Rv4(i) = xu
          Rv1(i-1) = E(i)
          Rv2(i-1) = D(i) - x1
          Rv3(i-1) = 0.0E0
          IF ( i/=q ) Rv3(i-1) = E(i+1)
          u = v - xu*Rv2(i-1)
          v = -xu*Rv3(i-1)
          CYCLE
        ENDIF
      ENDIF
      u = D(i) - x1 - xu*v
      IF ( i/=q ) v = E(i+1)
    ENDDO
    !
    IF ( u==0.0E0 ) u = eps3
    Rv1(q) = u
    Rv2(q) = 0.0E0
    Rv3(q) = 0.0E0
    DO
      !     .......... BACK SUBSTITUTION
      !                FOR I=Q STEP -1 UNTIL P DO -- ..........
      DO ii = p, q
        i = p + q - ii
        Rv6(i) = (Rv6(i)-u*Rv2(i)-v*Rv3(i))/Rv1(i)
        v = u
        u = Rv6(i)
      ENDDO
      !     .......... ORTHOGONALIZE WITH RESPECT TO PREVIOUS
      !                MEMBERS OF GROUP ..........
      IF ( group/=0 ) THEN
        j = r
        !
        DO jj = 1, group
          DO
            j = j - 1
            IF ( Ind(j)==tag ) THEN
              xu = 0.0E0
              !
              DO i = p, q
                xu = xu + Rv6(i)*Z(i,j)
              ENDDO
              !
              DO i = p, q
                Rv6(i) = Rv6(i) - xu*Z(i,j)
              ENDDO
              EXIT
            ENDIF
          ENDDO
          !
        ENDDO
      ENDIF
      !
      norm = 0.0E0
      !
      DO i = p, q
        norm = norm + ABS(Rv6(i))
      ENDDO
      !
      IF ( norm>=1.0E0 ) THEN
        !     .......... NORMALIZE SO THAT SUM OF SQUARES IS
        !                1 AND EXPAND TO FULL ORDER ..........
        u = 0.0E0
        !
        DO i = p, q
          u = u + Rv6(i)**2
        ENDDO
        !
        xu = 1.0E0/SQRT(u)
        EXIT
        !     .......... FORWARD SUBSTITUTION ..........
      ELSEIF ( its==5 ) THEN
        !     .......... SET ERROR -- NON-CONVERGED EIGENVECTOR ..........
        Ierr = -r
        xu = 0.0E0
        EXIT
      ELSE
        IF ( norm/=0.0E0 ) THEN
          xu = eps4/norm
          !
          DO i = p, q
            Rv6(i) = Rv6(i)*xu
          ENDDO
        ELSE
          Rv6(s) = eps4
          s = s + 1
          IF ( s>q ) s = p
        ENDIF
        !     .......... ELIMINATION OPERATIONS ON NEXT VECTOR
        !                ITERATE ..........
        DO i = ip, q
          u = Rv6(i)
          !     .......... IF RV1(I-1) .EQ. E(I), A ROW INTERCHANGE
          !                WAS PERFORMED EARLIER IN THE
          !                TRIANGULARIZATION PROCESS ..........
          IF ( Rv1(i-1)==E(i) ) THEN
            u = Rv6(i-1)
            Rv6(i-1) = Rv6(i)
          ENDIF
          Rv6(i) = u - Rv4(i)*Rv6(i-1)
        ENDDO
        !
        its = its + 1
      ENDIF
    ENDDO
    !
    150    DO i = 1, N
    Z(i,r) = 0.0E0
  ENDDO
  !
  DO i = p, q
    Z(i,r) = Rv6(i)*xu
  ENDDO
  !
  x0 = x1
ENDDO
!
IF ( q<N ) GOTO 100
  99999 CONTINUE
  END SUBROUTINE TINVIT
