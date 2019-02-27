!DECK IMTQLV
SUBROUTINE IMTQLV(N,D,E,E2,W,Ind,Ierr,Rv1)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  IMTQLV
  !***PURPOSE  Compute the eigenvalues of a symmetric tridiagonal matrix
  !            using the implicit QL method.  Eigenvectors may be computed
  !            later.
  !***LIBRARY   SLATEC (EISPACK)
  !***CATEGORY  D4A5, D4C2A
  !***TYPE      SINGLE PRECISION (IMTQLV-S)
  !***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
  !***AUTHOR  Smith, B. T., et al.
  !***DESCRIPTION
  !
  !     This subroutine is a variant of  IMTQL1  which is a translation of
  !     ALGOL procedure IMTQL1, NUM. MATH. 12, 377-383(1968) by Martin and
  !     Wilkinson, as modified in NUM. MATH. 15, 450(1970) by Dubrulle.
  !     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 241-248(1971).
  !
  !     This subroutine finds the eigenvalues of a SYMMETRIC TRIDIAGONAL
  !     matrix by the implicit QL method and associates with them
  !     their corresponding submatrix indices.
  !
  !     On INPUT
  !
  !        N is the order of the matrix.  N is an INTEGER variable.
  !
  !        D contains the diagonal elements of the symmetric tridiagonal
  !          matrix.  D is a one-dimensional REAL array, dimensioned D(N).
  !
  !        E contains the subdiagonal elements of the symmetric
  !          tridiagonal matrix in its last N-1 positions.  E(1) is
  !          arbitrary.  E is a one-dimensional REAL array, dimensioned
  !          E(N).
  !
  !        E2 contains the squares of the corresponding elements of E in
  !          its last N-1 positions.  E2(1) is arbitrary.  E2 is a one-
  !          dimensional REAL array, dimensioned E2(N).
  !
  !     On OUTPUT
  !
  !        D and E are unaltered.
  !
  !        Elements of E2, corresponding to elements of E regarded as
  !          negligible, have been replaced by zero causing the matrix to
  !          split into a direct sum of submatrices.  E2(1) is also set
  !          to zero.
  !
  !        W contains the eigenvalues in ascending order.  If an error
  !          exit is made, the eigenvalues are correct and ordered for
  !          indices 1, 2, ..., IERR-1, but may not be the smallest
  !          eigenvalues.  W is a one-dimensional REAL array, dimensioned
  !          W(N).
  !
  !        IND contains the submatrix indices associated with the
  !          corresponding eigenvalues in W -- 1 for eigenvalues belonging
  !          to the first submatrix from the top, 2 for those belonging to
  !          the second submatrix, etc.  IND is a one-dimensional REAL
  !          array, dimensioned IND(N).
  !
  !        IERR is an INTEGER flag set to
  !          Zero       for normal return,
  !          J          if the J-th eigenvalue has not been
  !                     determined after 30 iterations.
  !                     The eigenvalues should be correct for indices
  !                     1, 2, ..., IERR-1.  These eigenvalues are
  !                     ordered, but are not necessarily the smallest.
  !
  !        RV1 is a one-dimensional REAL array used for temporary storage,
  !          dimensioned RV1(N).
  !
  !     Calls PYTHAG(A,B) for sqrt(A**2 + B**2).
  !
  !     Questions and comments should be directed to B. S. Garbow,
  !     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
  !     ------------------------------------------------------------------
  !
  !***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
  !                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
  !                 system Routines - EISPACK Guide, Springer-Verlag,
  !                 1976.
  !***ROUTINES CALLED  PYTHAG
  !***REVISION HISTORY  (YYMMDD)
  !   760101  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  IMTQLV
  !
  INTEGER i, j, k, l, m, N, ii, mml, tag, Ierr
  REAL D(*), E(*), E2(*), W(*), Rv1(*)
  REAL b, c, f, g, p, r, s, s1, s2
  REAL PYTHAG
  INTEGER Ind(*)
  !
  !***FIRST EXECUTABLE STATEMENT  IMTQLV
  Ierr = 0
  k = 0
  tag = 0
  !
  DO i = 1, N
    W(i) = D(i)
    IF ( i/=1 ) Rv1(i-1) = E(i)
  ENDDO
  !
  E2(1) = 0.0E0
  Rv1(N) = 0.0E0
  !
  DO l = 1, N
    j = 0
    !     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
    50 CONTINUE
    DO m = l, N
      IF ( m==N ) EXIT
      s1 = ABS(W(m)) + ABS(W(m+1))
      s2 = s1 + ABS(Rv1(m))
      IF ( s2==s1 ) EXIT
      !     .......... GUARD AGAINST UNDERFLOWED ELEMENT OF E2 ..........
      IF ( E2(m+1)==0.0E0 ) GOTO 100
    ENDDO
    !
    IF ( m<=k ) GOTO 150
    IF ( m/=N ) E2(m+1) = 0.0E0
    100    k = m
    tag = tag + 1
    150    p = W(l)
    IF ( m==l ) THEN
      !     .......... ORDER EIGENVALUES ..........
      IF ( l/=1 ) THEN
        !     .......... FOR I=L STEP -1 UNTIL 2 DO -- ..........
        DO ii = 2, l
          i = l + 2 - ii
          IF ( p>=W(i-1) ) GOTO 160
          W(i) = W(i-1)
          Ind(i) = Ind(i-1)
        ENDDO
      ENDIF
      !
      i = 1
      160      W(i) = p
      Ind(i) = tag
    ELSE
      IF ( j==30 ) GOTO 200
      j = j + 1
      !     .......... FORM SHIFT ..........
      g = (W(l+1)-p)/(2.0E0*Rv1(l))
      r = PYTHAG(g,1.0E0)
      g = W(m) - p + Rv1(l)/(g+SIGN(r,g))
      s = 1.0E0
      c = 1.0E0
      p = 0.0E0
      mml = m - l
      !     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
      DO ii = 1, mml
        i = m - ii
        f = s*Rv1(i)
        b = c*Rv1(i)
        IF ( ABS(f)<ABS(g) ) THEN
          s = f/g
          r = SQRT(s*s+1.0E0)
          Rv1(i+1) = g*r
          c = 1.0E0/r
          s = s*c
        ELSE
          c = g/f
          r = SQRT(c*c+1.0E0)
          Rv1(i+1) = f*r
          s = 1.0E0/r
          c = c*s
        ENDIF
        g = W(i+1) - p
        r = (W(i)-g)*s + 2.0E0*c*b
        p = s*r
        W(i+1) = g + p
        g = c*r - b
      ENDDO
      !
      W(l) = W(l) - p
      Rv1(l) = g
      Rv1(m) = 0.0E0
      GOTO 50
    ENDIF
  ENDDO
  !
  GOTO 99999
  !     .......... SET ERROR -- NO CONVERGENCE TO AN
  !                EIGENVALUE AFTER 30 ITERATIONS ..........
  200  Ierr = l
  99999 CONTINUE
END SUBROUTINE IMTQLV
