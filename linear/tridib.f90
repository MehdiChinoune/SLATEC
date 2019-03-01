!DECK TRIDIB
SUBROUTINE TRIDIB(N,Eps1,D,E,E2,Lb,Ub,M11,M,W,Ind,Ierr,Rv4,Rv5)
  IMPLICIT NONE
  REAL R1MACH
  !***BEGIN PROLOGUE  TRIDIB
  !***PURPOSE  Compute the eigenvalues of a symmetric tridiagonal matrix
  !            in a given interval using Sturm sequencing.
  !***LIBRARY   SLATEC (EISPACK)
  !***CATEGORY  D4A5, D4C2A
  !***TYPE      SINGLE PRECISION (TRIDIB-S)
  !***KEYWORDS  EIGENVALUES OF A REAL SYMMETRIC MATRIX, EISPACK
  !***AUTHOR  Smith, B. T., et al.
  !***DESCRIPTION
  !
  !     This subroutine is a translation of the ALGOL procedure BISECT,
  !     NUM. MATH. 9, 386-393(1967) by Barth, Martin, and Wilkinson.
  !     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 249-256(1971).
  !
  !     This subroutine finds those eigenvalues of a TRIDIAGONAL
  !     SYMMETRIC matrix between specified boundary indices,
  !     using bisection.
  !
  !     On Input
  !
  !        N is the order of the matrix.  N is an INTEGER variable.
  !
  !        EPS1 is an absolute error tolerance for the computed eigen-
  !          values.  If the input EPS1 is non-positive, it is reset for
  !          each submatrix to a default value, namely, minus the product
  !          of the relative machine precision and the 1-norm of the
  !          submatrix.  EPS1 is a REAL variable.
  !
  !        D contains the diagonal elements of the symmetric tridiagonal
  !          matrix.  D is a one-dimensional REAL array, dimensioned D(N).
  !
  !        E contains the subdiagonal elements of the symmetric
  !          tridiagonal matrix in its last N-1 positions.  E(1) is
  !          arbitrary.  E is a one-dimensional REAL array, dimensioned
  !          E(N).
  !
  !        E2 contains the squares of the corresponding elements of E.
  !          E2(1) is arbitrary.  E2 is a one-dimensional REAL array,
  !          dimensioned E2(N).
  !
  !        M11 specifies the lower boundary index for the set of desired
  !          eigenvalues.  M11 is an INTEGER variable.
  !
  !        M specifies the number of eigenvalues desired.  The upper
  !          boundary index M22 is then obtained as M22=M11+M-1.
  !          M is an INTEGER variable.
  !
  !     On Output
  !
  !        EPS1 is unaltered unless it has been reset to its
  !          (last) default value.
  !
  !        D and E are unaltered.
  !
  !        Elements of E2, corresponding to elements of E regarded
  !          as negligible, have been replaced by zero causing the
  !          matrix to split into a direct sum of submatrices.
  !          E2(1) is also set to zero.
  !
  !        LB and UB define an interval containing exactly the desired
  !          eigenvalues.  LB and UB are REAL variables.
  !
  !        W contains, in its first M positions, the eigenvalues
  !          between indices M11 and M22 in ascending order.
  !          W is a one-dimensional REAL array, dimensioned W(M).
  !
  !        IND contains in its first M positions the submatrix indices
  !          associated with the corresponding eigenvalues in W --
  !          1 for eigenvalues belonging to the first submatrix from
  !          the top, 2 for those belonging to the second submatrix, etc.
  !          IND is an one-dimensional INTEGER array, dimensioned IND(M).
  !
  !        IERR is an INTEGER flag set to
  !          Zero       for normal return,
  !          3*N+1      if multiple eigenvalues at index M11 make
  !                     unique selection of LB impossible,
  !          3*N+2      if multiple eigenvalues at index M22 make
  !                     unique selection of UB impossible.
  !
  !        RV4 and RV5 are one-dimensional REAL arrays used for temporary
  !          storage of the lower and upper bounds for the eigenvalues in
  !          the bisection process.  RV4 and RV5 are dimensioned RV4(N)
  !          and RV5(N).
  !
  !     Note that subroutine TQL1, IMTQL1, or TQLRAT is generally faster
  !     than TRIDIB, if more than N/4 eigenvalues are to be found.
  !
  !     Questions and comments should be directed to B. S. Garbow,
  !     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
  !     ------------------------------------------------------------------
  !
  !***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
  !                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
  !                 system Routines - EISPACK Guide, Springer-Verlag,
  !                 1976.
  !***ROUTINES CALLED  R1MACH
  !***REVISION HISTORY  (YYMMDD)
  !   760101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  TRIDIB
  !
  INTEGER i, j, k, l, M, N, p, q, r, s, ii, m1, m2, M11, m22, &
    tag, Ierr, isturm
  REAL D(*), E(*), E2(*), W(*), Rv4(*), Rv5(*)
  REAL u, v, Lb, t1, t2, Ub, xu, x0, x1, Eps1, machep, s1, s2
  INTEGER Ind(*)
  LOGICAL first
  !
  SAVE first, machep
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  TRIDIB
  IF ( first ) machep = R1MACH(4)
  first = .FALSE.
  !
  Ierr = 0
  tag = 0
  xu = D(1)
  x0 = D(1)
  u = 0.0E0
  !     .......... LOOK FOR SMALL SUB-DIAGONAL ENTRIES AND DETERMINE AN
  !                INTERVAL CONTAINING ALL THE EIGENVALUES ..........
  DO i = 1, N
    x1 = u
    u = 0.0E0
    IF ( i/=N ) u = ABS(E(i+1))
    xu = MIN(D(i)-(x1+u),xu)
    x0 = MAX(D(i)+(x1+u),x0)
    IF ( i/=1 ) THEN
      s1 = ABS(D(i)) + ABS(D(i-1))
      s2 = s1 + ABS(E(i))
      IF ( s2>s1 ) CYCLE
    ENDIF
    E2(i) = 0.0E0
  ENDDO
  !
  x1 = MAX(ABS(xu),ABS(x0))*machep*N
  xu = xu - x1
  t1 = xu
  x0 = x0 + x1
  t2 = x0
  !     .......... DETERMINE AN INTERVAL CONTAINING EXACTLY
  !                THE DESIRED EIGENVALUES ..........
  p = 1
  q = N
  m1 = M11 - 1
  IF ( m1==0 ) GOTO 200
  isturm = 1
  100  v = x1
  x1 = xu + (x0-xu)*0.5E0
  IF ( x1/=v ) GOTO 700
  !     .......... SET ERROR -- INTERVAL CANNOT BE FOUND CONTAINING
  !                EXACTLY THE DESIRED EIGENVALUES ..........
  Ierr = 3*N + isturm
  GOTO 1000
  200  m22 = m1 + M
  IF ( m22/=N ) THEN
    x0 = t2
    isturm = 2
    GOTO 100
  ENDIF
  300  q = 0
  r = 0
  !     .......... ESTABLISH AND PROCESS NEXT SUBMATRIX, REFINING
  !                INTERVAL BY THE GERSCHGORIN BOUNDS ..........
  400 CONTINUE
  IF ( r==M ) GOTO 1000
  tag = tag + 1
  p = q + 1
  xu = D(p)
  x0 = D(p)
  u = 0.0E0
  !
  DO q = p, N
    x1 = u
    u = 0.0E0
    v = 0.0E0
    IF ( q/=N ) THEN
      u = ABS(E(q+1))
      v = E2(q+1)
    ENDIF
    xu = MIN(D(q)-(x1+u),xu)
    x0 = MAX(D(q)+(x1+u),x0)
    IF ( v==0.0E0 ) EXIT
  ENDDO
  !
  x1 = MAX(ABS(xu),ABS(x0))*machep
  IF ( Eps1<=0.0E0 ) Eps1 = -x1
  IF ( p/=q ) THEN
    x1 = x1*(q-p+1)
    Lb = MAX(t1,xu-x1)
    Ub = MIN(t2,x0+x1)
    x1 = Lb
    isturm = 3
    GOTO 700
  ELSE
    !     .......... CHECK FOR ISOLATED ROOT WITHIN INTERVAL ..........
    IF ( t1>D(p).OR.D(p)>=t2 ) GOTO 900
    m1 = p
    m2 = p
    Rv5(p) = D(p)
    GOTO 800
  ENDIF
  500  xu = Lb
  !     .......... FOR I=K STEP -1 UNTIL M1 DO -- ..........
  DO ii = m1, k
    i = m1 + k - ii
    IF ( xu<Rv4(i) ) THEN
      xu = Rv4(i)
      EXIT
    ENDIF
  ENDDO
  !
  IF ( x0>Rv5(k) ) x0 = Rv5(k)
  !     .......... NEXT BISECTION STEP ..........
  600  x1 = (xu+x0)*0.5E0
  s1 = ABS(xu) + ABS(x0) + ABS(Eps1)
  s2 = s1 + ABS(x0-xu)/2.0E0
  IF ( s2==s1 ) THEN
    !     .......... K-TH EIGENVALUE FOUND ..........
    Rv5(k) = x1
    k = k - 1
    IF ( k<m1 ) GOTO 800
    GOTO 500
  ENDIF
  700 CONTINUE
  DO
    !     .......... IN-LINE PROCEDURE FOR STURM SEQUENCE ..........
    s = p - 1
    u = 1.0E0
    !
    DO i = p, q
      IF ( u/=0.0E0 ) THEN
        v = E2(i)/u
      ELSE
        v = ABS(E(i))/machep
        IF ( E2(i)==0.0E0 ) v = 0.0E0
      ENDIF
      u = D(i) - x1 - v
      IF ( u<0.0E0 ) s = s + 1
    ENDDO
    !
    SELECT CASE (isturm)
      CASE (1)
        IF ( s<m1 ) THEN
          xu = x1
        ELSEIF ( s==m1 ) THEN
          xu = x1
          t1 = x1
          GOTO 200
        ELSE
          x0 = x1
        ENDIF
        GOTO 100
      CASE (2)
        IF ( s<m22 ) THEN
          xu = x1
        ELSEIF ( s==m22 ) THEN
          t2 = x1
          GOTO 300
        ELSE
          x0 = x1
        ENDIF
        GOTO 100
      CASE (3)
        m1 = s + 1
        x1 = Ub
        isturm = 4
      CASE (4)
        m2 = s
        IF ( m1>m2 ) GOTO 900
        !     .......... FIND ROOTS BY BISECTION ..........
        x0 = Ub
        isturm = 5
        !
        DO i = m1, m2
          Rv5(i) = Ub
          Rv4(i) = Lb
        ENDDO
        !     .......... LOOP FOR K-TH EIGENVALUE
        !                FOR K=M2 STEP -1 UNTIL M1 DO --
        !                (-DO- NOT USED TO LEGALIZE -COMPUTED GO TO-) ..........
        k = m2
        GOTO 500
      CASE DEFAULT
        !     .......... REFINE INTERVALS ..........
        IF ( s>=k ) THEN
          x0 = x1
        ELSE
          xu = x1
          IF ( s>=m1 ) THEN
            Rv4(s+1) = x1
            IF ( Rv5(s)>x1 ) Rv5(s) = x1
          ELSE
            Rv4(m1) = x1
          ENDIF
        ENDIF
        GOTO 600
    END SELECT
  ENDDO
  !     .......... ORDER EIGENVALUES TAGGED WITH THEIR
  !                SUBMATRIX ASSOCIATIONS ..........
  800  s = r
  r = r + m2 - m1 + 1
  j = 1
  k = m1
  !
  DO l = 1, r
    IF ( j<=s ) THEN
      IF ( k>m2 ) EXIT
      IF ( Rv5(k)>=W(l) ) THEN
        j = j + 1
        CYCLE
      ELSE
        !
        DO ii = j, s
          i = l + s - ii
          W(i+1) = W(i)
          Ind(i+1) = Ind(i)
        ENDDO
      ENDIF
    ENDIF
    !
    W(l) = Rv5(k)
    Ind(l) = tag
    k = k + 1
  ENDDO
  !
  900 CONTINUE
  IF ( q<N ) GOTO 400
  1000 Lb = t1
  Ub = t2
END SUBROUTINE TRIDIB
