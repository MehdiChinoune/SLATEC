!*==BISECT.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK BISECT
      SUBROUTINE BISECT(N,Eps1,D,E,E2,Lb,Ub,Mm,M,W,Ind,Ierr,Rv4,Rv5)
      IMPLICIT NONE
!*--BISECT5
!*** Start of declarations inserted by SPAG
      REAL R1MACH
!*** End of declarations inserted by SPAG
!***BEGIN PROLOGUE  BISECT
!***PURPOSE  Compute the eigenvalues of a symmetric tridiagonal matrix
!            in a given interval using Sturm sequencing.
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4A5, D4C2A
!***TYPE      SINGLE PRECISION (BISECT-S)
!***KEYWORDS  EIGENVALUES, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the bisection technique
!     in the ALGOL procedure TRISTURM by Peters and Wilkinson.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971).
!
!     This subroutine finds those eigenvalues of a TRIDIAGONAL
!     SYMMETRIC matrix which lie in a specified interval,
!     using bisection.
!
!     On INPUT
!
!        N is the order of the matrix.  N is an INTEGER variable.
!
!        EPS1 is an absolute error tolerance for the computed
!          eigenvalues.  If the input EPS1 is non-positive,
!          it is reset for each submatrix to a default value,
!          namely, minus the product of the relative machine
!          precision and the 1-norm of the submatrix.
!          EPS1 is a REAL variable.
!
!        D contains the diagonal elements of the input matrix.
!          D is a one-dimensional REAL array, dimensioned D(N).
!
!        E contains the subdiagonal elements of the input matrix
!          in its last N-1 positions.  E(1) is arbitrary.
!          E is a one-dimensional REAL array, dimensioned E(N).
!
!        E2 contains the squares of the corresponding elements of E.
!          E2(1) is arbitrary.  E2 is a one-dimensional REAL array,
!          dimensioned E2(N).
!
!        LB and UB define the interval to be searched for eigenvalues.
!          If LB is not less than UB, no eigenvalues will be found.
!          LB and UB are REAL variables.
!
!        MM should be set to an upper bound for the number of
!          eigenvalues in the interval.  WARNING - If more than
!          MM eigenvalues are determined to lie in the interval,
!          an error return is made with no eigenvalues found.
!          MM is an INTEGER variable.
!
!     On OUTPUT
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
!        M is the number of eigenvalues determined to lie in (LB,UB).
!          M is an INTEGER variable.
!
!        W contains the M eigenvalues in ascending order.
!          W is a one-dimensional REAL array, dimensioned W(MM).
!
!        IND contains in its first M positions the submatrix indices
!          associated with the corresponding eigenvalues in W --
!          1 for eigenvalues belonging to the first submatrix from
!          the top, 2 for those belonging to the second submatrix, etc.
!          IND is an one-dimensional INTEGER array, dimensioned IND(MM).
!
!        IERR is an INTEGER flag set to
!          Zero       for normal return,
!          3*N+1      if M exceeds MM.  In this case, M contains the
!                     number of eigenvalues determined to lie in
!                     (LB,UB).
!
!        RV4 and RV5 are one-dimensional REAL arrays used for temporary
!          storage, dimensioned RV4(N) and RV5(N).
!
!     The ALGOL procedure STURMCNT contained in TRISTURM
!     appears in BISECT in-line.
!
!     Note that subroutine TQL1 or IMTQL1 is generally faster than
!     BISECT, if more than N/4 eigenvalues are to be found.
!
!     Questions and comments should be directed to B. S. Garbow,
!     Applied Mathematics Division, ARGONNE NATIONAL LABORATORY
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
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  BISECT
!
      INTEGER i , j , k , l , M , N , p , q , r , s , ii , Mm , m1 , m2 , tag , 
     &        Ierr , isturm
      REAL D(*) , E(*) , E2(*) , W(*) , Rv4(*) , Rv5(*)
      REAL u , v , Lb , t1 , t2 , Ub , xu , x0 , x1 , Eps1 , machep , s1 , s2
      INTEGER Ind(*)
      LOGICAL first
!
      SAVE first , machep
      DATA first/.TRUE./
!***FIRST EXECUTABLE STATEMENT  BISECT
      IF ( first ) machep = R1MACH(4)
      first = .FALSE.
!
      Ierr = 0
      tag = 0
      t1 = Lb
      t2 = Ub
!     .......... LOOK FOR SMALL SUB-DIAGONAL ENTRIES ..........
      DO i = 1 , N
        IF ( i/=1 ) THEN
          s1 = ABS(D(i)) + ABS(D(i-1))
          s2 = s1 + ABS(E(i))
          IF ( s2>s1 ) CYCLE
        ENDIF
        E2(i) = 0.0E0
      ENDDO
!     .......... DETERMINE THE NUMBER OF EIGENVALUES
!                IN THE INTERVAL ..........
      p = 1
      q = N
      x1 = Ub
      isturm = 1
      GOTO 400
!     .......... ESTABLISH AND PROCESS NEXT SUBMATRIX, REFINING
!                INTERVAL BY THE GERSCHGORIN BOUNDS ..........
 100  IF ( r==M ) GOTO 700
      tag = tag + 1
      p = q + 1
      xu = D(p)
      x0 = D(p)
      u = 0.0E0
!
      DO q = p , N
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
        GOTO 400
      ELSE
!     .......... CHECK FOR ISOLATED ROOT WITHIN INTERVAL ..........
        IF ( t1>D(p).OR.D(p)>=t2 ) GOTO 600
        m1 = p
        m2 = p
        Rv5(p) = D(p)
        GOTO 500
      ENDIF
 200  xu = Lb
!     .......... FOR I=K STEP -1 UNTIL M1 DO -- ..........
      DO ii = m1 , k
        i = m1 + k - ii
        IF ( xu<Rv4(i) ) THEN
          xu = Rv4(i)
          EXIT
        ENDIF
      ENDDO
!
      IF ( x0>Rv5(k) ) x0 = Rv5(k)
!     .......... NEXT BISECTION STEP ..........
 300  x1 = (xu+x0)*0.5E0
      s1 = 2.0E0*(ABS(xu)+ABS(x0)+ABS(Eps1))
      s2 = s1 + ABS(x0-xu)
      IF ( s2==s1 ) THEN
!     .......... K-TH EIGENVALUE FOUND ..........
        Rv5(k) = x1
        k = k - 1
        IF ( k<m1 ) GOTO 500
        GOTO 200
      ENDIF
 400  DO
!     .......... IN-LINE PROCEDURE FOR STURM SEQUENCE ..........
        s = p - 1
        u = 1.0E0
!
        DO i = p , q
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
          M = s
          x1 = Lb
          isturm = 2
        CASE (2)
          M = M - s
          IF ( M>Mm ) THEN
!     .......... SET ERROR -- UNDERESTIMATE OF NUMBER OF
!                EIGENVALUES IN INTERVAL ..........
            Ierr = 3*N + 1
            GOTO 700
          ELSE
            q = 0
            r = 0
            GOTO 100
          ENDIF
        CASE (3)
          m1 = s + 1
          x1 = Ub
          isturm = 4
        CASE (4)
          m2 = s
          IF ( m1>m2 ) GOTO 600
!     .......... FIND ROOTS BY BISECTION ..........
          x0 = Ub
          isturm = 5
!
          DO i = m1 , m2
            Rv5(i) = Ub
            Rv4(i) = Lb
          ENDDO
!     .......... LOOP FOR K-TH EIGENVALUE
!                FOR K=M2 STEP -1 UNTIL M1 DO --
!                (-DO- NOT USED TO LEGALIZE -COMPUTED GO TO-) ..........
          k = m2
          GOTO 200
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
          GOTO 300
        END SELECT
      ENDDO
!     .......... ORDER EIGENVALUES TAGGED WITH THEIR
!                SUBMATRIX ASSOCIATIONS ..........
 500  s = r
      r = r + m2 - m1 + 1
      j = 1
      k = m1
!
      DO l = 1 , r
        IF ( j<=s ) THEN
          IF ( k>m2 ) EXIT
          IF ( Rv5(k)>=W(l) ) THEN
            j = j + 1
            CYCLE
          ELSE
!
            DO ii = j , s
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
 600  IF ( q<N ) GOTO 100
 700  Lb = t1
      Ub = t2
      END SUBROUTINE BISECT
