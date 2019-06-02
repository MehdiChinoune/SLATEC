!** TSTURM
SUBROUTINE TSTURM(Nm,N,Eps1,D,E,E2,Lb,Ub,Mm,M,W,Z,Ierr,Rv1,Rv2,Rv3,Rv4,Rv5,Rv6)
  !>
  !  Find those eigenvalues of a symmetric tridiagonal matrix
  !            in a given interval and their associated eigenvectors by
  !            Sturm sequencing.
  !***
  ! **Library:**   SLATEC (EISPACK)
  !***
  ! **Category:**  D4A5, D4C2A
  !***
  ! **Type:**      SINGLE PRECISION (TSTURM-S)
  !***
  ! **Keywords:**  EIGENVALUES, EIGENVECTORS, EISPACK
  !***
  ! **Author:**  Smith, B. T., et al.
  !***
  ! **Description:**
  !
  !     This subroutine finds those eigenvalues of a TRIDIAGONAL
  !     SYMMETRIC matrix which lie in a specified interval and their
  !     associated eigenvectors, using bisection and inverse iteration.
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
  !        EPS1 is an absolute error tolerance for the computed eigen-
  !          values.  It should be chosen so that the accuracy of these
  !          eigenvalues is commensurate with relative perturbations of
  !          the order of the relative machine precision in the matrix
  !          elements.  If the input EPS1 is non-positive, it is reset
  !          for each submatrix to a default value, namely, minus the
  !          product of the relative machine precision and the 1-norm of
  !          the submatrix.  EPS1 is a REAL variable.
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
  !        LB and UB define the interval to be searched for eigenvalues.
  !          If LB is not less than UB, no eigenvalues will be found.
  !          LB and UB are REAL variables.
  !
  !        MM should be set to an upper bound for the number of
  !          eigenvalues in the interval.  MM is an INTEGER variable.
  !          WARNING -  If more than MM eigenvalues are determined to lie
  !          in the interval, an error return is made with no values or
  !          vectors found.
  !
  !     On Output
  !
  !        EPS1 is unaltered unless it has been reset to its
  !          (last) default value.
  !
  !        D and E are unaltered.
  !
  !        Elements of E2, corresponding to elements of E regarded as
  !          negligible, have been replaced by zero causing the matrix to
  !          split into a direct sum of submatrices.  E2(1) is also set
  !          to zero.
  !
  !        M is the number of eigenvalues determined to lie in (LB,UB).
  !          M is an INTEGER variable.
  !
  !        W contains the M eigenvalues in ascending order if the matrix
  !          does not split.  If the matrix splits, the eigenvalues are
  !          in ascending order for each submatrix.  If a vector error
  !          exit is made, W contains those values already found.  W is a
  !          one-dimensional REAL array, dimensioned W(MM).
  !
  !        Z contains the associated set of orthonormal eigenvectors.
  !          If an error exit is made, Z contains those vectors already
  !          found.  Z is a one-dimensional REAL array, dimensioned
  !          Z(NM,MM).
  !
  !        IERR is an INTEGER flag set to
  !          Zero       for normal return,
  !          3*N+1      if M exceeds MM no eigenvalues or eigenvectors
  !                     are computed,
  !          4*N+J      if the eigenvector corresponding to the J-th
  !                     eigenvalue fails to converge in 5 iterations, then
  !                     the eigenvalues and eigenvectors in W and Z should
  !                     be correct for indices 1, 2, ..., J-1.
  !
  !        RV1, RV2, RV3, RV4, RV5, and RV6 are temporary storage arrays,
  !          dimensioned RV1(N), RV2(N), RV3(N), RV4(N), RV5(N), and
  !          RV6(N).
  !
  !     The ALGOL procedure STURMCNT contained in TRISTURM
  !     appears in TSTURM in-line.
  !
  !     Questions and comments should be directed to B. S. Garbow,
  !     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
  !     ------------------------------------------------------------------
  !
  !***
  ! **References:**  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
  !                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
  !                 system Routines - EISPACK Guide, Springer-Verlag,
  !                 1976.
  !***
  ! **Routines called:**  R1MACH

  !* REVISION HISTORY  (YYMMDD)
  !   760101  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  USE service, ONLY : R1MACH
  INTEGER :: i, j, k, M, N, p, q, r, s, ii, ip, jj, Mm, m1, m2, Nm, its, Ierr, &
    group, isturm
  REAL :: D(*), E(*), E2(*), W(*), Z(Nm,*), Rv1(*), Rv2(*), Rv3(*), Rv4(*), &
    Rv5(*), Rv6(*)
  REAL u, v, Lb, t1, t2, Ub, uk, xu, x0, x1, Eps1, eps2, eps3, eps4, norm, s1, s2
  !
  REAL, PARAMETER :: machep = R1MACH(4)
  !* FIRST EXECUTABLE STATEMENT  TSTURM
  Ierr = 0
  t1 = Lb
  t2 = Ub
  !     .......... LOOK FOR SMALL SUB-DIAGONAL ENTRIES ..........
  DO i = 1, N
    IF ( i/=1 ) THEN
      s1 = ABS(D(i)) + ABS(D(i-1))
      s2 = s1 + ABS(E(i))
      IF ( s2>s1 ) CYCLE
    END IF
    E2(i) = 0.0E0
  END DO
  !     .......... DETERMINE THE NUMBER OF EIGENVALUES
  !                IN THE INTERVAL ..........
  p = 1
  q = N
  x1 = Ub
  isturm = 1
  GOTO 400
  !     .......... ESTABLISH AND PROCESS NEXT SUBMATRIX, REFINING
  !                INTERVAL BY THE GERSCHGORIN BOUNDS ..........
  100 CONTINUE
  IF ( r==M ) GOTO 700
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
    END IF
    xu = MIN(D(q)-(x1+u),xu)
    x0 = MAX(D(q)+(x1+u),x0)
    IF ( v==0.0E0 ) EXIT
  END DO
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
    IF ( t1<=D(p).AND.D(p)<t2 ) THEN
      r = r + 1
      !
      DO i = 1, N
        Z(i,r) = 0.0E0
      END DO
      !
      W(r) = D(p)
      Z(p,r) = 1.0E0
    END IF
    GOTO 500
  END IF
  200  xu = Lb
  !     .......... FOR I=K STEP -1 UNTIL M1 DO -- ..........
  DO ii = m1, k
    i = m1 + k - ii
    IF ( xu<Rv4(i) ) THEN
      xu = Rv4(i)
      EXIT
    END IF
  END DO
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
    IF ( k>=m1 ) GOTO 200
    !     .......... FIND VECTORS BY INVERSE ITERATION ..........
    norm = ABS(D(p))
    ip = p + 1
    !
    DO i = ip, q
      norm = MAX(norm,ABS(D(i))+ABS(E(i)))
    END DO
    !     .......... EPS2 IS THE CRITERION FOR GROUPING,
    !                EPS3 REPLACES ZERO PIVOTS AND EQUAL
    !                ROOTS ARE MODIFIED BY EPS3,
    !                EPS4 IS TAKEN VERY SMALL TO AVOID OVERFLOW ..........
    eps2 = 1.0E-3*norm
    uk = SQRT(REAL(q-p+5))
    eps3 = uk*machep*norm
    eps4 = uk*eps3
    uk = eps4/SQRT(uk)
    group = 0
    s = p
    !
    DO k = m1, m2
      r = r + 1
      its = 1
      W(r) = Rv5(k)
      x1 = Rv5(k)
      !     .......... LOOK FOR CLOSE OR COINCIDENT ROOTS ..........
      IF ( k/=m1 ) THEN
        IF ( x1-x0>=eps2 ) group = -1
        group = group + 1
        IF ( x1<=x0 ) x1 = x0 + eps3
      END IF
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
            xu = u/E(i)
            Rv4(i) = xu
            Rv1(i-1) = E(i)
            Rv2(i-1) = D(i) - x1
            Rv3(i-1) = 0.0E0
            IF ( i/=q ) Rv3(i-1) = E(i+1)
            u = v - xu*Rv2(i-1)
            v = -xu*Rv3(i-1)
            CYCLE
          END IF
        END IF
        u = D(i) - x1 - xu*v
        IF ( i/=q ) v = E(i+1)
      END DO
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
        END DO
        !     .......... ORTHOGONALIZE WITH RESPECT TO PREVIOUS
        !                MEMBERS OF GROUP ..........
        IF ( group/=0 ) THEN
          !
          DO jj = 1, group
            j = r - group - 1 + jj
            xu = 0.0E0
            !
            DO i = p, q
              xu = xu + Rv6(i)*Z(i,j)
            END DO
            !
            DO i = p, q
              Rv6(i) = Rv6(i) - xu*Z(i,j)
            END DO
            !
          END DO
        END IF
        !
        norm = 0.0E0
        !
        DO i = p, q
          norm = norm + ABS(Rv6(i))
        END DO
        !
        IF ( norm>=1.0E0 ) THEN
          !     .......... NORMALIZE SO THAT SUM OF SQUARES IS
          !                1 AND EXPAND TO FULL ORDER ..........
          u = 0.0E0
          !
          DO i = p, q
            u = u + Rv6(i)**2
          END DO
          !
          xu = 1.0E0/SQRT(u)
          !
          DO i = 1, N
            Z(i,r) = 0.0E0
          END DO
          !
          DO i = p, q
            Z(i,r) = Rv6(i)*xu
          END DO
          !
          x0 = x1
          EXIT
        ELSE
          !     .......... FORWARD SUBSTITUTION ..........
          IF ( its==5 ) GOTO 600
          IF ( norm/=0.0E0 ) THEN
            xu = eps4/norm
            !
            DO i = p, q
              Rv6(i) = Rv6(i)*xu
            END DO
          ELSE
            Rv6(s) = eps4
            s = s + 1
            IF ( s>q ) s = p
          END IF
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
            END IF
            Rv6(i) = u - Rv4(i)*Rv6(i-1)
          END DO
          !
          its = its + 1
        END IF
      END DO
    END DO
    GOTO 500
  END IF
  400 CONTINUE
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
      END IF
      u = D(i) - x1 - v
      IF ( u<0.0E0 ) s = s + 1
    END DO
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
        END IF
      CASE (3)
        m1 = s + 1
        x1 = Ub
        isturm = 4
      CASE (4)
        m2 = s
        IF ( m1>m2 ) EXIT
        !     .......... FIND ROOTS BY BISECTION ..........
        x0 = Ub
        isturm = 5
        !
        DO i = m1, m2
          Rv5(i) = Ub
          Rv4(i) = Lb
        END DO
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
          END IF
        END IF
        GOTO 300
    END SELECT
  END DO
  !
  500 CONTINUE
  IF ( q>=N ) GOTO 700
  GOTO 100
  !     .......... SET ERROR -- NON-CONVERGED EIGENVECTOR ..........
  600  Ierr = 4*N + r
  700  Lb = t1
  Ub = t2
END SUBROUTINE TSTURM
