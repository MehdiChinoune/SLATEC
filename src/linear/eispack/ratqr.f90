!** RATQR
SUBROUTINE RATQR(N,Eps1,D,E,E2,M,W,Ind,Bd,Type,Idef,Ierr)
  !>
  !  Compute the largest or smallest eigenvalues of a symmetric
  !            tridiagonal matrix using the rational QR method with Newton
  !            correction.
  !***
  ! **Library:**   SLATEC (EISPACK)
  !***
  ! **Category:**  D4A5, D4C2A
  !***
  ! **Type:**      SINGLE PRECISION (RATQR-S)
  !***
  ! **Keywords:**  EIGENVALUES, EIGENVECTORS, EISPACK
  !***
  ! **Author:**  Smith, B. T., et al.
  !***
  ! **Description:**
  !
  !     This subroutine is a translation of the ALGOL procedure RATQR,
  !     NUM. MATH. 11, 264-272(1968) by REINSCH and BAUER.
  !     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 257-265(1971).
  !
  !     This subroutine finds the algebraically smallest or largest
  !     eigenvalues of a SYMMETRIC TRIDIAGONAL matrix by the
  !     rational QR method with Newton corrections.
  !
  !     On Input
  !
  !        N is the order of the matrix.  N is an INTEGER variable.
  !
  !        EPS1 is a theoretical absolute error tolerance for the
  !          computed eigenvalues.  If the input EPS1 is non-positive, or
  !          indeed smaller than its default value, it is reset at each
  !          iteration to the respective default value, namely, the
  !          product of the relative machine precision and the magnitude
  !          of the current eigenvalue iterate.  The theoretical absolute
  !          error in the K-th eigenvalue is usually not greater than
  !          K times EPS1.  EPS1 is a REAL variable.
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
  !        M is the number of eigenvalues to be found.  M is an INTEGER
  !          variable.
  !
  !        IDEF should be set to 1 if the input matrix is known to be
  !          positive definite, to -1 if the input matrix is known to
  !          be negative definite, and to 0 otherwise.  IDEF is an
  !          INTEGER variable.
  !
  !        TYPE should be set to .TRUE. if the smallest eigenvalues are
  !          to be found, and to .FALSE. if the largest eigenvalues are
  !          to be found.  TYPE is a LOGICAL variable.
  !
  !     On Output
  !
  !        EPS1 is unaltered unless it has been reset to its
  !          (last) default value.
  !
  !        D and E are unaltered (unless W overwrites D).
  !
  !        Elements of E2, corresponding to elements of E regarded as
  !          negligible, have been replaced by zero causing the matrix
  !          to split into a direct sum of submatrices.  E2(1) is set
  !          to 0.0e0 if the smallest eigenvalues have been found, and
  !          to 2.0e0 if the largest eigenvalues have been found.  E2
  !          is otherwise unaltered (unless overwritten by BD).
  !
  !        W contains the M algebraically smallest eigenvalues in
  !          ascending order, or the M largest eigenvalues in descending
  !          order.  If an error exit is made because of an incorrect
  !          specification of IDEF, no eigenvalues are found.  If the
  !          Newton iterates for a particular eigenvalue are not monotone,
  !          the best estimate obtained is returned and IERR is set.
  !          W is a one-dimensional REAL array, dimensioned W(N).  W need
  !          not be distinct from D.
  !
  !        IND contains in its first M positions the submatrix indices
  !          associated with the corresponding eigenvalues in W --
  !          1 for eigenvalues belonging to the first submatrix from
  !          the top, 2 for those belonging to the second submatrix, etc.
  !          IND is an one-dimensional INTEGER array, dimensioned IND(N).
  !
  !        BD contains refined bounds for the theoretical errors of the
  !          corresponding eigenvalues in W.  These bounds are usually
  !          within the tolerance specified by EPS1.  BD is a one-
  !          dimensional REAL array, dimensioned BD(N).  BD need not be
  !          distinct from E2.
  !
  !        IERR is an INTEGER flag set to
  !          Zero       for normal return,
  !          6*N+1      if  IDEF  is set to 1 and  TYPE  to .TRUE.
  !                     when the matrix is NOT positive definite, or
  !                     if  IDEF  is set to -1 and  TYPE  to .FALSE.
  !                     when the matrix is NOT negative definite,
  !                     no eigenvalues are computed, or
  !                     M is greater than N,
  !          5*N+K      if successive iterates to the K-th eigenvalue
  !                     are NOT monotone increasing, where K refers
  !                     to the last such occurrence.
  !
  !     Note that subroutine TRIDIB is generally faster and more
  !     accurate than RATQR if the eigenvalues are clustered.
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
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  USE service, ONLY : R1MACH
  INTEGER i, j, k, M, N, ii, jj, k1, Idef, Ierr, jdef
  REAL(SP) D(*), E(*), E2(*), W(*), Bd(*)
  REAL(SP) f, p, q, r, s, ep, qp, err, tot, Eps1, delta
  INTEGER Ind(*)
  LOGICAL :: Type
  !
  REAL(SP), PARAMETER :: machep = R1MACH(4)
  !* FIRST EXECUTABLE STATEMENT  RATQR
  !
  Ierr = 0
  jdef = Idef
  !     .......... COPY D ARRAY INTO W ..........
  DO i = 1, N
    W(i) = D(i)
  END DO
  !
  IF ( .NOT.(Type) ) THEN
    j = 1
    GOTO 200
  END IF
  100  err = 0.0E0
  s = 0.0E0
  !     .......... LOOK FOR SMALL SUB-DIAGONAL ENTRIES AND DEFINE
  !                INITIAL SHIFT FROM LOWER GERSCHGORIN BOUND.
  !                COPY E2 ARRAY INTO BD ..........
  tot = W(1)
  q = 0.0E0
  j = 0
  !
  DO i = 1, N
    p = q
    IF ( i==1 ) THEN
      E2(i) = 0.0E0
    ELSEIF ( p<=machep*(ABS(D(i))+ABS(D(i-1))) ) THEN
      E2(i) = 0.0E0
    END IF
    Bd(i) = E2(i)
    !     .......... COUNT ALSO IF ELEMENT OF E2 HAS UNDERFLOWED ..........
    IF ( E2(i)==0.0E0 ) j = j + 1
    Ind(i) = j
    q = 0.0E0
    IF ( i/=N ) q = ABS(E(i+1))
    tot = MIN(W(i)-p-q,tot)
  END DO
  !
  IF ( jdef==1.AND.tot<0.0E0 ) THEN
    tot = 0.0E0
  ELSE
    !
    DO i = 1, N
      W(i) = W(i) - tot
      !
    END DO
  END IF
  !
  DO k = 1, M
    DO
      !     .......... NEXT QR TRANSFORMATION ..........
      tot = tot + s
      delta = W(N) - s
      i = N
      f = ABS(machep*tot)
      IF ( Eps1<f ) Eps1 = f
      IF ( delta>Eps1 ) THEN
        !     .......... REPLACE SMALL SUB-DIAGONAL SQUARES BY ZERO
        !                TO REDUCE THE INCIDENCE OF UNDERFLOWS ..........
        IF ( k/=N ) THEN
          k1 = k + 1
          DO j = k1, N
            IF ( Bd(j)<=(machep*(W(j)+W(j-1)))**2 ) Bd(j) = 0.0E0
          END DO
        END IF
        !
        f = Bd(N)/delta
        qp = delta + f
        p = 1.0E0
        IF ( k/=N ) THEN
          k1 = N - k
          !     .......... FOR I=N-1 STEP -1 UNTIL K DO -- ..........
          DO ii = 1, k1
            i = N - ii
            q = W(i) - s - f
            r = q/qp
            p = p*r + 1.0E0
            ep = f*r
            W(i+1) = qp + ep
            delta = q - ep
            IF ( delta>Eps1 ) THEN
              f = Bd(i)/q
              qp = delta + f
              Bd(i+1) = qp*ep
            ELSE
              IF ( delta>=(-Eps1) ) GOTO 150
              GOTO 300
            END IF
          END DO
        END IF
        !
        W(k) = qp
        s = qp/p
        IF ( tot+s<=tot ) THEN
          !     .......... SET ERROR -- IRREGULAR END OF ITERATION.
          !                DEFLATE MINIMUM DIAGONAL ELEMENT ..........
          Ierr = 5*N + k
          s = 0.0E0
          delta = qp
          !
          DO j = k, N
            IF ( W(j)<=delta ) THEN
              i = j
              delta = W(j)
            END IF
          END DO
          EXIT
        END IF
      ELSE
        IF ( delta>=(-Eps1) ) EXIT
        GOTO 300
      END IF
    END DO
    !     .......... CONVERGENCE ..........
    150  IF ( i<N ) Bd(i+1) = Bd(i)*f/qp
    ii = Ind(i)
    IF ( i/=k ) THEN
      k1 = i - k
      !     .......... FOR J=I-1 STEP -1 UNTIL K DO -- ..........
      DO jj = 1, k1
        j = i - jj
        W(j+1) = W(j) - s
        Bd(j+1) = Bd(j)
        Ind(j+1) = Ind(j)
      END DO
    END IF
    !
    W(k) = tot
    err = err + ABS(delta)
    Bd(k) = err
    Ind(k) = ii
  END DO
  !
  IF ( Type ) RETURN
  f = Bd(1)
  E2(1) = 2.0E0
  Bd(1) = f
  j = 2
  !     .......... NEGATE ELEMENTS OF W FOR LARGEST VALUES ..........
  200 CONTINUE
  DO i = 1, N
    W(i) = -W(i)
  END DO
  !
  jdef = -jdef
  IF ( j==1 ) GOTO 100
  IF ( j==2 ) RETURN
  !     .......... SET ERROR -- IDEF SPECIFIED INCORRECTLY ..........
  300  Ierr = 6*N + 1
  RETURN
END SUBROUTINE RATQR
