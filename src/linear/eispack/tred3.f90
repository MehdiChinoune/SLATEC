!** TRED3
SUBROUTINE TRED3(N,Nv,A,D,E,E2)
  IMPLICIT NONE
  !>
  !***
  !  Reduce a real symmetric matrix stored in packed form to
  !            symmetric tridiagonal matrix using orthogonal
  !            transformations.
  !***
  ! **Library:**   SLATEC (EISPACK)
  !***
  ! **Category:**  D4C1B1
  !***
  ! **Type:**      SINGLE PRECISION (TRED3-S)
  !***
  ! **Keywords:**  EIGENVALUES, EIGENVECTORS, EISPACK
  !***
  ! **Author:**  Smith, B. T., et al.
  !***
  ! **Description:**
  !
  !     This subroutine is a translation of the ALGOL procedure TRED3,
  !     NUM. MATH. 11, 181-195(1968) by Martin, Reinsch, and Wilkinson.
  !     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
  !
  !     This subroutine reduces a REAL SYMMETRIC matrix, stored as
  !     a one-dimensional array, to a symmetric tridiagonal matrix
  !     using orthogonal similarity transformations.
  !
  !     On Input
  !
  !        N is the order of the matrix A.  N is an INTEGER variable.
  !
  !        NV is an INTEGER variable set equal to the dimension of the
  !          array A as specified in the calling program.  NV must not
  !          be less than  N*(N+1)/2.
  !
  !        A contains the lower triangle, stored row-wise, of the real
  !          symmetric packed matrix.  A is a one-dimensional REAL
  !          array, dimensioned A(NV).
  !
  !     On Output
  !
  !        A contains information about the orthogonal transformations
  !          used in the reduction in its first N*(N+1)/2 positions.
  !
  !        D contains the diagonal elements of the symmetric tridiagonal
  !          matrix.  D is a one-dimensional REAL array, dimensioned D(N).
  !
  !        E contains the subdiagonal elements of the symmetric
  !          tridiagonal matrix in its last N-1 positions.  E(1) is set
  !          to zero.  E is a one-dimensional REAL array, dimensioned
  !          E(N).
  !
  !        E2 contains the squares of the corresponding elements of E.
  !          E2 may coincide with E if the squares are not needed.
  !          E2 is a one-dimensional REAL array, dimensioned E2(N).
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
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   760101  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)

  !
  INTEGER i, j, k, l, N, ii, iz, jk, Nv
  REAL A(*), D(*), E(*), E2(*)
  REAL f, g, h, hh, scalee
  !
  !     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
  !* FIRST EXECUTABLE STATEMENT  TRED3
  DO ii = 1, N
    i = N + 1 - ii
    l = i - 1
    iz = (i*l)/2
    h = 0.0E0
    scalee = 0.0E0
    IF ( l>=1 ) THEN
      !     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
      DO k = 1, l
        iz = iz + 1
        D(k) = A(iz)
        scalee = scalee + ABS(D(k))
      END DO
      !
      IF ( scalee/=0.0E0 ) THEN
        !
        DO k = 1, l
          D(k) = D(k)/scalee
          h = h + D(k)*D(k)
        END DO
        !
        E2(i) = scalee*scalee*h
        f = D(l)
        g = -SIGN(SQRT(h),f)
        E(i) = scalee*g
        h = h - f*g
        D(l) = f - g
        A(iz) = scalee*D(l)
        IF ( l/=1 ) THEN
          f = 0.0E0
          !
          DO j = 1, l
            g = 0.0E0
            jk = (j*(j-1))/2
            !     .......... FORM ELEMENT OF A*U ..........
            DO k = 1, l
              jk = jk + 1
              IF ( k>j ) jk = jk + k - 2
              g = g + A(jk)*D(k)
            END DO
            !     .......... FORM ELEMENT OF P ..........
            E(j) = g/h
            f = f + E(j)*D(j)
          END DO
          !
          hh = f/(h+h)
          jk = 0
          !     .......... FORM REDUCED A ..........
          DO j = 1, l
            f = D(j)
            g = E(j) - hh*f
            E(j) = g
            !
            DO k = 1, j
              jk = jk + 1
              A(jk) = A(jk) - f*E(k) - g*D(k)
            END DO
          END DO
        END IF
        GOTO 50
      END IF
    END IF
    E(i) = 0.0E0
    E2(i) = 0.0E0
    !
    50  D(i) = A(iz+1)
    A(iz+1) = scalee*SQRT(h)
  END DO
  !
END SUBROUTINE TRED3
