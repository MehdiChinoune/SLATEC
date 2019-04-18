!** ORTHES
SUBROUTINE ORTHES(Nm,N,Low,Igh,A,Ort)
  !>
  !  Reduce a real general matrix to upper Hessenberg form
  !            using orthogonal similarity transformations.
  !***
  ! **Library:**   SLATEC (EISPACK)
  !***
  ! **Category:**  D4C1B2
  !***
  ! **Type:**      SINGLE PRECISION (ORTHES-S, CORTH-C)
  !***
  ! **Keywords:**  EIGENVALUES, EIGENVECTORS, EISPACK
  !***
  ! **Author:**  Smith, B. T., et al.
  !***
  ! **Description:**
  !
  !     This subroutine is a translation of the ALGOL procedure ORTHES,
  !     NUM. MATH. 12, 349-368(1968) by Martin and Wilkinson.
  !     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
  !
  !     Given a REAL GENERAL matrix, this subroutine
  !     reduces a submatrix situated in rows and columns
  !     LOW through IGH to upper Hessenberg form by
  !     orthogonal similarity transformations.
  !
  !     On INPUT
  !
  !        NM must be set to the row dimension of the two-dimensional
  !          array parameter, A, as declared in the calling program
  !          dimension statement.  NM is an INTEGER variable.
  !
  !        N is the order of the matrix A.  N is an INTEGER variable.
  !          N must be less than or equal to NM.
  !
  !        LOW and IGH are two INTEGER variables determined by the
  !          balancing subroutine  BALANC.  If  BALANC  has not been
  !          used, set LOW=1 and IGH equal to the order of the matrix, N.
  !
  !        A contains the general matrix to be reduced to upper
  !          Hessenberg form.  A is a two-dimensional REAL array,
  !          dimensioned A(NM,N).
  !
  !     On OUTPUT
  !
  !        A contains the upper Hessenberg matrix.  Some information about
  !          the orthogonal transformations used in the reduction
  !          is stored in the remaining triangle under the Hessenberg
  !          matrix.
  !
  !        ORT contains further information about the orthogonal trans-
  !          formations used in the reduction.  Only elements LOW+1
  !          through IGH are used.  ORT is a one-dimensional REAL array,
  !          dimensioned ORT(IGH).
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
  INTEGER i, j, m, N, ii, jj, la, mp, Nm, Igh, kp1, Low
  REAL A(Nm,*), Ort(*)
  REAL f, g, h, scalee
  !
  !* FIRST EXECUTABLE STATEMENT  ORTHES
  la = Igh - 1
  kp1 = Low + 1
  IF ( la>=kp1 ) THEN
    !
    DO m = kp1, la
      h = 0.0E0
      Ort(m) = 0.0E0
      scalee = 0.0E0
      !     .......... SCALE COLUMN (ALGOL TOL THEN NOT NEEDED) ..........
      DO i = m, Igh
        scalee = scalee + ABS(A(i,m-1))
      END DO
      !
      IF ( scalee/=0.0E0 ) THEN
        mp = m + Igh
        !     .......... FOR I=IGH STEP -1 UNTIL M DO -- ..........
        DO ii = m, Igh
          i = mp - ii
          Ort(i) = A(i,m-1)/scalee
          h = h + Ort(i)*Ort(i)
        END DO
        !
        g = -SIGN(SQRT(h),Ort(m))
        h = h - Ort(m)*g
        Ort(m) = Ort(m) - g
        !     .......... FORM (I-(U*UT)/H) * A ..........
        DO j = m, N
          f = 0.0E0
          !     .......... FOR I=IGH STEP -1 UNTIL M DO -- ..........
          DO ii = m, Igh
            i = mp - ii
            f = f + Ort(i)*A(i,j)
          END DO
          !
          f = f/h
          !
          DO i = m, Igh
            A(i,j) = A(i,j) - f*Ort(i)
          END DO
          !
        END DO
        !     .......... FORM (I-(U*UT)/H)*A*(I-(U*UT)/H) ..........
        DO i = 1, Igh
          f = 0.0E0
          !     .......... FOR J=IGH STEP -1 UNTIL M DO -- ..........
          DO jj = m, Igh
            j = mp - jj
            f = f + Ort(j)*A(i,j)
          END DO
          !
          f = f/h
          !
          DO j = m, Igh
            A(i,j) = A(i,j) - f*Ort(j)
          END DO
          !
        END DO
        !
        Ort(m) = scalee*Ort(m)
        A(m,m-1) = scalee*g
      END IF
    END DO
  END IF
  !
END SUBROUTINE ORTHES
