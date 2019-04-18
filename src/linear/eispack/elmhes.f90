!** ELMHES
SUBROUTINE ELMHES(Nm,N,Low,Igh,A,Intt)
  !>
  !  Reduce a real general matrix to upper Hessenberg form
  !            using stabilized elementary similarity transformations.
  !***
  ! **Library:**   SLATEC (EISPACK)
  !***
  ! **Category:**  D4C1B2
  !***
  ! **Type:**      SINGLE PRECISION (ELMHES-S, COMHES-C)
  !***
  ! **Keywords:**  EIGENVALUES, EIGENVECTORS, EISPACK
  !***
  ! **Author:**  Smith, B. T., et al.
  !***
  ! **Description:**
  !
  !     This subroutine is a translation of the ALGOL procedure ELMHES,
  !     NUM. MATH. 12, 349-368(1968) by Martin and Wilkinson.
  !     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
  !
  !     Given a REAL GENERAL matrix, this subroutine
  !     reduces a submatrix situated in rows and columns
  !     LOW through IGH to upper Hessenberg form by
  !     stabilized elementary similarity transformations.
  !
  !     On INPUT
  !
  !        NM must be set to the row dimension of the two-dimensional
  !          array parameter, A, as declared in the calling program
  !          dimension statement.  NM is an INTEGER variable.
  !
  !        N is the order of the matrix, A.  N is an INTEGER variable.
  !          N must be less than or equal to NM.
  !
  !        LOW and IGH are two INTEGER variables determined by the
  !          balancing subroutine  BALANC.  If  BALANC  has not been
  !          used, set LOW=1 and IGH equal to the order of the matrix, N.
  !
  !        A contains the input matrix.  A is a two-dimensional REAL
  !          array, dimensioned A(NM,N).
  !
  !     On OUTPUT
  !
  !        A contains the upper Hessenberg matrix.  The multipliers which
  !          were used in the reduction are stored in the remaining
  !          triangle under the Hessenberg matrix.
  !
  !        INT contains information on the rows and columns interchanged
  !          in the reduction.  Only elements LOW through IGH are used.
  !          INT is a one-dimensional INTEGER array, dimensioned INT(IGH).
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
  INTEGER i, j, m, N, la, Nm, Igh, kp1, Low, mm1, mp1
  REAL A(Nm,*)
  REAL x, y
  INTEGER Intt(*)
  !
  !* FIRST EXECUTABLE STATEMENT  ELMHES
  la = Igh - 1
  kp1 = Low + 1
  IF ( la>=kp1 ) THEN
    !
    DO m = kp1, la
      mm1 = m - 1
      x = 0.0E0
      i = m
      !
      DO j = m, Igh
        IF ( ABS(A(j,mm1))>ABS(x) ) THEN
          x = A(j,mm1)
          i = j
        END IF
      END DO
      !
      Intt(m) = i
      IF ( i/=m ) THEN
        !    .......... INTERCHANGE ROWS AND COLUMNS OF A ..........
        DO j = mm1, N
          y = A(i,j)
          A(i,j) = A(m,j)
          A(m,j) = y
        END DO
        !
        DO j = 1, Igh
          y = A(j,i)
          A(j,i) = A(j,m)
          A(j,m) = y
        END DO
      END IF
      !    .......... END INTERCHANGE ..........
      IF ( x/=0.0E0 ) THEN
        mp1 = m + 1
        !
        DO i = mp1, Igh
          y = A(i,mm1)
          IF ( y/=0.0E0 ) THEN
            y = y/x
            A(i,mm1) = y
            !
            DO j = m, N
              A(i,j) = A(i,j) - y*A(m,j)
            END DO
            !
            DO j = 1, Igh
              A(j,m) = A(j,m) + y*A(j,i)
            END DO
          END IF
          !
        END DO
      END IF
      !
    END DO
  END IF
  !
END SUBROUTINE ELMHES
