!DECK COMHES
SUBROUTINE COMHES(Nm,N,Low,Igh,Ar,Ai,Int)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  COMHES
  !***PURPOSE  Reduce a complex general matrix to complex upper Hessenberg
  !            form using stabilized elementary similarity
  !            transformations.
  !***LIBRARY   SLATEC (EISPACK)
  !***CATEGORY  D4C1B2
  !***TYPE      COMPLEX (ELMHES-S, COMHES-C)
  !***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
  !***AUTHOR  Smith, B. T., et al.
  !***DESCRIPTION
  !
  !     This subroutine is a translation of the ALGOL procedure COMHES,
  !     NUM. MATH. 12, 349-368(1968) by Martin and Wilkinson.
  !     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
  !
  !     Given a COMPLEX GENERAL matrix, this subroutine
  !     reduces a submatrix situated in rows and columns
  !     LOW through IGH to upper Hessenberg form by
  !     stabilized elementary similarity transformations.
  !
  !     On INPUT
  !
  !        NM must be set to the row dimension of the two-dimensional
  !          array parameters, AR and AI, as declared in the calling
  !          program dimension statement.  NM is an INTEGER variable.
  !
  !        N is the order of the matrix A=(AR,AI).  N is an INTEGER
  !          variable.  N must be less than or equal to NM.
  !
  !        LOW and IGH are two INTEGER variables determined by the
  !          balancing subroutine  CBAL.  If  CBAL  has not been used,
  !          set LOW=1 and IGH equal to the order of the matrix, N.
  !
  !        AR and AI contain the real and imaginary parts, respectively,
  !          of the complex input matrix.  AR and AI are two-dimensional
  !          REAL arrays, dimensioned AR(NM,N) and AI(NM,N).
  !
  !     On OUTPUT
  !
  !        AR and AI contain the real and imaginary parts, respectively,
  !          of the upper Hessenberg matrix.  The multipliers which
  !          were used in the reduction are stored in the remaining
  !          triangles under the Hessenberg matrix.
  !
  !        INT contains information on the rows and columns
  !          interchanged in the reduction.  Only elements LOW through
  !          IGH are used.  INT is a one-dimensional INTEGER array,
  !          dimensioned INT(IGH).
  !
  !     Calls CDIV for complex division.
  !
  !     Questions and comments should be directed to B. S. Garbow,
  !     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
  !     ------------------------------------------------------------------
  !
  !***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
  !                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
  !                 system Routines - EISPACK Guide, Springer-Verlag,
  !                 1976.
  !***ROUTINES CALLED  CDIV
  !***REVISION HISTORY  (YYMMDD)
  !   760101  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  COMHES
  !
  INTEGER i, j, m, N, la, Nm, Igh, kp1, Low, mm1, mp1
  REAL Ar(Nm,*), Ai(Nm,*)
  REAL xr, xi, yr, yi
  INTEGER Int(*)
  !
  !***FIRST EXECUTABLE STATEMENT  COMHES
  la = Igh - 1
  kp1 = Low + 1
  IF ( la>=kp1 ) THEN
    !
    DO m = kp1, la
      mm1 = m - 1
      xr = 0.0E0
      xi = 0.0E0
      i = m
      !
      DO j = m, Igh
        IF ( ABS(Ar(j,mm1))+ABS(Ai(j,mm1))>ABS(xr)+ABS(xi) ) THEN
          xr = Ar(j,mm1)
          xi = Ai(j,mm1)
          i = j
        ENDIF
      ENDDO
      !
      Int(m) = i
      IF ( i/=m ) THEN
        !     .......... INTERCHANGE ROWS AND COLUMNS OF AR AND AI ..........
        DO j = mm1, N
          yr = Ar(i,j)
          Ar(i,j) = Ar(m,j)
          Ar(m,j) = yr
          yi = Ai(i,j)
          Ai(i,j) = Ai(m,j)
          Ai(m,j) = yi
        ENDDO
        !
        DO j = 1, Igh
          yr = Ar(j,i)
          Ar(j,i) = Ar(j,m)
          Ar(j,m) = yr
          yi = Ai(j,i)
          Ai(j,i) = Ai(j,m)
          Ai(j,m) = yi
        ENDDO
      ENDIF
      !     .......... END INTERCHANGE ..........
      IF ( xr/=0.0E0.OR.xi/=0.0E0 ) THEN
        mp1 = m + 1
        !
        DO i = mp1, Igh
          yr = Ar(i,mm1)
          yi = Ai(i,mm1)
          IF ( yr/=0.0E0.OR.yi/=0.0E0 ) THEN
            CALL CDIV(yr,yi,xr,xi,yr,yi)
            Ar(i,mm1) = yr
            Ai(i,mm1) = yi
            !
            DO j = m, N
              Ar(i,j) = Ar(i,j) - yr*Ar(m,j) + yi*Ai(m,j)
              Ai(i,j) = Ai(i,j) - yr*Ai(m,j) - yi*Ar(m,j)
            ENDDO
            !
            DO j = 1, Igh
              Ar(j,m) = Ar(j,m) + yr*Ar(j,i) - yi*Ai(j,i)
              Ai(j,m) = Ai(j,m) + yr*Ai(j,i) + yi*Ar(j,i)
            ENDDO
          ENDIF
          !
        ENDDO
      ENDIF
      !
    ENDDO
  ENDIF
  !
END SUBROUTINE COMHES
