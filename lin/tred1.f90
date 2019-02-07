!*==TRED1.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK TRED1
SUBROUTINE TRED1(Nm,N,A,D,E,E2)
  IMPLICIT NONE
  !*--TRED15
  !***BEGIN PROLOGUE  TRED1
  !***PURPOSE  Reduce a real symmetric matrix to symmetric tridiagonal
  !            matrix using orthogonal similarity transformations.
  !***LIBRARY   SLATEC (EISPACK)
  !***CATEGORY  D4C1B1
  !***TYPE      SINGLE PRECISION (TRED1-S)
  !***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
  !***AUTHOR  Smith, B. T., et al.
  !***DESCRIPTION
  !
  !     This subroutine is a translation of the ALGOL procedure TRED1,
  !     NUM. MATH. 11, 181-195(1968) by Martin, Reinsch, and Wilkinson.
  !     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
  !
  !     This subroutine reduces a REAL SYMMETRIC matrix
  !     to a symmetric tridiagonal matrix using
  !     orthogonal similarity transformations.
  !
  !     On Input
  !
  !        NM must be set to the row dimension of the two-dimensional
  !          array parameter, A, as declared in the calling program
  !          dimension statement.  NM is an INTEGER variable.
  !
  !        N is the order of the matrix A.  N is an INTEGER variable.
  !          N must be less than or equal to NM.
  !
  !        A contains the real symmetric input matrix.  Only the lower
  !          triangle of the matrix need be supplied.  A is a two-
  !          dimensional REAL array, dimensioned A(NM,N).
  !
  !     On Output
  !
  !        A contains information about the orthogonal transformations
  !          used in the reduction in its strict lower triangle.  The
  !          full upper triangle of A is unaltered.
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
  !***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
  !                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
  !                 system Routines - EISPACK Guide, Springer-Verlag,
  !                 1976.
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   760101  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  TRED1
  !
  INTEGER i , j , k , l , N , ii , Nm , jp1
  REAL A(Nm,*) , D(*) , E(*) , E2(*)
  REAL f , g , h , scale
  !
  !***FIRST EXECUTABLE STATEMENT  TRED1
  DO i = 1 , N
    D(i) = A(i,i)
  ENDDO
  !     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
  DO ii = 1 , N
    i = N + 1 - ii
    l = i - 1
    h = 0.0E0
    scale = 0.0E0
    IF ( l>=1 ) THEN
      !     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
      DO k = 1 , l
        scale = scale + ABS(A(i,k))
      ENDDO
      !
      IF ( scale/=0.0E0 ) THEN
        !
        DO k = 1 , l
          A(i,k) = A(i,k)/scale
          h = h + A(i,k)*A(i,k)
        ENDDO
        !
        E2(i) = scale*scale*h
        f = A(i,l)
        g = -SIGN(SQRT(h),f)
        E(i) = scale*g
        h = h - f*g
        A(i,l) = f - g
        IF ( l/=1 ) THEN
          f = 0.0E0
          !
          DO j = 1 , l
            g = 0.0E0
            !     .......... FORM ELEMENT OF A*U ..........
            DO k = 1 , j
              g = g + A(j,k)*A(i,k)
            ENDDO
            !
            jp1 = j + 1
            IF ( l>=jp1 ) THEN
              !
              DO k = jp1 , l
                g = g + A(k,j)*A(i,k)
              ENDDO
            ENDIF
            !     .......... FORM ELEMENT OF P ..........
            E(j) = g/h
            f = f + E(j)*A(i,j)
          ENDDO
          !
          h = f/(h+h)
          !     .......... FORM REDUCED A ..........
          DO j = 1 , l
            f = A(i,j)
            g = E(j) - h*f
            E(j) = g
            !
            DO k = 1 , j
              A(j,k) = A(j,k) - f*E(k) - g*A(i,k)
            ENDDO
          ENDDO
        ENDIF
        !
        DO k = 1 , l
          A(i,k) = scale*A(i,k)
        ENDDO
        GOTO 50
      ENDIF
    ENDIF
    E(i) = 0.0E0
    E2(i) = 0.0E0
    !
    50     h = D(i)
    D(i) = A(i,i)
    A(i,i) = h
  ENDDO
  !
END SUBROUTINE TRED1
