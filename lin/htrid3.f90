!DECK HTRID3
SUBROUTINE HTRID3(Nm,N,A,D,E,E2,Tau)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  HTRID3
  !***PURPOSE  Reduce a complex Hermitian (packed) matrix to a real
  !            symmetric tridiagonal matrix by unitary similarity
  !            transformations.
  !***LIBRARY   SLATEC (EISPACK)
  !***CATEGORY  D4C1B1
  !***TYPE      SINGLE PRECISION (HTRID3-S)
  !***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
  !***AUTHOR  Smith, B. T., et al.
  !***DESCRIPTION
  !
  !     This subroutine is a translation of a complex analogue of
  !     the ALGOL procedure TRED3, NUM. MATH. 11, 181-195(1968)
  !     by Martin, Reinsch, and Wilkinson.
  !     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
  !
  !     This subroutine reduces a COMPLEX HERMITIAN matrix, stored as
  !     a single square array, to a real symmetric tridiagonal matrix
  !     using unitary similarity transformations.
  !
  !     On INPUT
  !
  !        NM must be set to the row dimension of the two-dimensional
  !          array parameter, A, as declared in the calling program
  !          dimension statement.  NM is an INTEGER variable.
  !
  !        N is the order of the matrix.  N is an INTEGER variable.
  !          N must be less than or equal to NM.
  !
  !        A contains the lower triangle of the complex Hermitian input
  !          matrix.  The real parts of the matrix elements are stored
  !          in the full lower triangle of A, and the imaginary parts
  !          are stored in the transposed positions of the strict upper
  !          triangle of A.  No storage is required for the zero
  !          imaginary parts of the diagonal elements.  A is a two-
  !          dimensional REAL array, dimensioned A(NM,N).
  !
  !     On OUTPUT
  !
  !        A contains some information about the unitary transformations
  !          used in the reduction.
  !
  !        D contains the diagonal elements of the real symmetric
  !          tridiagonal matrix.  D is a one-dimensional REAL array,
  !          dimensioned D(N).
  !
  !        E contains the subdiagonal elements of the real tridiagonal
  !          matrix in its last N-1 positions.  E(1) is set to zero.
  !          E is a one-dimensional REAL array, dimensioned E(N).
  !
  !        E2 contains the squares of the corresponding elements of E.
  !          E2(1) is set to zero.  E2 may coincide with E if the squares
  !          are not needed.  E2 is a one-dimensional REAL array,
  !          dimensioned E2(N).
  !
  !        TAU contains further information about the transformations.
  !          TAU is a one-dimensional REAL array, dimensioned TAU(2,N).
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
  !***END PROLOGUE  HTRID3
  !
  INTEGER i, j, k, l, N, ii, Nm, jm1, jp1
  REAL A(Nm,*), D(*), E(*), E2(*), Tau(2,*)
  REAL f, g, h, fi, gi, hh, si, scale
  REAL PYTHAG
  !
  !***FIRST EXECUTABLE STATEMENT  HTRID3
  Tau(1,N) = 1.0E0
  Tau(2,N) = 0.0E0
  !     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
  DO ii = 1, N
    i = N + 1 - ii
    l = i - 1
    h = 0.0E0
    scale = 0.0E0
    IF ( l>=1 ) THEN
      !     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
      DO k = 1, l
        scale = scale + ABS(A(i,k)) + ABS(A(k,i))
      ENDDO
      !
      IF ( scale/=0.0E0 ) THEN
        !
        DO k = 1, l
          A(i,k) = A(i,k)/scale
          A(k,i) = A(k,i)/scale
          h = h + A(i,k)*A(i,k) + A(k,i)*A(k,i)
        ENDDO
        !
        E2(i) = scale*scale*h
        g = SQRT(h)
        E(i) = scale*g
        f = PYTHAG(A(i,l),A(l,i))
        !     .......... FORM NEXT DIAGONAL ELEMENT OF MATRIX T ..........
        IF ( f==0.0E0 ) THEN
          Tau(1,l) = -Tau(1,i)
          si = Tau(2,i)
          A(i,l) = g
          GOTO 50
        ELSE
          Tau(1,l) = (A(l,i)*Tau(2,i)-A(i,l)*Tau(1,i))/f
          si = (A(i,l)*Tau(2,i)+A(l,i)*Tau(1,i))/f
          h = h + f*g
          g = 1.0E0 + g/f
          A(i,l) = g*A(i,l)
          A(l,i) = g*A(l,i)
          IF ( l/=1 ) GOTO 50
          GOTO 100
        ENDIF
      ELSE
        Tau(1,l) = 1.0E0
        Tau(2,l) = 0.0E0
      ENDIF
    ENDIF
    E(i) = 0.0E0
    E2(i) = 0.0E0
    GOTO 150
    50     f = 0.0E0
    !
    DO j = 1, l
      g = 0.0E0
      gi = 0.0E0
      IF ( j/=1 ) THEN
        jm1 = j - 1
        !     .......... FORM ELEMENT OF A*U ..........
        DO k = 1, jm1
          g = g + A(j,k)*A(i,k) + A(k,j)*A(k,i)
          gi = gi - A(j,k)*A(k,i) + A(k,j)*A(i,k)
        ENDDO
      ENDIF
      !
      g = g + A(j,j)*A(i,j)
      gi = gi - A(j,j)*A(j,i)
      jp1 = j + 1
      IF ( l>=jp1 ) THEN
        !
        DO k = jp1, l
          g = g + A(k,j)*A(i,k) - A(j,k)*A(k,i)
          gi = gi - A(k,j)*A(k,i) - A(j,k)*A(i,k)
        ENDDO
      ENDIF
      !     .......... FORM ELEMENT OF P ..........
      E(j) = g/h
      Tau(2,j) = gi/h
      f = f + E(j)*A(i,j) - Tau(2,j)*A(j,i)
    ENDDO
    !
    hh = f/(h+h)
    !     .......... FORM REDUCED A ..........
    DO j = 1, l
      f = A(i,j)
      g = E(j) - hh*f
      E(j) = g
      fi = -A(j,i)
      gi = Tau(2,j) - hh*fi
      Tau(2,j) = -gi
      A(j,j) = A(j,j) - 2.0E0*(f*g+fi*gi)
      IF ( j/=1 ) THEN
        jm1 = j - 1
        !
        DO k = 1, jm1
          A(j,k) = A(j,k) - f*E(k) - g*A(i,k) + fi*Tau(2,k) + gi*A(k,i)
          A(k,j) = A(k,j) - f*Tau(2,k) - g*A(k,i) - fi*E(k) - gi*A(i,k)
        ENDDO
      ENDIF
      !
    ENDDO
    !
    100    DO k = 1, l
    A(i,k) = scale*A(i,k)
    A(k,i) = scale*A(k,i)
  ENDDO
  !
  Tau(2,l) = -si
  150    D(i) = A(i,i)
  A(i,i) = scale*SQRT(h)
ENDDO
!
END SUBROUTINE HTRID3
