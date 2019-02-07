!*==TRED2.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK TRED2
SUBROUTINE TRED2(Nm,N,A,D,E,Z)
  IMPLICIT NONE
  !*--TRED25
  !***BEGIN PROLOGUE  TRED2
  !***PURPOSE  Reduce a real symmetric matrix to a symmetric tridiagonal
  !            matrix using and accumulating orthogonal transformations.
  !***LIBRARY   SLATEC (EISPACK)
  !***CATEGORY  D4C1B1
  !***TYPE      SINGLE PRECISION (TRED2-S)
  !***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
  !***AUTHOR  Smith, B. T., et al.
  !***DESCRIPTION
  !
  !     This subroutine is a translation of the ALGOL procedure TRED2,
  !     NUM. MATH. 11, 181-195(1968) by Martin, Reinsch, and Wilkinson.
  !     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
  !
  !     This subroutine reduces a REAL SYMMETRIC matrix to a
  !     symmetric tridiagonal matrix using and accumulating
  !     orthogonal similarity transformations.
  !
  !     On Input
  !
  !        NM must be set to the row dimension of the two-dimensional
  !          array parameters, A and Z, as declared in the calling
  !          program dimension statement.  NM is an INTEGER variable.
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
  !        D contains the diagonal elements of the symmetric tridiagonal
  !          matrix.  D is a one-dimensional REAL array, dimensioned D(N).
  !
  !        E contains the subdiagonal elements of the symmetric
  !          tridiagonal matrix in its last N-1 positions.  E(1) is set
  !          to zero.  E is a one-dimensional REAL array, dimensioned
  !          E(N).
  !
  !        Z contains the orthogonal transformation matrix produced in
  !          the reduction.  Z is a two-dimensional REAL array,
  !          dimensioned Z(NM,N).
  !
  !        A and Z may coincide.  If distinct, A is unaltered.
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
  !***END PROLOGUE  TRED2
  !
  INTEGER i, j, k, l, N, ii, Nm, jp1
  REAL A(Nm,*), D(*), E(*), Z(Nm,*)
  REAL f, g, h, hh, scale
  !
  !***FIRST EXECUTABLE STATEMENT  TRED2
  DO i = 1, N
    !
    DO j = 1, i
      Z(i,j) = A(i,j)
    ENDDO
  ENDDO
  !
  IF ( N/=1 ) THEN
    !     .......... FOR I=N STEP -1 UNTIL 2 DO -- ..........
    DO ii = 2, N
      i = N + 2 - ii
      l = i - 1
      h = 0.0E0
      scale = 0.0E0
      IF ( l<2 ) THEN
        E(i) = Z(i,l)
      ELSE
        !     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........
        DO k = 1, l
          scale = scale + ABS(Z(i,k))
        ENDDO
        !
        IF ( scale/=0.0E0 ) THEN
          !
          DO k = 1, l
            Z(i,k) = Z(i,k)/scale
            h = h + Z(i,k)*Z(i,k)
          ENDDO
          !
          f = Z(i,l)
          g = -SIGN(SQRT(h),f)
          E(i) = scale*g
          h = h - f*g
          Z(i,l) = f - g
          f = 0.0E0
          !
          DO j = 1, l
            Z(j,i) = Z(i,j)/h
            g = 0.0E0
            !     .......... FORM ELEMENT OF A*U ..........
            DO k = 1, j
              g = g + Z(j,k)*Z(i,k)
            ENDDO
            !
            jp1 = j + 1
            IF ( l>=jp1 ) THEN
              !
              DO k = jp1, l
                g = g + Z(k,j)*Z(i,k)
              ENDDO
            ENDIF
            !     .......... FORM ELEMENT OF P ..........
            E(j) = g/h
            f = f + E(j)*Z(i,j)
          ENDDO
          !
          hh = f/(h+h)
          !     .......... FORM REDUCED A ..........
          DO j = 1, l
            f = Z(i,j)
            g = E(j) - hh*f
            E(j) = g
            !
            DO k = 1, j
              Z(j,k) = Z(j,k) - f*E(k) - g*Z(i,k)
            ENDDO
          ENDDO
        ELSE
          E(i) = Z(i,l)
        ENDIF
      ENDIF
      !
      D(i) = h
    ENDDO
  ENDIF
  !
  D(1) = 0.0E0
  E(1) = 0.0E0
  !     .......... ACCUMULATION OF TRANSFORMATION MATRICES ..........
  DO i = 1, N
    l = i - 1
    IF ( D(i)/=0.0E0 ) THEN
      !
      DO j = 1, l
        g = 0.0E0
        !
        DO k = 1, l
          g = g + Z(i,k)*Z(k,j)
        ENDDO
        !
        DO k = 1, l
          Z(k,j) = Z(k,j) - g*Z(k,i)
        ENDDO
      ENDDO
    ENDIF
    !
    D(i) = Z(i,i)
    Z(i,i) = 1.0E0
    IF ( l>=1 ) THEN
      !
      DO j = 1, l
        Z(i,j) = 0.0E0
        Z(j,i) = 0.0E0
      ENDDO
    ENDIF
    !
  ENDDO
  !
END SUBROUTINE TRED2
