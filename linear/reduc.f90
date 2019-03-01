!DECK REDUC
SUBROUTINE REDUC(Nm,N,A,B,Dl,Ierr)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  REDUC
  !***PURPOSE  Reduce a generalized symmetric eigenproblem to a standard
  !            symmetric eigenproblem using Cholesky factorization.
  !***LIBRARY   SLATEC (EISPACK)
  !***CATEGORY  D4C1C
  !***TYPE      SINGLE PRECISION (REDUC-S)
  !***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
  !***AUTHOR  Smith, B. T., et al.
  !***DESCRIPTION
  !
  !     This subroutine is a translation of the ALGOL procedure REDUC1,
  !     NUM. MATH. 11, 99-110(1968) by Martin and Wilkinson.
  !     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 303-314(1971).
  !
  !     This subroutine reduces the generalized SYMMETRIC eigenproblem
  !     Ax=(LAMBDA)Bx, where B is POSITIVE DEFINITE, to the standard
  !     symmetric eigenproblem using the Cholesky factorization of B.
  !
  !     On Input
  !
  !        NM must be set to the row dimension of the two-dimensional
  !          array parameters, A and B, as declared in the calling
  !          program dimension statement.  NM is an INTEGER variable.
  !
  !        N is the order of the matrices A and B.  If the Cholesky
  !          factor L of B is already available, N should be prefixed
  !          with a minus sign.  N is an INTEGER variable.
  !
  !        A and B contain the real symmetric input matrices.  Only
  !          the full upper triangles of the matrices need be supplied.
  !          If N is negative, the strict lower triangle of B contains,
  !          instead, the strict lower triangle of its Cholesky factor L.
  !          A and B are two-dimensional REAL arrays, dimensioned A(NM,N)
  !          and B(NM,N).
  !
  !       DL contains, if N is negative, the diagonal elements of L.
  !          DL is a one-dimensional REAL array, dimensioned DL(N).
  !
  !     On Output
  !
  !        A contains in its full lower triangle the full lower triangle
  !          of the symmetric matrix derived from the reduction to the
  !          standard form.  The strict upper triangle of A is unaltered.
  !
  !        B contains in its strict lower triangle the strict lower
  !          triangle of its Cholesky factor L.  The full upper triangle
  !          of B is unaltered.
  !
  !        DL contains the diagonal elements of L.
  !
  !        IERR is an INTEGER flag set to
  !          Zero       for normal return,
  !          7*N+1      if B is not positive definite.
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
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  REDUC
  !
  INTEGER i, j, k, N, i1, j1, Nm, nn, Ierr
  REAL A(Nm,*), B(Nm,*), Dl(*)
  REAL x, y
  !
  !***FIRST EXECUTABLE STATEMENT  REDUC
  Ierr = 0
  nn = ABS(N)
  IF ( N>=0 ) THEN
    !     .......... FORM L IN THE ARRAYS B AND DL ..........
    DO i = 1, N
      i1 = i - 1
      !
      DO j = i, N
        x = B(i,j)
        IF ( i/=1 ) THEN
          !
          DO k = 1, i1
            x = x - B(i,k)*B(j,k)
          ENDDO
        ENDIF
        !
        IF ( j/=i ) THEN
          B(j,i) = x/y
        ELSE
          IF ( x<=0.0E0 ) GOTO 100
          y = SQRT(x)
          Dl(i) = y
        ENDIF
      ENDDO
    ENDDO
  ENDIF
  !     .......... FORM THE TRANSPOSE OF THE UPPER TRIANGLE OF INV(L)*A
  !                IN THE LOWER TRIANGLE OF THE ARRAY A ..........
  DO i = 1, nn
    i1 = i - 1
    y = Dl(i)
    !
    DO j = i, nn
      x = A(i,j)
      IF ( i/=1 ) THEN
        !
        DO k = 1, i1
          x = x - B(i,k)*A(j,k)
        ENDDO
      ENDIF
      !
      A(j,i) = x/y
    ENDDO
  ENDDO
  !     .......... PRE-MULTIPLY BY INV(L) AND OVERWRITE ..........
  DO j = 1, nn
    j1 = j - 1
    !
    DO i = j, nn
      x = A(i,j)
      IF ( i/=j ) THEN
        i1 = i - 1
        !
        DO k = j, i1
          x = x - A(k,j)*B(i,k)
        ENDDO
      ENDIF
      !
      IF ( j/=1 ) THEN
        !
        DO k = 1, j1
          x = x - A(j,k)*B(i,k)
        ENDDO
      ENDIF
      !
      A(i,j) = x/Dl(i)
    ENDDO
  ENDDO
  !
  RETURN
  !     .......... SET ERROR -- B IS NOT POSITIVE DEFINITE ..........
  100  Ierr = 7*N + 1
  RETURN
END SUBROUTINE REDUC
