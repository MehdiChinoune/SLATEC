!DECK FIGI
SUBROUTINE FIGI(Nm,N,T,D,E,E2,Ierr)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  FIGI
  !***PURPOSE  Transforms certain real non-symmetric tridiagonal matrix
  !            to symmetric tridiagonal matrix.
  !***LIBRARY   SLATEC (EISPACK)
  !***CATEGORY  D4C1C
  !***TYPE      SINGLE PRECISION (FIGI-S)
  !***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
  !***AUTHOR  Smith, B. T., et al.
  !***DESCRIPTION
  !
  !     Given a NONSYMMETRIC TRIDIAGONAL matrix such that the products
  !     of corresponding pairs of off-diagonal elements are all
  !     non-negative, this subroutine reduces it to a symmetric
  !     tridiagonal matrix with the same eigenvalues.  If, further,
  !     a zero product only occurs when both factors are zero,
  !     the reduced matrix is similar to the original matrix.
  !
  !     On INPUT
  !
  !        NM must be set to the row dimension of the two-dimensional
  !          array parameter, T, as declared in the calling program
  !          dimension statement.  NM is an INTEGER variable.
  !
  !        N is the order of the matrix T.  N is an INTEGER variable.
  !          N must be less than or equal to NM.
  !
  !        T contains the nonsymmetric matrix.  Its subdiagonal is
  !          stored in the last N-1 positions of the first column,
  !          its diagonal in the N positions of the second column,
  !          and its superdiagonal in the first N-1 positions of
  !          the third column.  T(1,1) and T(N,3) are arbitrary.
  !          T is a two-dimensional REAL array, dimensioned T(NM,3).
  !
  !     On OUTPUT
  !
  !        T is unaltered.
  !
  !        D contains the diagonal elements of the tridiagonal symmetric
  !          matrix.  D is a one-dimensional REAL array, dimensioned D(N).
  !
  !        E contains the subdiagonal elements of the tridiagonal
  !          symmetric matrix in its last N-1 positions.  E(1) is not set.
  !          E is a one-dimensional REAL array, dimensioned E(N).
  !
  !        E2 contains the squares of the corresponding elements of E.
  !          E2 may coincide with E if the squares are not needed.
  !          E2 is a one-dimensional REAL array, dimensioned E2(N).
  !
  !        IERR is an INTEGER flag set to
  !          Zero       for normal return,
  !          N+I        if T(I,1)*T(I-1,3) is negative and a symmetric
  !                     matrix cannot be produced with FIGI,
  !          -(3*N+I)   if T(I,1)*T(I-1,3) is zero with one factor
  !                     non-zero.  In this case, the eigenvectors of
  !                     the symmetric matrix are not simply related
  !                     to those of  T  and should not be sought.
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
  !***END PROLOGUE  FIGI
  !
  INTEGER i, N, Nm, Ierr
  REAL T(Nm,3), D(*), E(*), E2(*)
  !
  !***FIRST EXECUTABLE STATEMENT  FIGI
  Ierr = 0
  !
  DO i = 1, N
    IF ( i/=1 ) THEN
      E2(i) = T(i,1)*T(i-1,3)
      IF ( E2(i)<0 ) GOTO 100
      IF ( E2(i)==0 ) THEN
        !     .......... SET ERROR -- PRODUCT OF SOME PAIR OF OFF-DIAGONAL
        !                ELEMENTS IS ZERO WITH ONE MEMBER NON-ZERO ..........
        IF ( T(i,1)/=0.0E0.OR.T(i-1,3)/=0.0E0 ) Ierr = -(3*N+i)
      ENDIF
      E(i) = SQRT(E2(i))
    ENDIF
    D(i) = T(i,2)
  ENDDO
  !
  RETURN
  !     .......... SET ERROR -- PRODUCT OF SOME PAIR OF OFF-DIAGONAL
  !                ELEMENTS IS NEGATIVE ..........
  100  Ierr = N + i
  RETURN
END SUBROUTINE FIGI
