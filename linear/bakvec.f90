!DECK BAKVEC
SUBROUTINE BAKVEC(Nm,N,T,E,M,Z,Ierr)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  BAKVEC
  !***PURPOSE  Form the eigenvectors of a certain real non-symmetric
  !            tridiagonal matrix from a symmetric tridiagonal matrix
  !            output from FIGI.
  !***LIBRARY   SLATEC (EISPACK)
  !***CATEGORY  D4C4
  !***TYPE      SINGLE PRECISION (BAKVEC-S)
  !***KEYWORDS  EIGENVECTORS, EISPACK
  !***AUTHOR  Smith, B. T., et al.
  !***DESCRIPTION
  !
  !     This subroutine forms the eigenvectors of a NONSYMMETRIC
  !     TRIDIAGONAL matrix by back transforming those of the
  !     corresponding symmetric matrix determined by  FIGI.
  !
  !     On INPUT
  !
  !        NM must be set to the row dimension of the two-dimensional
  !          array parameters, T and Z, as declared in the calling
  !          program dimension statement.  NM is an INTEGER variable.
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
  !        E contains the subdiagonal elements of the symmetric
  !          matrix in its last N-1 positions.  E(1) is arbitrary.
  !          E is a one-dimensional REAL array, dimensioned E(N).
  !
  !        M is the number of eigenvectors to be back transformed.
  !          M is an INTEGER variable.
  !
  !        Z contains the eigenvectors to be back transformed
  !          in its first M columns.  Z is a two-dimensional REAL
  !          array, dimensioned Z(NM,M).
  !
  !     On OUTPUT
  !
  !        T is unaltered.
  !
  !        E is destroyed.
  !
  !        Z contains the transformed eigenvectors in its first M columns.
  !
  !        IERR is an INTEGER flag set to
  !          Zero       for normal return,
  !          2*N+I      if E(I) is zero with T(I,1) or T(I-1,3) non-zero.
  !                     In this case, the symmetric matrix is not similar
  !                     to the original matrix, and the eigenvectors
  !                     cannot be found by this program.
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
  !***END PROLOGUE  BAKVEC
  !
  INTEGER i, j, M, N, Nm, Ierr
  REAL T(Nm,3), E(*), Z(Nm,*)
  !
  !***FIRST EXECUTABLE STATEMENT  BAKVEC
  Ierr = 0
  IF ( M/=0 ) THEN
    E(1) = 1.0E0
    IF ( N/=1 ) THEN
      !
      DO i = 2, N
        IF ( E(i)/=0.0E0 ) THEN
          E(i) = E(i-1)*E(i)/T(i-1,3)
        ELSE
          IF ( T(i,1)/=0.0E0.OR.T(i-1,3)/=0.0E0 ) GOTO 50
          E(i) = 1.0E0
        ENDIF
      ENDDO
      !
      DO j = 1, M
        !
        DO i = 2, N
          Z(i,j) = Z(i,j)*E(i)
        ENDDO
        !
      ENDDO
    ENDIF
    RETURN
    !     .......... SET ERROR -- EIGENVECTORS CANNOT BE
    !                FOUND BY THIS PROGRAM ..........
    50     Ierr = 2*N + i
  ENDIF
  RETURN
END SUBROUTINE BAKVEC
