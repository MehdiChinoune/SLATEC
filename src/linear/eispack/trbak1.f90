!** TRBAK1
SUBROUTINE TRBAK1(Nm,N,A,E,M,Z)
  !>
  !***
  !  Form the eigenvectors of real symmetric matrix from
  !            the eigenvectors of a symmetric tridiagonal matrix formed
  !            by TRED1.
  !***
  ! **Library:**   SLATEC (EISPACK)
  !***
  ! **Category:**  D4C4
  !***
  ! **Type:**      SINGLE PRECISION (TRBAK1-S)
  !***
  ! **Keywords:**  EIGENVECTORS OF A REAL SYMMETRIC MATRIX, EISPACK
  !***
  ! **Author:**  Smith, B. T., et al.
  !***
  ! **Description:**
  !
  !     This subroutine is a translation of the ALGOL procedure TRBAK1,
  !     NUM. MATH. 11, 181-195(1968) by Martin, Reinsch, and Wilkinson.
  !     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
  !
  !     This subroutine forms the eigenvectors of a REAL SYMMETRIC
  !     matrix by back transforming those of the corresponding
  !     symmetric tridiagonal matrix determined by  TRED1.
  !
  !     On Input
  !
  !        NM must be set to the row dimension of the two-dimensional
  !          array parameters, A and Z, as declared in the calling
  !          program dimension statement.  NM is an INTEGER variable.
  !
  !        N is the order of the matrix.  N is an INTEGER variable.
  !          N must be less than or equal to NM.
  !
  !        A contains information about the orthogonal transformations
  !          used in the reduction by  TRED1  in its strict lower
  !          triangle.  A is a two-dimensional REAL array, dimensioned
  !          A(NM,N).
  !
  !        E contains the subdiagonal elements of the tridiagonal matrix
  !          in its last N-1 positions.  E(1) is arbitrary.  These
  !          elements provide the remaining information about the
  !          orthogonal transformations.  E is a one-dimensional REAL
  !          array, dimensioned E(N).
  !
  !        M is the number of columns of Z to be back transformed.
  !          M is an INTEGER variable.
  !
  !        Z contains the eigenvectors to be back transformed in its
  !          first M columns.  Z is a two-dimensional REAL array,
  !          dimensioned Z(NM,M).
  !
  !     On Output
  !
  !        Z contains the transformed eigenvectors in its first M columns.
  !
  !     Note that TRBAK1 preserves vector Euclidean norms.
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
  INTEGER i, j, k, l, M, N, Nm
  REAL A(Nm,*), E(*), Z(Nm,*)
  REAL s
  !
  !* FIRST EXECUTABLE STATEMENT  TRBAK1
  IF ( M/=0 ) THEN
    IF ( N/=1 ) THEN
      !
      DO i = 2, N
        l = i - 1
        IF ( E(i)/=0.0E0 ) THEN
          !
          DO j = 1, M
            s = 0.0E0
            !
            DO k = 1, l
              s = s + A(i,k)*Z(k,j)
            END DO
            !     .......... DIVISOR BELOW IS NEGATIVE OF H FORMED IN TRED1.
            !                DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
            s = (s/A(i,l))/E(i)
            !
            DO k = 1, l
              Z(k,j) = Z(k,j) + s*A(i,k)
            END DO
            !
          END DO
        END IF
        !
      END DO
    END IF
  END IF
  !
END SUBROUTINE TRBAK1
