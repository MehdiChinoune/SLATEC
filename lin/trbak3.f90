!*==TRBAK3.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK TRBAK3
SUBROUTINE TRBAK3(Nm,N,Nv,A,M,Z)
  IMPLICIT NONE
  !*--TRBAK35
  !***BEGIN PROLOGUE  TRBAK3
  !***PURPOSE  Form the eigenvectors of a real symmetric matrix from the
  !            eigenvectors of a symmetric tridiagonal matrix formed
  !            by TRED3.
  !***LIBRARY   SLATEC (EISPACK)
  !***CATEGORY  D4C4
  !***TYPE      SINGLE PRECISION (TRBAK3-S)
  !***KEYWORDS  EIGENVECTORS OF A REAL SYMMETRIC MATRIX, EISPACK
  !***AUTHOR  Smith, B. T., et al.
  !***DESCRIPTION
  !
  !     This subroutine is a translation of the ALGOL procedure TRBAK3,
  !     NUM. MATH. 11, 181-195(1968) by Martin, Reinsch, and Wilkinson.
  !     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
  !
  !     This subroutine forms the eigenvectors of a REAL SYMMETRIC
  !     matrix by back transforming those of the corresponding
  !     symmetric tridiagonal matrix determined by  TRED3.
  !
  !     On Input
  !
  !        NM must be set to the row dimension of the two-dimensional
  !          array parameter, Z, as declared in the calling program
  !          dimension statement.  NM is an INTEGER variable.
  !
  !        N is the order of the matrix.  N is an INTEGER variable.
  !          N must be less than or equal to NM.
  !
  !        NV is an INTEGER variable set equal to the dimension of the
  !          array A as specified in the calling program.  NV must not
  !          be less than  N*(N+1)/2.
  !
  !        A contains information about the orthogonal transformations
  !          used in the reduction by  TRED3  in its first N*(N+1)/2
  !          positions.  A is a one-dimensional REAL array, dimensioned
  !          A(NV).
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
  !     Note that TRBAK3 preserves vector Euclidean norms.
  !
  !     Questions and comments should be directed to b. s. Garbow,
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
  !***END PROLOGUE  TRBAK3
  !
  INTEGER i , j , k , l , M , N , ik , iz , Nm , Nv
  REAL A(*) , Z(Nm,*)
  REAL h , s
  !
  !***FIRST EXECUTABLE STATEMENT  TRBAK3
  IF ( M/=0 ) THEN
    IF ( N/=1 ) THEN
      !
      DO i = 2 , N
        l = i - 1
        iz = (i*l)/2
        ik = iz + i
        h = A(ik)
        IF ( h/=0.0E0 ) THEN
          !
          DO j = 1 , M
            s = 0.0E0
            ik = iz
            !
            DO k = 1 , l
              ik = ik + 1
              s = s + A(ik)*Z(k,j)
            ENDDO
            !     .......... DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
            s = (s/h)/h
            ik = iz
            !
            DO k = 1 , l
              ik = ik + 1
              Z(k,j) = Z(k,j) - s*A(ik)
            ENDDO
            !
          ENDDO
        ENDIF
        !
      ENDDO
    ENDIF
  ENDIF
  !
END SUBROUTINE TRBAK3
