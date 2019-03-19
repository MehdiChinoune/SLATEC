!** REBAKB
SUBROUTINE REBAKB(Nm,N,B,Dl,M,Z)
  IMPLICIT NONE
  !>
  !***
  !  Form the eigenvectors of a generalized symmetric
  !            eigensystem from the eigenvectors of derived matrix output
  !            from REDUC2.
  !***
  ! **Library:**   SLATEC (EISPACK)
  !***
  ! **Category:**  D4C4
  !***
  ! **Type:**      SINGLE PRECISION (REBAKB-S)
  !***
  ! **Keywords:**  EIGENVALUES, EIGENVECTORS, EISPACK
  !***
  ! **Author:**  Smith, B. T., et al.
  !***
  ! **Description:**
  !
  !     This subroutine is a translation of the ALGOL procedure REBAKB,
  !     NUM. MATH. 11, 99-110(1968) by Martin and Wilkinson.
  !     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 303-314(1971).
  !
  !     This subroutine forms the eigenvectors of a generalized
  !     SYMMETRIC eigensystem by back transforming those of the
  !     derived symmetric matrix determined by  REDUC2.
  !
  !     On Input
  !
  !        NM must be set to the row dimension of the two-dimensional
  !          array parameters, B and Z, as declared in the calling
  !          program dimension statement.  NM is an INTEGER variable.
  !
  !        N is the order of the matrix system.  N is an INTEGER
  !          variable.  N must be less than or equal to NM.
  !
  !        B contains information about the similarity transformation
  !          (Cholesky decomposition) used in the reduction by  REDUC2
  !          in its strict lower triangle.  B is a two-dimensional
  !          REAL array, dimensioned B(NM,N).
  !
  !        DL contains further information about the transformation.
  !          DL is a one-dimensional REAL array, dimensioned DL(N).
  !
  !        M is the number of eigenvectors to be back transformed.
  !          M is an INTEGER variable.
  !
  !        Z contains the eigenvectors to be back transformed in its
  !          first M columns.  Z is a two-dimensional REAL array
  !          dimensioned Z(NM,M).
  !
  !     On Output
  !
  !        Z contains the transformed eigenvectors in its first
  !          M columns.
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
  INTEGER i, j, k, M, N, i1, ii, Nm
  REAL B(Nm,*), Dl(*), Z(Nm,*)
  REAL x
  !
  !* FIRST EXECUTABLE STATEMENT  REBAKB
  IF ( M/=0 ) THEN
    !
    DO j = 1, M
      !     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
      DO ii = 1, N
        i1 = N - ii
        i = i1 + 1
        x = Dl(i)*Z(i,j)
        IF ( i/=1 ) THEN
          !
          DO k = 1, i1
            x = x + B(i,k)*Z(k,j)
          ENDDO
        ENDIF
        !
        Z(i,j) = x
      ENDDO
    ENDDO
  ENDIF
  !
END SUBROUTINE REBAKB
