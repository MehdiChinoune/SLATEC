!** BALBAK
SUBROUTINE BALBAK(Nm,N,Low,Igh,Scalee,M,Z)
  !>
  !***
  !  Form the eigenvectors of a real general matrix from the
  !            eigenvectors of matrix output from BALANC.
  !***
  ! **Library:**   SLATEC (EISPACK)
  !***
  ! **Category:**  D4C4
  !***
  ! **Type:**      SINGLE PRECISION (BALBAK-S, CBABK2-C)
  !***
  ! **Keywords:**  EIGENVECTORS, EISPACK
  !***
  ! **Author:**  Smith, B. T., et al.
  !***
  ! **Description:**
  !
  !     This subroutine is a translation of the ALGOL procedure BALBAK,
  !     NUM. MATH. 13, 293-304(1969) by Parlett and Reinsch.
  !     HANDBOOK FOR AUTO. COMP., Vol.II-LINEAR ALGEBRA, 315-326(1971).
  !
  !     This subroutine forms the eigenvectors of a REAL GENERAL
  !     matrix by back transforming those of the corresponding
  !     balanced matrix determined by  BALANC.
  !
  !     On INPUT
  !
  !        NM must be set to the row dimension of the two-dimensional
  !          array parameter, Z, as declared in the calling program
  !          dimension statement.  NM is an INTEGER variable.
  !
  !        N is the number of components of the vectors in matrix Z.
  !          N is an INTEGER variable.  N must be less than or equal
  !          to NM.
  !
  !        LOW and IGH are INTEGER variables determined by  BALANC.
  !
  !        SCALE contains information determining the permutations and
  !          scaling factors used by  BALANC.  SCALE is a one-dimensional
  !          REAL array, dimensioned SCALE(N).
  !
  !        M is the number of columns of Z to be back transformed.
  !          M is an INTEGER variable.
  !
  !        Z contains the real and imaginary parts of the eigen-
  !          vectors to be back transformed in its first M columns.
  !          Z is a two-dimensional REAL array, dimensioned Z(NM,M).
  !
  !     On OUTPUT
  !
  !        Z contains the real and imaginary parts of the
  !          transformed eigenvectors in its first M columns.
  !
  !     Questions and comments should be directed to B. S. Garbow,
  !     Applied Mathematics Division, ARGONNE NATIONAL LABORATORY
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
  INTEGER i, j, k, M, N, ii, Nm, Igh, Low
  REAL Scalee(*), Z(Nm,*)
  REAL s
  !
  !* FIRST EXECUTABLE STATEMENT  BALBAK
  IF ( M/=0 ) THEN
    IF ( Igh/=Low ) THEN
      !
      DO i = Low, Igh
        s = Scalee(i)
        !     .......... LEFT HAND EIGENVECTORS ARE BACK TRANSFORMED
        !                IF THE FOREGOING STATEMENT IS REPLACED BY
        !                S=1.0E0/SCALE(I). ..........
        DO j = 1, M
          Z(i,j) = Z(i,j)*s
        END DO
        !
      END DO
    END IF
    !     ......... FOR I=LOW-1 STEP -1 UNTIL 1,
    !               IGH+1 STEP 1 UNTIL N DO -- ..........
    DO ii = 1, N
      i = ii
      IF ( i<Low.OR.i>Igh ) THEN
        IF ( i<Low ) i = Low - ii
        k = INT(Scalee(i))
        IF ( k/=i ) THEN
          !
          DO j = 1, M
            s = Z(i,j)
            Z(i,j) = Z(k,j)
            Z(k,j) = s
          END DO
        END IF
      END IF
      !
    END DO
  END IF
  !
END SUBROUTINE BALBAK
