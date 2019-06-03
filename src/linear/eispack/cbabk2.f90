!** CBABK2
SUBROUTINE CBABK2(Nm,N,Low,Igh,Scalee,M,Zr,Zi)
  !>
  !  Form the eigenvectors of a complex general matrix from the
  !            eigenvectors of matrix output from CBAL.
  !***
  ! **Library:**   SLATEC (EISPACK)
  !***
  ! **Category:**  D4C4
  !***
  ! **Type:**      COMPLEX (BALBAK-S, CBABK2-C)
  !***
  ! **Keywords:**  EIGENVECTORS, EISPACK
  !***
  ! **Author:**  Smith, B. T., et al.
  !***
  ! **Description:**
  !
  !     This subroutine is a translation of the ALGOL procedure
  !     CBABK2, which is a complex version of BALBAK,
  !     NUM. MATH. 13, 293-304(1969) by Parlett and Reinsch.
  !     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 315-326(1971).
  !
  !     This subroutine forms the eigenvectors of a COMPLEX GENERAL
  !     matrix by back transforming those of the corresponding
  !     balanced matrix determined by  CBAL.
  !
  !     On INPUT
  !
  !        NM must be set to the row dimension of the two-dimensional
  !          array parameters, ZR and ZI, as declared in the calling
  !          program dimension statement.  NM is an INTEGER variable.
  !
  !        N is the order of the matrix Z=(ZR,ZI).  N is an INTEGER
  !          variable.  N must be less than or equal to NM.
  !
  !        LOW and IGH are INTEGER variables determined by  CBAL.
  !
  !        SCALE contains information determining the permutations and
  !          scaling factors used by  CBAL.  SCALE is a one-dimensional
  !          REAL array, dimensioned SCALE(N).
  !
  !        M is the number of eigenvectors to be back transformed.
  !          M is an INTEGER variable.
  !
  !        ZR and ZI contain the real and imaginary parts, respectively,
  !          of the eigenvectors to be back transformed in their first
  !          M columns.  ZR and ZI are two-dimensional REAL arrays,
  !          dimensioned ZR(NM,M) and ZI(NM,M).
  !
  !     On OUTPUT
  !
  !        ZR and ZI contain the real and imaginary parts,
  !          respectively, of the transformed eigenvectors
  !          in their first M columns.
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
  INTEGER i, j, k, M, N, ii, Nm, Igh, Low
  REAL(SP) Scalee(*), Zr(Nm,*), Zi(Nm,*)
  REAL(SP) s
  !
  !* FIRST EXECUTABLE STATEMENT  CBABK2
  IF ( M/=0 ) THEN
    IF ( Igh/=Low ) THEN
      !
      DO i = Low, Igh
        s = Scalee(i)
        !     .......... LEFT HAND EIGENVECTORS ARE BACK TRANSFORMED
        !                IF THE FOREGOING STATEMENT IS REPLACED BY
        !                S=1.0E0/SCALE(I). ..........
        DO j = 1, M
          Zr(i,j) = Zr(i,j)*s
          Zi(i,j) = Zi(i,j)*s
        END DO
        !
      END DO
    END IF
    !     .......... FOR I=LOW-1 STEP -1 UNTIL 1,
    !                IGH+1 STEP 1 UNTIL N DO -- ..........
    DO ii = 1, N
      i = ii
      IF ( i<Low.OR.i>Igh ) THEN
        IF ( i<Low ) i = Low - ii
        k = INT(Scalee(i))
        IF ( k/=i ) THEN
          !
          DO j = 1, M
            s = Zr(i,j)
            Zr(i,j) = Zr(k,j)
            Zr(k,j) = s
            s = Zi(i,j)
            Zi(i,j) = Zi(k,j)
            Zi(k,j) = s
          END DO
        END IF
      END IF
      !
    END DO
  END IF
  !
END SUBROUTINE CBABK2
