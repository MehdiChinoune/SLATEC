!** COMBAK
SUBROUTINE COMBAK(Nm,Low,Igh,Ar,Ai,Intt,M,Zr,Zi)
  IMPLICIT NONE
  !>
  !***
  !  Form the eigenvectors of a complex general matrix from the
  !            eigenvectors of a upper Hessenberg matrix output from
  !            COMHES.
  !***
  ! **Library:**   SLATEC (EISPACK)
  !***
  ! **Category:**  D4C4
  !***
  ! **Type:**      COMPLEX (ELMBAK-S, COMBAK-C)
  !***
  ! **Keywords:**  EIGENVALUES, EIGENVECTORS, EISPACK
  !***
  ! **Author:**  Smith, B. T., et al.
  !***
  ! **Description:**
  !
  !     This subroutine is a translation of the ALGOL procedure COMBAK,
  !     NUM. MATH. 12, 349-368(1968) by Martin and Wilkinson.
  !     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
  !
  !     This subroutine forms the eigenvectors of a COMPLEX GENERAL
  !     matrix by back transforming those of the corresponding
  !     upper Hessenberg matrix determined by  COMHES.
  !
  !     On INPUT
  !
  !        NM must be set to the row dimension of the two-dimensional
  !          array parameters, AR, AI, ZR and ZI, as declared in the
  !          calling program dimension statement.  NM is an INTEGER
  !          variable.
  !
  !        LOW and IGH are two INTEGER variables determined by the
  !          balancing subroutine  CBAL.  If  CBAL  has not been used,
  !          set LOW=1 and IGH equal to the order of the matrix.
  !
  !        AR and AI contain the multipliers which were used in the
  !           reduction by  COMHES  in their lower triangles below
  !           the subdiagonal.  AR and AI are two-dimensional REAL
  !           arrays, dimensioned AR(NM,IGH) and AI(NM,IGH).
  !
  !        INT contains information on the rows and columns
  !          interchanged in the reduction by  COMHES.  Only
  !          elements LOW through IGH are used.  INT is a
  !          one-dimensional INTEGER array, dimensioned INT(IGH).
  !
  !        M is the number of eigenvectors to be back transformed.
  !          M is an INTEGER variable.
  !
  !        ZR and ZI contain the real and imaginary parts, respectively,
  !          of the eigenvectors to be back transformed in their first M
  !          columns.  ZR and ZI are two-dimensional REAL arrays,
  !          dimensioned ZR(NM,M) and ZI(NM,M).
  !
  !     On OUTPUT
  !
  !        ZR and ZI contain the real and imaginary parts, respectively,
  !          of the transformed eigenvectors in their first M columns.
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
  INTEGER i, j, M, la, mm, mp, Nm, Igh, kp1, Low, mp1
  REAL Ar(Nm,*), Ai(Nm,*), Zr(Nm,*), Zi(Nm,*)
  REAL xr, xi
  INTEGER Intt(*)
  !
  !* FIRST EXECUTABLE STATEMENT  COMBAK
  IF ( M/=0 ) THEN
    la = Igh - 1
    kp1 = Low + 1
    IF ( la>=kp1 ) THEN
      !     .......... FOR MP=IGH-1 STEP -1 UNTIL LOW+1 DO -- ..........
      DO mm = kp1, la
        mp = Low + Igh - mm
        mp1 = mp + 1
        !
        DO i = mp1, Igh
          xr = Ar(i,mp-1)
          xi = Ai(i,mp-1)
          IF ( xr/=0.0E0.OR.xi/=0.0E0 ) THEN
            !
            DO j = 1, M
              Zr(i,j) = Zr(i,j) + xr*Zr(mp,j) - xi*Zi(mp,j)
              Zi(i,j) = Zi(i,j) + xr*Zi(mp,j) + xi*Zr(mp,j)
            END DO
          END IF
          !
        END DO
        !
        i = Intt(mp)
        IF ( i/=mp ) THEN
          !
          DO j = 1, M
            xr = Zr(i,j)
            Zr(i,j) = Zr(mp,j)
            Zr(mp,j) = xr
            xi = Zi(i,j)
            Zi(i,j) = Zi(mp,j)
            Zi(mp,j) = xi
          END DO
        END IF
        !
      END DO
    END IF
  END IF
  !
END SUBROUTINE COMBAK
