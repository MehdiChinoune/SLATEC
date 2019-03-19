!** ELMBAK
SUBROUTINE ELMBAK(Nm,Low,Igh,A,Int,M,Z)
  IMPLICIT NONE
  !>
  !***
  !  Form the eigenvectors of a real general matrix from the
  !            eigenvectors of the upper Hessenberg matrix output from
  !            ELMHES.
  !***
  ! **Library:**   SLATEC (EISPACK)
  !***
  ! **Category:**  D4C4
  !***
  ! **Type:**      SINGLE PRECISION (ELMBAK-S, COMBAK-C)
  !***
  ! **Keywords:**  EIGENVALUES, EIGENVECTORS, EISPACK
  !***
  ! **Author:**  Smith, B. T., et al.
  !***
  ! **Description:**
  !
  !     This subroutine is a translation of the ALGOL procedure ELMBAK,
  !     NUM. MATH. 12, 349-368(1968) by Martin and Wilkinson.
  !     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
  !
  !     This subroutine forms the eigenvectors of a REAL GENERAL
  !     matrix by back transforming those of the corresponding
  !     upper Hessenberg matrix determined by  ELMHES.
  !
  !     On INPUT
  !
  !        NM must be set to the row dimension of the two-dimensional
  !          array parameters, A and Z, as declared in the calling
  !          program dimension statement.  NM is an INTEGER variable.
  !
  !        LOW and IGH are two INTEGER variables determined by the
  !          balancing subroutine  BALANC.  If  BALANC  has not been
  !          used, set LOW=1 and IGH equal to the order of the matrix.
  !
  !        A contains the multipliers which were used in the reduction
  !          by  ELMHES  in its lower triangle below the subdiagonal.
  !          A is a two-dimensional REAL array, dimensioned A(NM,IGH).
  !
  !        INT contains information on the rows and columns interchanged
  !          in the reduction by  ELMHES.  Only elements LOW through IGH
  !          are used.  INT is a one-dimensional INTEGER array,
  !          dimensioned INT(IGH).
  !
  !        M is the number of columns of Z to be back transformed.
  !          M is an INTEGER variable.
  !
  !        Z contains the real and imaginary parts of the eigenvectors
  !          to be back transformed in its first M columns.  Z is a
  !          two-dimensional REAL array, dimensioned Z(NM,M).
  !
  !     On OUTPUT
  !
  !        Z contains the real and imaginary parts of the transformed
  !          eigenvectors in its first M columns.
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
  REAL A(Nm,*), Z(Nm,*)
  REAL x
  INTEGER Int(*)
  !
  !* FIRST EXECUTABLE STATEMENT  ELMBAK
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
          x = A(i,mp-1)
          IF ( x/=0.0E0 ) THEN
            !
            DO j = 1, M
              Z(i,j) = Z(i,j) + x*Z(mp,j)
            ENDDO
          ENDIF
          !
        ENDDO
        !
        i = Int(mp)
        IF ( i/=mp ) THEN
          !
          DO j = 1, M
            x = Z(i,j)
            Z(i,j) = Z(mp,j)
            Z(mp,j) = x
          ENDDO
        ENDIF
        !
      ENDDO
    ENDIF
  ENDIF
  !
END SUBROUTINE ELMBAK
