!DECK ORTBAK
SUBROUTINE ORTBAK(Nm,Low,Igh,A,Ort,M,Z)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  ORTBAK
  !***PURPOSE  Form the eigenvectors of a general real matrix from the
  !            eigenvectors of the upper Hessenberg matrix output from
  !            ORTHES.
  !***LIBRARY   SLATEC (EISPACK)
  !***CATEGORY  D4C4
  !***TYPE      SINGLE PRECISION (ORTBAK-S, CORTB-C)
  !***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
  !***AUTHOR  Smith, B. T., et al.
  !***DESCRIPTION
  !
  !     This subroutine is a translation of the ALGOL procedure ORTBAK,
  !     NUM. MATH. 12, 349-368(1968) by Martin and Wilkinson.
  !     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
  !
  !     This subroutine forms the eigenvectors of a REAL GENERAL
  !     matrix by back transforming those of the corresponding
  !     upper Hessenberg matrix determined by  ORTHES.
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
  !        A contains some information about the orthogonal trans-
  !          formations used in the reduction to Hessenberg form by
  !          ORTHES  in its strict lower triangle.  A is a two-dimensional
  !          REAL array, dimensioned A(NM,IGH).
  !
  !        ORT contains further information about the orthogonal trans-
  !          formations used in the reduction by  ORTHES.  Only elements
  !          LOW through IGH are used.  ORT is a one-dimensional REAL
  !          array, dimensioned ORT(IGH).
  !
  !        M is the number of columns of Z to be back transformed.
  !          M is an INTEGER variable.
  !
  !        Z contains the real and imaginary parts of the eigenvectors to
  !          be back transformed in its first M columns.  Z is a two-
  !          dimensional REAL array, dimensioned Z(NM,M).
  !
  !     On OUTPUT
  !
  !        Z contains the real and imaginary parts of the transformed
  !          eigenvectors in its first M columns.
  !
  !        ORT has been used for temporary storage as is not restored.
  !
  !     NOTE that ORTBAK preserves vector Euclidean norms.
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
  !***END PROLOGUE  ORTBAK
  !
  INTEGER i, j, M, la, mm, mp, Nm, Igh, kp1, Low, mp1
  REAL A(Nm,*), Ort(*), Z(Nm,*)
  REAL g
  !
  !***FIRST EXECUTABLE STATEMENT  ORTBAK
  IF ( M/=0 ) THEN
    la = Igh - 1
    kp1 = Low + 1
    IF ( la>=kp1 ) THEN
      !     .......... FOR MP=IGH-1 STEP -1 UNTIL LOW+1 DO -- ..........
      DO mm = kp1, la
        mp = Low + Igh - mm
        IF ( A(mp,mp-1)/=0.0E0 ) THEN
          mp1 = mp + 1
          !
          DO i = mp1, Igh
            Ort(i) = A(i,mp-1)
          ENDDO
          !
          DO j = 1, M
            g = 0.0E0
            !
            DO i = mp, Igh
              g = g + Ort(i)*Z(i,j)
            ENDDO
            !     .......... DIVISOR BELOW IS NEGATIVE OF H FORMED IN ORTHES.
            !                DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
            g = (g/Ort(mp))/A(mp,mp-1)
            !
            DO i = mp, Igh
              Z(i,j) = Z(i,j) + g*Ort(i)
            ENDDO
            !
          ENDDO
        ENDIF
        !
      ENDDO
    ENDIF
  ENDIF
  !
END SUBROUTINE ORTBAK
