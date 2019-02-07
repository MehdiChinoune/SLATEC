!*==CORTB.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK CORTB
SUBROUTINE CORTB(Nm,Low,Igh,Ar,Ai,Ortr,Orti,M,Zr,Zi)
  IMPLICIT NONE
  !*--CORTB5
  !***BEGIN PROLOGUE  CORTB
  !***PURPOSE  Form the eigenvectors of a complex general matrix from
  !            eigenvectors of upper Hessenberg matrix output from
  !            CORTH.
  !***LIBRARY   SLATEC (EISPACK)
  !***CATEGORY  D4C4
  !***TYPE      COMPLEX (ORTBAK-S, CORTB-C)
  !***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
  !***AUTHOR  Smith, B. T., et al.
  !***DESCRIPTION
  !
  !     This subroutine is a translation of a complex analogue of
  !     the ALGOL procedure ORTBAK, NUM. MATH. 12, 349-368(1968)
  !     by Martin and Wilkinson.
  !     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
  !
  !     This subroutine forms the eigenvectors of a COMPLEX GENERAL
  !     matrix by back transforming those of the corresponding
  !     upper Hessenberg matrix determined by  CORTH.
  !
  !     On INPUT
  !
  !        NM must be set to the row dimension of the two-dimensional
  !          array parameters, AR, AI, ZR, and ZI, as declared in the
  !          calling program dimension statement.  NM is an INTEGER
  !          variable.
  !
  !        LOW and IGH are two INTEGER variables determined by the
  !          balancing subroutine  CBAL.  If  CBAL  has not been used,
  !          set LOW=1 and IGH equal to the order of the matrix.
  !
  !        AR and AI contain information about the unitary trans-
  !          formations used in the reduction by  CORTH  in their
  !          strict lower triangles.  AR and AI are two-dimensional
  !          REAL arrays, dimensioned AR(NM,IGH) and AI(NM,IGH).
  !
  !        ORTR and ORTI contain further information about the unitary
  !          transformations used in the reduction by  CORTH.  Only
  !          elements LOW through IGH are used.  ORTR and ORTI are
  !          one-dimensional REAL arrays, dimensioned ORTR(IGH) and
  !          ORTI(IGH).
  !
  !        M is the number of columns of Z=(ZR,ZI) to be back transformed.
  !          M is an INTEGER variable.
  !
  !        ZR and ZI contain the real and imaginary parts, respectively,
  !          of the eigenvectors to be back transformed in their first
  !          M columns.  ZR and ZI are two-dimensional REAL arrays,
  !          dimensioned ZR(NM,M) and ZI(NM,M).
  !
  !     On OUTPUT
  !
  !        ZR and ZI contain the real and imaginary parts, respectively,
  !          of the transformed eigenvectors in their first M columns.
  !
  !        ORTR and ORTI have been altered.
  !
  !     Note that CORTB preserves vector Euclidean norms.
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
  !***END PROLOGUE  CORTB
  !
  INTEGER i , j , M , la , mm , mp , Nm , Igh , kp1 , Low , mp1
  REAL Ar(Nm,*) , Ai(Nm,*) , Ortr(*) , Orti(*)
  REAL Zr(Nm,*) , Zi(Nm,*)
  REAL h , gi , gr
  !
  !***FIRST EXECUTABLE STATEMENT  CORTB
  IF ( M/=0 ) THEN
    la = Igh - 1
    kp1 = Low + 1
    IF ( la>=kp1 ) THEN
      !     .......... FOR MP=IGH-1 STEP -1 UNTIL LOW+1 DO -- ..........
      DO mm = kp1 , la
        mp = Low + Igh - mm
        IF ( Ar(mp,mp-1)/=0.0E0.OR.Ai(mp,mp-1)/=0.0E0 ) THEN
          !     .......... H BELOW IS NEGATIVE OF H FORMED IN CORTH ..........
          h = Ar(mp,mp-1)*Ortr(mp) + Ai(mp,mp-1)*Orti(mp)
          mp1 = mp + 1
          !
          DO i = mp1 , Igh
            Ortr(i) = Ar(i,mp-1)
            Orti(i) = Ai(i,mp-1)
          ENDDO
          !
          DO j = 1 , M
            gr = 0.0E0
            gi = 0.0E0
            !
            DO i = mp , Igh
              gr = gr + Ortr(i)*Zr(i,j) + Orti(i)*Zi(i,j)
              gi = gi + Ortr(i)*Zi(i,j) - Orti(i)*Zr(i,j)
            ENDDO
            !
            gr = gr/h
            gi = gi/h
            !
            DO i = mp , Igh
              Zr(i,j) = Zr(i,j) + gr*Ortr(i) - gi*Orti(i)
              Zi(i,j) = Zi(i,j) + gr*Orti(i) + gi*Ortr(i)
            ENDDO
            !
          ENDDO
        ENDIF
        !
      ENDDO
    ENDIF
  ENDIF
  !
END SUBROUTINE CORTB
