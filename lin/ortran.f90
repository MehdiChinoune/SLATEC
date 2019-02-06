!*==ORTRAN.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK ORTRAN
      SUBROUTINE ORTRAN(Nm,N,Low,Igh,A,Ort,Z)
      IMPLICIT NONE
!*--ORTRAN5
!***BEGIN PROLOGUE  ORTRAN
!***PURPOSE  Accumulate orthogonal similarity transformations in the
!            reduction of real general matrix by ORTHES.
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C4
!***TYPE      SINGLE PRECISION (ORTRAN-S)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure ORTRANS,
!     NUM. MATH. 16, 181-204(1970) by Peters and Wilkinson.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
!
!     This subroutine accumulates the orthogonal similarity
!     transformations used in the reduction of a REAL GENERAL
!     matrix to upper Hessenberg form by  ORTHES.
!
!     On INPUT
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameters, A and Z, as declared in the calling
!          program dimension statement.  NM is an INTEGER variable.
!
!        N is the order of the matrix A.  N is an INTEGER variable.
!          N must be less than or equal to NM.
!
!        LOW and IGH are two INTEGER variables determined by the
!          balancing subroutine  BALANC.  If  BALANC  has not been
!          used, set LOW=1 and IGH equal to the order of the matrix, N.
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
!     On OUTPUT
!
!        Z contains the transformation matrix produced in the reduction
!          by  ORTHES  to the upper Hessenberg form.  Z is a two-
!          dimensional REAL array, dimensioned Z(NM,N).
!
!        ORT has been used for temporary storage as is not restored.
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
!***END PROLOGUE  ORTRAN
!
      INTEGER i , j , N , kl , mm , mp , Nm , Igh , Low , mp1
      REAL A(Nm,*) , Ort(*) , Z(Nm,*)
      REAL g
!
!     .......... INITIALIZE Z TO IDENTITY MATRIX ..........
!***FIRST EXECUTABLE STATEMENT  ORTRAN
      DO i = 1 , N
!
        DO j = 1 , N
          Z(i,j) = 0.0E0
        ENDDO
!
        Z(i,i) = 1.0E0
      ENDDO
!
      kl = Igh - Low - 1
      IF ( kl>=1 ) THEN
!     .......... FOR MP=IGH-1 STEP -1 UNTIL LOW+1 DO -- ..........
        DO mm = 1 , kl
          mp = Igh - mm
          IF ( A(mp,mp-1)/=0.0E0 ) THEN
            mp1 = mp + 1
!
            DO i = mp1 , Igh
              Ort(i) = A(i,mp-1)
            ENDDO
!
            DO j = mp , Igh
              g = 0.0E0
!
              DO i = mp , Igh
                g = g + Ort(i)*Z(i,j)
              ENDDO
!     .......... DIVISOR BELOW IS NEGATIVE OF H FORMED IN ORTHES.
!                DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
              g = (g/Ort(mp))/A(mp,mp-1)
!
              DO i = mp , Igh
                Z(i,j) = Z(i,j) + g*Ort(i)
              ENDDO
!
            ENDDO
          ENDIF
!
        ENDDO
      ENDIF
!
      END SUBROUTINE ORTRAN
