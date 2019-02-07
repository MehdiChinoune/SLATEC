!*==SGBDI.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK SGBDI
SUBROUTINE SGBDI(Abd,Lda,N,Ml,Mu,Ipvt,Det)
  IMPLICIT NONE
  !*--SGBDI5
  !***BEGIN PROLOGUE  SGBDI
  !***PURPOSE  Compute the determinant of a band matrix using the factors
  !            computed by SGBCO or SGBFA.
  !***LIBRARY   SLATEC (LINPACK)
  !***CATEGORY  D3A2
  !***TYPE      SINGLE PRECISION (SGBDI-S, DGBDI-D, CGBDI-C)
  !***KEYWORDS  BANDED, DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK,
  !             MATRIX
  !***AUTHOR  Moler, C. B., (U. of New Mexico)
  !***DESCRIPTION
  !
  !     SGBDI computes the determinant of a band matrix
  !     using the factors computed by SBGCO or SGBFA.
  !     If the inverse is needed, use SGBSL  N  times.
  !
  !     On Entry
  !
  !        ABD     REAL(LDA, N)
  !                the output from SBGCO or SGBFA.
  !
  !        LDA     INTEGER
  !                the leading dimension of the array  ABD .
  !
  !        N       INTEGER
  !                the order of the original matrix.
  !
  !        ML      INTEGER
  !                number of diagonals below the main diagonal.
  !
  !        MU      INTEGER
  !                number of diagonals above the main diagonal.
  !
  !        IPVT    INTEGER(N)
  !                the pivot vector from SBGCO or SGBFA.
  !
  !     On Return
  !
  !        DET     REAL(2)
  !                determinant of original matrix.
  !                Determinant = DET(1) * 10.0**DET(2)
  !                with  1.0 .LE. ABS(DET(1)) .LT. 10.0
  !                or  DET(1) = 0.0 .
  !
  !***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***ROUTINES CALLED  (NONE)
  !***REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  SGBDI
  INTEGER Lda , N , Ml , Mu , Ipvt(*)
  REAL Abd(Lda,*) , Det(2)
  !
  REAL ten
  INTEGER i , m
  !***FIRST EXECUTABLE STATEMENT  SGBDI
  m = Ml + Mu + 1
  Det(1) = 1.0E0
  Det(2) = 0.0E0
  ten = 10.0E0
  DO i = 1 , N
    IF ( Ipvt(i)/=i ) Det(1) = -Det(1)
    Det(1) = Abd(m,i)*Det(1)
    IF ( Det(1)==0.0E0 ) EXIT
    DO WHILE ( ABS(Det(1))<1.0E0 )
      Det(1) = ten*Det(1)
      Det(2) = Det(2) - 1.0E0
    ENDDO
    DO WHILE ( ABS(Det(1))>=ten )
      Det(1) = Det(1)/ten
      Det(2) = Det(2) + 1.0E0
    ENDDO
  ENDDO
END SUBROUTINE SGBDI
