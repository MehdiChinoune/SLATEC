!DECK DGBDI
SUBROUTINE DGBDI(Abd,Lda,N,Ml,Mu,Ipvt,Det)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  DGBDI
  !***PURPOSE  Compute the determinant of a band matrix using the factors
  !            computed by DGBCO or DGBFA.
  !***LIBRARY   SLATEC (LINPACK)
  !***CATEGORY  D3A2
  !***TYPE      DOUBLE PRECISION (SGBDI-S, DGBDI-D, CGBDI-C)
  !***KEYWORDS  BANDED, DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK,
  !             MATRIX
  !***AUTHOR  Moler, C. B., (U. of New Mexico)
  !***DESCRIPTION
  !
  !     DGBDI computes the determinant of a band matrix
  !     using the factors computed by DGBCO or DGBFA.
  !     If the inverse is needed, use DGBSL  N  times.
  !
  !     On Entry
  !
  !        ABD     DOUBLE PRECISION(LDA, N)
  !                the output from DGBCO or DGBFA.
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
  !                the pivot vector from DGBCO or DGBFA.
  !
  !     On Return
  !
  !        DET     DOUBLE PRECISION(2)
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
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  !***END PROLOGUE  DGBDI
  INTEGER Lda, N, Ml, Mu, Ipvt(*)
  REAL(8) :: Abd(Lda,*), Det(2)
  !
  REAL(8) :: ten
  INTEGER i, m
  !***FIRST EXECUTABLE STATEMENT  DGBDI
  m = Ml + Mu + 1
  Det(1) = 1.0D0
  Det(2) = 0.0D0
  ten = 10.0D0
  DO i = 1, N
    IF ( Ipvt(i)/=i ) Det(1) = -Det(1)
    Det(1) = Abd(m,i)*Det(1)
    IF ( Det(1)==0.0D0 ) EXIT
    DO WHILE ( ABS(Det(1))<1.0D0 )
      Det(1) = ten*Det(1)
      Det(2) = Det(2) - 1.0D0
    ENDDO
    DO WHILE ( ABS(Det(1))>=ten )
      Det(1) = Det(1)/ten
      Det(2) = Det(2) + 1.0D0
    ENDDO
  ENDDO
END SUBROUTINE DGBDI
