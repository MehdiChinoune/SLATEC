!*==DPBDI.f90  processed by SPAG 6.72Dc at 10:58 on  6 Feb 2019
!DECK DPBDI
SUBROUTINE DPBDI(Abd,Lda,N,M,Det)
  IMPLICIT NONE
  !*--DPBDI5
  !***BEGIN PROLOGUE  DPBDI
  !***PURPOSE  Compute the determinant of a symmetric positive definite
  !            band matrix using the factors computed by DPBCO or DPBFA.
  !***LIBRARY   SLATEC (LINPACK)
  !***CATEGORY  D3B2
  !***TYPE      DOUBLE PRECISION (SPBDI-S, DPBDI-D, CPBDI-C)
  !***KEYWORDS  BANDED, DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK,
  !             MATRIX, POSITIVE DEFINITE
  !***AUTHOR  Moler, C. B., (U. of New Mexico)
  !***DESCRIPTION
  !
  !     DPBDI computes the determinant
  !     of a double precision symmetric positive definite band matrix
  !     using the factors computed by DPBCO or DPBFA.
  !     If the inverse is needed, use DPBSL  N  times.
  !
  !     On Entry
  !
  !        ABD     DOUBLE PRECISION(LDA, N)
  !                the output from DPBCO or DPBFA.
  !
  !        LDA     INTEGER
  !                the leading dimension of the array  ABD .
  !
  !        N       INTEGER
  !                the order of the matrix  A .
  !
  !        M       INTEGER
  !                the number of diagonals above the main diagonal.
  !
  !     On Return
  !
  !        DET     DOUBLE PRECISION(2)
  !                determinant of original matrix in the form
  !                DETERMINANT = DET(1) * 10.0**DET(2)
  !                with  1.0 .LE. DET(1) .LT. 10.0
  !                or  DET(1) .EQ. 0.0 .
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
  !***END PROLOGUE  DPBDI
  INTEGER Lda , N , M
  REAL(8) :: Abd(Lda,*)
  REAL(8) :: Det(2)
  !
  REAL(8) :: s
  INTEGER i
  !***FIRST EXECUTABLE STATEMENT  DPBDI
  !
  !     COMPUTE DETERMINANT
  !
  Det(1) = 1.0D0
  Det(2) = 0.0D0
  s = 10.0D0
  DO i = 1 , N
    Det(1) = Abd(M+1,i)**2*Det(1)
    IF ( Det(1)==0.0D0 ) EXIT
    DO WHILE ( Det(1)<1.0D0 )
      Det(1) = s*Det(1)
      Det(2) = Det(2) - 1.0D0
    ENDDO
    DO WHILE ( Det(1)>=s )
      Det(1) = Det(1)/s
      Det(2) = Det(2) + 1.0D0
    ENDDO
  ENDDO
END SUBROUTINE DPBDI
