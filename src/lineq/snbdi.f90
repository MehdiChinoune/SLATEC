!** SNBDI
SUBROUTINE SNBDI(Abe,Lda,N,Ml,Mu,Ipvt,Det)
  IMPLICIT NONE
  !>
  !***
  !  Compute the determinant of a band matrix using the factors
  !            computed by SNBCO or SNBFA.
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  D3A2
  !***
  ! **Type:**      SINGLE PRECISION (SNBDI-S, DNBDI-D, CNBDI-C)
  !***
  ! **Keywords:**  BANDED, DETERMINANT, LINEAR EQUATIONS, NONSYMMETRIC
  !***
  ! **Author:**  Voorhees, E. A., (LANL)
  !***
  ! **Description:**
  !
  !     SNBDI computes the determinant of a band matrix
  !     using the factors computed by SNBCO or SNBFA.
  !     If the inverse is needed, use SNBSL  N  times.
  !
  !     On Entry
  !
  !        ABE     REAL(LDA, NC)
  !                the output from SNBCO or SNBFA.
  !                NC must be .GE. 2*ML+MU+1 .
  !
  !        LDA     INTEGER
  !                the leading dimension of the array  ABE .
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
  !                the pivot vector from SNBCO or SNBFA.
  !
  !     On Return
  !
  !        DET     REAL(2)
  !                determinant of original matrix.
  !                Determinant = DET(1) * 10.0**DET(2)
  !                with  1.0 .LE. ABS(DET(1)) .LT. 10.0
  !                or  DET(1) = 0.0 .
  !
  !***
  ! **References:**  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   800725  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  
  INTEGER Lda, N, Ml, Mu, Ipvt(*)
  REAL Abe(Lda,*), Det(2)
  !
  REAL ten
  INTEGER i
  !* FIRST EXECUTABLE STATEMENT  SNBDI
  Det(1) = 1.0E0
  Det(2) = 0.0E0
  ten = 10.0E0
  DO i = 1, N
    IF ( Ipvt(i)/=i ) Det(1) = -Det(1)
    Det(1) = Abe(i,Ml+1)*Det(1)
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
END SUBROUTINE SNBDI
