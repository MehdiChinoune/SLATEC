!** SPBDI
SUBROUTINE SPBDI(Abd,Lda,N,M,Det)
  !>
  !  Compute the determinant of a symmetric positive definite
  !            band matrix using the factors computed by SPBCO or SPBFA.
  !***
  ! **Library:**   SLATEC (LINPACK)
  !***
  ! **Category:**  D3B2
  !***
  ! **Type:**      SINGLE PRECISION (SPBDI-S, DPBDI-D, CPBDI-C)
  !***
  ! **Keywords:**  BANDED, DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK,
  !             MATRIX, POSITIVE DEFINITE
  !***
  ! **Author:**  Moler, C. B., (U. of New Mexico)
  !***
  ! **Description:**
  !
  !     SPBDI computes the determinant
  !     of a real symmetric positive definite band matrix
  !     using the factors computed by SPBCO or SPBFA.
  !     If the inverse is needed, use SPBSL  N  times.
  !
  !     On Entry
  !
  !        ABD     REAL(LDA, N)
  !                the output from SPBCO or SPBFA.
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
  !        DET     REAL(2)
  !                determinant of original matrix in the form
  !                Determinant = DET(1) * 10.0**DET(2)
  !                with  1.0 .LE. DET(1) .LT. 10.0
  !                or  DET(1) .EQ. 0.0 .
  !
  !***
  ! **References:**  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
  !                 Stewart, LINPACK Users' Guide, SIAM, 1979.
  !***
  ! **Routines called:**  (NONE)

  !* REVISION HISTORY  (YYMMDD)
  !   780814  DATE WRITTEN
  !   890831  Modified array declarations.  (WRB)
  !   890831  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920501  Reformatted the REFERENCES section.  (WRB)
  
  INTEGER Lda, N, M
  REAL(SP) Abd(Lda,*)
  REAL(SP) Det(2)
  !
  REAL(SP) s
  INTEGER i
  !* FIRST EXECUTABLE STATEMENT  SPBDI
  !
  !     COMPUTE DETERMINANT
  !
  Det(1) = 1.0E0
  Det(2) = 0.0E0
  s = 10.0E0
  DO i = 1, N
    Det(1) = Abd(M+1,i)**2*Det(1)
    IF ( Det(1)==0.0E0 ) EXIT
    DO WHILE ( Det(1)<1.0E0 )
      Det(1) = s*Det(1)
      Det(2) = Det(2) - 1.0E0
    END DO
    DO WHILE ( Det(1)>=s )
      Det(1) = Det(1)/s
      Det(2) = Det(2) + 1.0E0
    END DO
  END DO
END SUBROUTINE SPBDI
