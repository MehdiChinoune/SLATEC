!** CGBDI
SUBROUTINE CGBDI(Abd,Lda,N,Ml,Mu,Ipvt,Det)
  !>
  !  Compute the determinant of a complex band matrix using the
  !            factors from CGBCO or CGBFA.
  !***
  ! **Library:**   SLATEC (LINPACK)
  !***
  ! **Category:**  D3C2
  !***
  ! **Type:**      COMPLEX (SGBDI-S, DGBDI-D, CGBDI-C)
  !***
  ! **Keywords:**  BANDED, DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK,
  !             MATRIX
  !***
  ! **Author:**  Moler, C. B., (U. of New Mexico)
  !***
  ! **Description:**
  !
  !     CGBDI computes the determinant of a band matrix
  !     using the factors computed by CGBCO or CGBFA.
  !     If the inverse is needed, use CGBSL  N  times.
  !
  !     On Entry
  !
  !        ABD     COMPLEX(LDA, N)
  !                the output from CGBCO or CGBFA.
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
  !                the pivot vector from CGBCO or CGBFA.
  !
  !     On Return
  !
  !        DET     COMPLEX(2)
  !                determinant of original matrix.
  !                Determinant = DET(1) * 10.0**DET(2)
  !                with  1.0 .LE. CABS1(DET(1)) .LT. 10.0
  !                or  DET(1) = 0.0 .
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
  
  INTEGER Lda, N, Ml, Mu, Ipvt(*)
  COMPLEX(SP) Abd(Lda,*), Det(2)
  !
  REAL(SP) ten
  INTEGER i, m
  !* FIRST EXECUTABLE STATEMENT  CGBDI
  m = Ml + Mu + 1
  Det(1) = (1.0E0,0.0E0)
  Det(2) = (0.0E0,0.0E0)
  ten = 10.0E0
  DO i = 1, N
    IF ( Ipvt(i)/=i ) Det(1) = -Det(1)
    Det(1) = Abd(m,i)*Det(1)
    IF ( CABS1(Det(1))==0.0E0 ) EXIT
    DO WHILE ( CABS1(Det(1))<1.0E0 )
      Det(1) = CMPLX(ten,0.0E0)*Det(1)
      Det(2) = Det(2) - (1.0E0,0.0E0)
    END DO
    DO WHILE ( CABS1(Det(1))>=ten )
      Det(1) = Det(1)/CMPLX(ten,0.0E0)
      Det(2) = Det(2) + (1.0E0,0.0E0)
    END DO
  END DO
END SUBROUTINE CGBDI
