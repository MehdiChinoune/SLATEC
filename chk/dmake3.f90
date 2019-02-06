!*==DMAKE3.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DMAKE3
      SUBROUTINE DMAKE3(Type,Uplo,Diag,M,N,A,Nmax,Aa,Lda,Reset,Transl)
      IMPLICIT NONE
!*--DMAKE35
!***BEGIN PROLOGUE  DMAKE3
!***SUBSIDIARY
!***PURPOSE  Generate values for an M by N matrix A.
!***LIBRARY   SLATEC (BLAS)
!***AUTHOR  Dongarra, J. J., (ANL)
!           Duff, I., (AERE)
!           Du Croz, J., (NAG)
!           Hammarling, S., (NAG)
!***DESCRIPTION
!
!  Generates values for an M by N matrix A within the bandwidth
!  defined by KL and KU.
!  Stores the values in the array AA in the data structure required
!  by the routine, with unwanted elements set to rogue value.
!
!  TYPE is 'GE', 'SY' or 'TR'.
!
!  Auxiliary routine for test program for Level 3 Blas.
!***REFERENCES  (NONE)
!***ROUTINES CALLED  DBEG
!***REVISION HISTORY  (YYMMDD)
!   890208  DATE WRITTEN
!   910620  Modified to meet SLATEC code and prologue standards.  (BKS)
!***END PROLOGUE  DMAKE3
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
      DOUBLE PRECISION ROGUE
      PARAMETER (ROGUE=-1.0D10)
!     .. Scalar Arguments ..
      DOUBLE PRECISION Transl
      INTEGER Lda , M , N , Nmax
      LOGICAL Reset
      CHARACTER*1 Diag , Uplo
      CHARACTER*2 Type
!     .. Array Arguments ..
      DOUBLE PRECISION A(Nmax,*) , Aa(*)
!     .. Local Scalars ..
      INTEGER i , ibeg , iend , j
      LOGICAL gen , lower , sym , tri , unit , upper
!     .. External Functions ..
      DOUBLE PRECISION DBEG
      EXTERNAL DBEG
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN
!***FIRST EXECUTABLE STATEMENT  DMAKE3
      gen = Type=='GE'
      sym = Type=='SY'
      tri = Type=='TR'
      upper = (sym.OR.tri) .AND. Uplo=='U'
      lower = (sym.OR.tri) .AND. Uplo=='L'
      unit = tri .AND. Diag=='U'
!
!     Generate data in array A.
!
      DO j = 1 , N
        DO i = 1 , M
          IF ( gen.OR.(upper.AND.i<=j).OR.(lower.AND.i>=j) ) THEN
            A(i,j) = DBEG(Reset) + Transl
            IF ( i/=j ) THEN
!                 Set some elements to zero
              IF ( N>3.AND.j==N/2 ) A(i,j) = ZERO
              IF ( sym ) THEN
                A(j,i) = A(i,j)
              ELSEIF ( tri ) THEN
                A(j,i) = ZERO
              ENDIF
            ENDIF
          ENDIF
        ENDDO
        IF ( tri ) A(j,j) = A(j,j) + ONE
        IF ( unit ) A(j,j) = ONE
      ENDDO
!
!     Store elements in array AS in data structure required by routine.
!
      IF ( Type=='GE' ) THEN
        DO j = 1 , N
          DO i = 1 , M
            Aa(i+(j-1)*Lda) = A(i,j)
          ENDDO
          DO i = M + 1 , Lda
            Aa(i+(j-1)*Lda) = ROGUE
          ENDDO
        ENDDO
      ELSEIF ( Type=='SY'.OR.Type=='TR' ) THEN
        DO j = 1 , N
          IF ( upper ) THEN
            ibeg = 1
            IF ( unit ) THEN
              iend = j - 1
            ELSE
              iend = j
            ENDIF
          ELSE
            IF ( unit ) THEN
              ibeg = j + 1
            ELSE
              ibeg = j
            ENDIF
            iend = N
          ENDIF
          DO i = 1 , ibeg - 1
            Aa(i+(j-1)*Lda) = ROGUE
          ENDDO
          DO i = ibeg , iend
            Aa(i+(j-1)*Lda) = A(i,j)
          ENDDO
          DO i = iend + 1 , Lda
            Aa(i+(j-1)*Lda) = ROGUE
          ENDDO
        ENDDO
      ENDIF
!
!     End of DMAKE3.
!
      END SUBROUTINE DMAKE3
