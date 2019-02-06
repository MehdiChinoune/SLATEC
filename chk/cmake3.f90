!*==CMAKE3.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK CMAKE3
      SUBROUTINE CMAKE3(Type,Uplo,Diag,M,N,A,Nmax,Aa,Lda,Reset,Transl)
      IMPLICIT NONE
!*--CMAKE35
!***BEGIN PROLOGUE  CMAKE3
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
!  TYPE is 'GE', 'HE', 'SY', OR 'TR'.
!
!  Auxiliary routine for test program for Level 3 Blas.
!***REFERENCES  (NONE)
!***ROUTINES CALLED  CBEG
!***REVISION HISTORY  (YYMMDD)
!   890208  DATE WRITTEN
!   910620  Modified to meet SLATEC code and prologue standards.  (BKS)
!***END PROLOGUE  CMAKE3
!     .. Parameters ..
      COMPLEX ZERO , ONE
      PARAMETER (ZERO=(0.0,0.0),ONE=(1.0,0.0))
      COMPLEX ROGUE
      PARAMETER (ROGUE=(-1.0E10,1.0E10))
      REAL RZERO
      PARAMETER (RZERO=0.0)
      REAL RROGUE
      PARAMETER (RROGUE=-1.0E10)
!     .. Scalar Arguments ..
      COMPLEX Transl
      INTEGER Lda , M , N , Nmax
      LOGICAL Reset
      CHARACTER*1 Diag , Uplo
      CHARACTER*2 Type
!     .. Array Arguments ..
      COMPLEX A(Nmax,*) , Aa(*)
!     .. Local Scalars ..
      INTEGER i , ibeg , iend , j , jj
      LOGICAL gen , lower , sym , tri , unit , upper , her
!     .. External Functions ..
      COMPLEX CBEG
      EXTERNAL CBEG
!     .. Intrinsic Functions ..
      INTRINSIC CMPLX , CONJG , REAL
!***FIRST EXECUTABLE STATEMENT  CMAKE3
      gen = Type=='GE'
      her = Type=='HE'
      sym = Type=='SY'
      tri = Type=='TR'
      upper = (her.OR.sym.OR.tri) .AND. Uplo=='U'
      lower = (her.OR.sym.OR.tri) .AND. Uplo=='L'
      unit = tri .AND. Diag=='U'
!
!     Generate data in array A.
!
      DO j = 1 , N
        DO i = 1 , M
          IF ( gen.OR.(upper.AND.i<=j).OR.(lower.AND.i>=j) ) THEN
            A(i,j) = CBEG(Reset) + Transl
            IF ( i/=j ) THEN
!                 Set some elements to zero
              IF ( N>3.AND.j==N/2 ) A(i,j) = ZERO
              IF ( her ) THEN
                A(j,i) = CONJG(A(i,j))
              ELSEIF ( sym ) THEN
                A(j,i) = A(i,j)
              ELSEIF ( tri ) THEN
                A(j,i) = ZERO
              ENDIF
            ENDIF
          ENDIF
        ENDDO
        IF ( her ) A(j,j) = CMPLX(REAL(A(j,j)),RZERO)
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
      ELSEIF ( Type=='HE'.OR.Type=='SY'.OR.Type=='TR' ) THEN
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
          IF ( her ) THEN
            jj = j + (j-1)*Lda
            Aa(jj) = CMPLX(REAL(Aa(jj)),RROGUE)
          ENDIF
        ENDDO
      ENDIF
!
!     End of CMAKE3.
!
      END SUBROUTINE CMAKE3
