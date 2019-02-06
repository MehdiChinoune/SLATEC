!*==DMAKE2.f90  processed by SPAG 6.72Dc at 10:52 on  6 Feb 2019
!DECK DMAKE2
      SUBROUTINE DMAKE2(Type,Uplo,Diag,M,N,A,Nmax,Aa,Lda,Kl,Ku,Reset,Transl)
      IMPLICIT NONE
!*--DMAKE25
!***BEGIN PROLOGUE  DMAKE2
!***SUBSIDIARY
!***PURPOSE  Generate values for an M by N matrix A.
!***LIBRARY   SLATEC (BLAS)
!***AUTHOR  Du Croz, J. J., (NAG)
!           Hanson, R. J., (SNLA)
!***DESCRIPTION
!
!  Generates values for an M by N matrix A within the bandwidth
!  defined by KL and KU.
!  Stores the values in the array AA in the data structure required
!  by the routine, with unwanted elements set to rogue value.
!
!  TYPE is 'GE', 'GB', 'SY', 'SB', 'SP', 'TR', 'TB' OR 'TP'.
!
!  Auxiliary routine for test program for Level 2 Blas.
!***REFERENCES  (NONE)
!***ROUTINES CALLED  DBEG
!***REVISION HISTORY  (YYMMDD)
!   870810  DATE WRITTEN
!   910620  Modified to meet SLATEC code and prologue standards.  (BKS)
!***END PROLOGUE  DMAKE2
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
      DOUBLE PRECISION ROGUE
      PARAMETER (ROGUE=-1.0D10)
!     .. Scalar Arguments ..
      DOUBLE PRECISION Transl
      INTEGER Kl , Ku , Lda , M , N , Nmax
      LOGICAL Reset
      CHARACTER*1 Diag , Uplo
      CHARACTER*2 Type
!     .. Array Arguments ..
      DOUBLE PRECISION A(Nmax,*) , Aa(*)
!     .. Local Scalars ..
      INTEGER i , i1 , i2 , i3 , ibeg , iend , ioff , j , kk
      LOGICAL gen , lower , sym , tri , unit , upper
!     .. External Functions ..
      DOUBLE PRECISION DBEG
      EXTERNAL DBEG
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN
!***FIRST EXECUTABLE STATEMENT  DMAKE2
      gen = Type(1:1)=='G'
      sym = Type(1:1)=='S'
      tri = Type(1:1)=='T'
      upper = (sym.OR.tri) .AND. Uplo=='U'
      lower = (sym.OR.tri) .AND. Uplo=='L'
      unit = tri .AND. Diag=='U'
!
!     Generate data in array A.
!
      DO j = 1 , N
        DO i = 1 , M
          IF ( gen.OR.(upper.AND.i<=j).OR.(lower.AND.i>=j) ) THEN
            IF ( (i<=j.AND.j-i<=Ku).OR.(i>=j.AND.i-j<=Kl) ) THEN
              A(i,j) = DBEG(Reset) + Transl
            ELSE
              A(i,j) = ZERO
            ENDIF
            IF ( i/=j ) THEN
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
      ELSEIF ( Type=='GB' ) THEN
        DO j = 1 , N
          DO i1 = 1 , Ku + 1 - j
            Aa(i1+(j-1)*Lda) = ROGUE
          ENDDO
          DO i2 = i1 , MIN(Kl+Ku+1,Ku+1+M-j)
            Aa(i2+(j-1)*Lda) = A(i2+j-Ku-1,j)
          ENDDO
          DO i3 = i2 , Lda
            Aa(i3+(j-1)*Lda) = ROGUE
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
      ELSEIF ( Type=='SB'.OR.Type=='TB' ) THEN
        DO j = 1 , N
          IF ( upper ) THEN
            kk = Kl + 1
            ibeg = MAX(1,Kl+2-j)
            IF ( unit ) THEN
              iend = Kl
            ELSE
              iend = Kl + 1
            ENDIF
          ELSE
            kk = 1
            IF ( unit ) THEN
              ibeg = 2
            ELSE
              ibeg = 1
            ENDIF
            iend = MIN(Kl+1,1+M-j)
          ENDIF
          DO i = 1 , ibeg - 1
            Aa(i+(j-1)*Lda) = ROGUE
          ENDDO
          DO i = ibeg , iend
            Aa(i+(j-1)*Lda) = A(i+j-kk,j)
          ENDDO
          DO i = iend + 1 , Lda
            Aa(i+(j-1)*Lda) = ROGUE
          ENDDO
        ENDDO
      ELSEIF ( Type=='SP'.OR.Type=='TP' ) THEN
        ioff = 0
        DO j = 1 , N
          IF ( upper ) THEN
            ibeg = 1
            iend = j
          ELSE
            ibeg = j
            iend = N
          ENDIF
          DO i = ibeg , iend
            ioff = ioff + 1
            Aa(ioff) = A(i,j)
            IF ( i==j ) THEN
              IF ( unit ) Aa(ioff) = ROGUE
            ENDIF
          ENDDO
        ENDDO
      ENDIF
!
!     End of DMAKE2.
!
      END SUBROUTINE DMAKE2
