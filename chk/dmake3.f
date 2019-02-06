*DECK DMAKE3
      SUBROUTINE DMAKE3 (TYPE, UPLO, DIAG, M, N, A, NMAX, AA, LDA,
     $   RESET, TRANSL)
C***BEGIN PROLOGUE  DMAKE3
C***SUBSIDIARY
C***PURPOSE  Generate values for an M by N matrix A.
C***LIBRARY   SLATEC (BLAS)
C***AUTHOR  Dongarra, J. J., (ANL)
C           Duff, I., (AERE)
C           Du Croz, J., (NAG)
C           Hammarling, S., (NAG)
C***DESCRIPTION
C
C  Generates values for an M by N matrix A within the bandwidth
C  defined by KL and KU.
C  Stores the values in the array AA in the data structure required
C  by the routine, with unwanted elements set to rogue value.
C
C  TYPE is 'GE', 'SY' or 'TR'.
C
C  Auxiliary routine for test program for Level 3 Blas.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DBEG
C***REVISION HISTORY  (YYMMDD)
C   890208  DATE WRITTEN
C   910620  Modified to meet SLATEC code and prologue standards.  (BKS)
C***END PROLOGUE  DMAKE3
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
      DOUBLE PRECISION   ROGUE
      PARAMETER          ( ROGUE = -1.0D10 )
C     .. Scalar Arguments ..
      DOUBLE PRECISION   TRANSL
      INTEGER            LDA, M, N, NMAX
      LOGICAL            RESET
      CHARACTER*1        DIAG, UPLO
      CHARACTER*2        TYPE
C     .. Array Arguments ..
      DOUBLE PRECISION   A( NMAX, * ), AA( * )
C     .. Local Scalars ..
      INTEGER            I, IBEG, IEND,  J
      LOGICAL            GEN, LOWER, SYM, TRI, UNIT, UPPER
C     .. External Functions ..
      DOUBLE PRECISION   DBEG
      EXTERNAL           DBEG
C     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
C***FIRST EXECUTABLE STATEMENT  DMAKE3
      GEN = TYPE.EQ.'GE'
      SYM = TYPE.EQ.'SY'
      TRI = TYPE.EQ.'TR'
      UPPER = ( SYM.OR.TRI ).AND.UPLO.EQ.'U'
      LOWER = ( SYM.OR.TRI ).AND.UPLO.EQ.'L'
      UNIT = TRI.AND.DIAG.EQ.'U'
C
C     Generate data in array A.
C
      DO 20 J = 1, N
         DO 10 I = 1, M
            IF( GEN.OR.( UPPER.AND.I.LE.J ).OR.( LOWER.AND.I.GE.J ) )
     $          THEN
               A( I, J ) = DBEG( RESET ) + TRANSL
               IF( I.NE.J )THEN
C                 Set some elements to zero
                  IF( N.GT.3.AND.J.EQ.N/2 )
     $               A( I, J ) = ZERO
                  IF( SYM )THEN
                     A( J, I ) = A( I, J )
                  ELSE IF( TRI )THEN
                     A( J, I ) = ZERO
                  END IF
               END IF
            END IF
   10    CONTINUE
         IF( TRI )
     $      A( J, J ) = A( J, J ) + ONE
         IF( UNIT )
     $      A( J, J ) = ONE
   20 CONTINUE
C
C     Store elements in array AS in data structure required by routine.
C
      IF( TYPE.EQ.'GE' )THEN
         DO 50 J = 1, N
            DO 30 I = 1, M
               AA( I + ( J - 1 )*LDA ) = A( I, J )
   30       CONTINUE
            DO 40 I = M + 1, LDA
               AA( I + ( J - 1 )*LDA ) = ROGUE
   40       CONTINUE
   50    CONTINUE
      ELSE IF( TYPE.EQ.'SY'.OR.TYPE.EQ.'TR' )THEN
         DO 90 J = 1, N
            IF( UPPER )THEN
               IBEG = 1
               IF( UNIT )THEN
                  IEND = J - 1
               ELSE
                  IEND = J
               END IF
            ELSE
               IF( UNIT )THEN
                  IBEG = J + 1
               ELSE
                  IBEG = J
               END IF
               IEND = N
            END IF
            DO 60 I = 1, IBEG - 1
               AA( I + ( J - 1 )*LDA ) = ROGUE
   60       CONTINUE
            DO 70 I = IBEG, IEND
               AA( I + ( J - 1 )*LDA ) = A( I, J )
   70       CONTINUE
            DO 80 I = IEND + 1, LDA
               AA( I + ( J - 1 )*LDA ) = ROGUE
   80       CONTINUE
   90    CONTINUE
      END IF
      RETURN
C
C     End of DMAKE3.
C
      END
