*DECK CMAKE2
      SUBROUTINE CMAKE2 (TYPE, UPLO, DIAG, M, N, A, NMAX, AA, LDA, KL,
     $   KU, RESET, TRANSL)
C***BEGIN PROLOGUE  CMAKE2
C***SUBSIDIARY
C***PURPOSE  Generate values for an M by N matrix A.
C***LIBRARY   SLATEC (BLAS)
C***AUTHOR  Du Croz, J. J., (NAG)
C           Hanson, R. J., (SNLA)
C***DESCRIPTION
C
C  Generates values for an M by N matrix A within the bandwidth
C  defined by KL and KU.
C  Stores the values in the array AA in the data structure required
C  by the routine, with unwanted elements set to rogue value.
C
C  TYPE is 'GE', 'GB', 'SY', 'SB', 'SP', 'TR', 'TB' OR 'TP'.
C
C  Auxiliary routine for test program for Level 2 Blas.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CBEG
C***REVISION HISTORY  (YYMMDD)
C   870810  DATE WRITTEN
C   910620  Modified to meet SLATEC code and prologue standards.  (BKS)
C***END PROLOGUE  CMAKE2
C     .. Parameters ..
      COMPLEX            ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0, 0.0 ), ONE = ( 1.0, 0.0 ) )
      COMPLEX            ROGUE
      PARAMETER          ( ROGUE = ( -1.0E10, 1.0E10 ) )
      REAL               RZERO
      PARAMETER          ( RZERO = 0.0 )
      REAL               RROGUE
      PARAMETER          ( RROGUE = -1.0E10 )
C     .. Scalar Arguments ..
      COMPLEX            TRANSL
      INTEGER            KL, KU, LDA, M, N, NMAX
      LOGICAL            RESET
      CHARACTER*1        DIAG, UPLO
      CHARACTER*2        TYPE
C     .. Array Arguments ..
      COMPLEX            A( NMAX, * ), AA( * )
C     .. Local Scalars ..
      INTEGER            I, I1, I2, I3, IBEG, IEND, IOFF, J, JJ, KK
      LOGICAL            GEN, LOWER, SYM, TRI, UNIT, UPPER
C     .. External Functions ..
      COMPLEX            CBEG
      EXTERNAL           CBEG
C     .. Intrinsic Functions ..
      INTRINSIC          CMPLX, CONJG, MAX, MIN, REAL
C***FIRST EXECUTABLE STATEMENT  CMAKE2
      GEN = TYPE( 1: 1 ).EQ.'G'
      SYM = TYPE( 1: 1 ).EQ.'H'
      TRI = TYPE( 1: 1 ).EQ.'T'
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
               IF( ( I.LE.J.AND.J - I.LE.KU ).OR.
     $             ( I.GE.J.AND.I - J.LE.KL ) )THEN
                  A( I, J ) = CBEG( RESET ) + TRANSL
               ELSE
                  A( I, J ) = ZERO
               END IF
               IF( I.NE.J )THEN
                  IF( SYM )THEN
                     A( J, I ) = CONJG( A( I, J ) )
                  ELSE IF( TRI )THEN
                     A( J, I ) = ZERO
                  END IF
               END IF
            END IF
   10    CONTINUE
         IF( SYM )
     $      A( J, J ) = CMPLX( REAL( A( J, J ) ), RZERO )
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
      ELSE IF( TYPE.EQ.'GB' )THEN
         DO 90 J = 1, N
            DO 60 I1 = 1, KU + 1 - J
               AA( I1 + ( J - 1 )*LDA ) = ROGUE
   60       CONTINUE
            DO 70 I2 = I1, MIN( KL + KU + 1, KU + 1 + M - J )
               AA( I2 + ( J - 1 )*LDA ) = A( I2 + J - KU - 1, J )
   70       CONTINUE
            DO 80 I3 = I2, LDA
               AA( I3 + ( J - 1 )*LDA ) = ROGUE
   80       CONTINUE
   90    CONTINUE
      ELSE IF( TYPE.EQ.'HE'.OR.TYPE.EQ.'TR' )THEN
         DO 130 J = 1, N
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
            DO 100 I = 1, IBEG - 1
               AA( I + ( J - 1 )*LDA ) = ROGUE
  100       CONTINUE
            DO 110 I = IBEG, IEND
               AA( I + ( J - 1 )*LDA ) = A( I, J )
  110       CONTINUE
            DO 120 I = IEND + 1, LDA
               AA( I + ( J - 1 )*LDA ) = ROGUE
  120       CONTINUE
            IF( SYM )THEN
               JJ = J + ( J - 1 )*LDA
               AA( JJ ) = CMPLX( REAL( AA( JJ ) ), RROGUE )
            END IF
  130    CONTINUE
      ELSE IF( TYPE.EQ.'HB'.OR.TYPE.EQ.'TB' )THEN
         DO 170 J = 1, N
            IF( UPPER )THEN
               KK = KL + 1
               IBEG = MAX( 1, KL + 2 - J )
               IF( UNIT )THEN
                  IEND = KL
               ELSE
                  IEND = KL + 1
               END IF
            ELSE
               KK = 1
               IF( UNIT )THEN
                  IBEG = 2
               ELSE
                  IBEG = 1
               END IF
               IEND = MIN( KL + 1, 1 + M - J )
            END IF
            DO 140 I = 1, IBEG - 1
               AA( I + ( J - 1 )*LDA ) = ROGUE
  140       CONTINUE
            DO 150 I = IBEG, IEND
               AA( I + ( J - 1 )*LDA ) = A( I + J - KK, J )
  150       CONTINUE
            DO 160 I = IEND + 1, LDA
               AA( I + ( J - 1 )*LDA ) = ROGUE
  160       CONTINUE
            IF( SYM )THEN
               JJ = KK + ( J - 1 )*LDA
               AA( JJ ) = CMPLX( REAL( AA( JJ ) ), RROGUE )
            END IF
  170    CONTINUE
      ELSE IF( TYPE.EQ.'HP'.OR.TYPE.EQ.'TP' )THEN
         IOFF = 0
         DO 190 J = 1, N
            IF( UPPER )THEN
               IBEG = 1
               IEND = J
            ELSE
               IBEG = J
               IEND = N
            END IF
            DO 180 I = IBEG, IEND
               IOFF = IOFF + 1
               AA( IOFF ) = A( I, J )
               IF( I.EQ.J )THEN
                  IF( UNIT )
     $               AA( IOFF ) = ROGUE
                  IF( SYM )
     $               AA( IOFF ) = CMPLX( REAL( AA( IOFF ) ), RROGUE )
               END IF
  180       CONTINUE
  190    CONTINUE
      END IF
      RETURN
C
C     End of CMAKE2.
C
      END
