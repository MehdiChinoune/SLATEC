*DECK DMMCH
      SUBROUTINE DMMCH (TRANSA, TRANSB, M, N, KK, ALPHA, A, LDA, B, LDB,
     $   BETA, C, LDC, CT, G, CC, LDCC, EPS, ERR, FTL, NOUT, MV, KPRINT)
C***BEGIN PROLOGUE  DMMCH
C***SUBSIDIARY
C***PURPOSE  Check the results of the computational tests.
C***LIBRARY   SLATEC (BLAS)
C***AUTHOR  Dongarra, J. J., (ANL)
C           Duff, I., (AERE)
C           Du Croz, J., (NAG)
C           Hammarling, S., (NAG)
C***DESCRIPTION
C
C  Checks the results of the computational tests.
C
C  Auxiliary routine for test program for Level 3 Blas.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   890208  DATE WRITTEN
C   910620  Modified to meet SLATEC code and prologue standards.  (BKS)
C***END PROLOGUE  DMMCH
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      LOGICAL            FTL
      DOUBLE PRECISION   ALPHA, BETA, EPS, ERR
      INTEGER            KK, KPRINT, LDA, LDB, LDC, LDCC, M, N, NOUT
      LOGICAL            MV
      CHARACTER*1        TRANSA, TRANSB
C     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * ),
     $                   CC( LDCC, * ), CT( * ), G( * )
C     .. Local Scalars ..
      DOUBLE PRECISION   ERRI
      INTEGER            I, J, K
      LOGICAL            TRANA, TRANB
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
C***FIRST EXECUTABLE STATEMENT  DMMCH
      TRANA = TRANSA.EQ.'T'.OR.TRANSA.EQ.'C'
      TRANB = TRANSB.EQ.'T'.OR.TRANSB.EQ.'C'
C
C     Compute expected result, one column at a time, in CT using data
C     in A, B and C.
C     Compute gauges in G.
C
      DO 120 J = 1, N
C
         DO 10 I = 1, M
            CT( I ) = ZERO
            G( I ) = ZERO
   10    CONTINUE
         IF( .NOT.TRANA.AND..NOT.TRANB )THEN
            DO 30 K = 1, KK
               DO 20 I = 1, M
                  CT( I ) = CT( I ) + A( I, K )*B( K, J )
                  G( I ) = G( I ) + ABS( A( I, K ) )*ABS( B( K, J ) )
   20          CONTINUE
   30       CONTINUE
         ELSE IF( TRANA.AND..NOT.TRANB )THEN
            DO 50 K = 1, KK
               DO 40 I = 1, M
                  CT( I ) = CT( I ) + A( K, I )*B( K, J )
                  G( I ) = G( I ) + ABS( A( K, I ) )*ABS( B( K, J ) )
   40          CONTINUE
   50       CONTINUE
         ELSE IF( .NOT.TRANA.AND.TRANB )THEN
            DO 70 K = 1, KK
               DO 60 I = 1, M
                  CT( I ) = CT( I ) + A( I, K )*B( J, K )
                  G( I ) = G( I ) + ABS( A( I, K ) )*ABS( B( J, K ) )
   60          CONTINUE
   70       CONTINUE
         ELSE IF( TRANA.AND.TRANB )THEN
            DO 90 K = 1, KK
               DO 80 I = 1, M
                  CT( I ) = CT( I ) + A( K, I )*B( J, K )
                  G( I ) = G( I ) + ABS( A( K, I ) )*ABS( B( J, K ) )
   80          CONTINUE
   90       CONTINUE
         END IF
         DO 100 I = 1, M
            CT( I ) = ALPHA*CT( I ) + BETA*C( I, J )
            G( I ) = ABS( ALPHA )*G( I ) + ABS( BETA )*ABS( C( I, J ) )
  100    CONTINUE
C
C        Compute the error ratio for this result.
C
         ERR = ZERO
         DO 110 I = 1, M
            ERRI = ABS( CT( I ) - CC( I, J ) )/EPS
            IF( G( I ).NE.ZERO )
     $         ERRI = ERRI/G( I )
            ERR = MAX( ERR, ERRI )
            IF( ERR*SQRT( EPS ).GE.ONE ) THEN
            FTL = .TRUE.
             IF (KPRINT .GE. 2) THEN
             WRITE( NOUT, FMT = 9999 )
             DO 140 K = 1, M
                IF( MV )THEN
                   WRITE( NOUT, FMT = 9998 )K, CT( K ), CC( K, J )
                ELSE
                   WRITE( NOUT, FMT = 9998 )K, CC( K, J ), CT( K )
                END IF
  140        CONTINUE
             ENDIF
           ENDIF
  110    CONTINUE
  120 CONTINUE
      RETURN
C
 9999 FORMAT( ' ******* FATAL ERROR - COMPUTED RESULT IS LESS THAN HAL',
     $      'F ACCURATE *******', /'           EXPECTED RESULT   COMPU',
     $      'TED RESULT' )
 9998 FORMAT( 1X, I7, 2G18.6 )
C
C     End of DMMCH.
C
      END
