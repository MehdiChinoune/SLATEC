*DECK CMMCH
      SUBROUTINE CMMCH (TRANSA, TRANSB, M, N, KK, ALPHA, A, LDA, B, LDB,
     $   BETA, C, LDC, CT, G, CC, LDCC, EPS, ERR, FTL, NOUT, MV, KPRINT)
C***BEGIN PROLOGUE  CMMCH
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
C***END PROLOGUE  CMMCH
C     .. Parameters ..
      COMPLEX            ZERO
      PARAMETER          ( ZERO = ( 0.0, 0.0 ) )
      REAL               RZERO, RONE
      PARAMETER          ( RZERO = 0.0, RONE = 1.0 )
C     .. Scalar Arguments ..
      LOGICAL            FTL
      COMPLEX            ALPHA, BETA
      REAL               EPS, ERR
      INTEGER            KK, KPRINT, LDA, LDB, LDC, LDCC, M, N, NOUT
      LOGICAL            MV
      CHARACTER*1        TRANSA, TRANSB
C     .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * ), C( LDC, * ),
     $                   CC( LDCC, * ), CT( * )
      REAL               G( * )
C     .. Local Scalars ..
      COMPLEX            CL
      REAL               ERRI
      INTEGER            I, J, K
      LOGICAL            CTRANA, CTRANB, TRANA, TRANB
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, AIMAG, CONJG, MAX, REAL, SQRT
C     .. Statement Functions ..
      REAL               ABS1
C     .. Statement Function definitions ..
      ABS1( CL ) = ABS( REAL( CL ) ) + ABS( AIMAG( CL ) )
C***FIRST EXECUTABLE STATEMENT  CMMCH
      TRANA = TRANSA.EQ.'T'.OR.TRANSA.EQ.'C'
      TRANB = TRANSB.EQ.'T'.OR.TRANSB.EQ.'C'
      CTRANA = TRANSA.EQ.'C'
      CTRANB = TRANSB.EQ.'C'
C
C     Compute expected result, one column at a time, in CT using data
C     in A, B and C.
C     Compute gauges in G.
C
      DO 220 J = 1, N
C
         DO 10 I = 1, M
            CT( I ) = ZERO
            G( I ) = RZERO
   10    CONTINUE
         IF( .NOT.TRANA.AND..NOT.TRANB )THEN
            DO 30 K = 1, KK
               DO 20 I = 1, M
                  CT( I ) = CT( I ) + A( I, K )*B( K, J )
                  G( I ) = G( I ) + ABS1( A( I, K ) )*ABS1( B( K, J ) )
   20          CONTINUE
   30       CONTINUE
         ELSE IF( TRANA.AND..NOT.TRANB )THEN
            IF( CTRANA )THEN
               DO 50 K = 1, KK
                  DO 40 I = 1, M
                     CT( I ) = CT( I ) + CONJG( A( K, I ) )*B( K, J )
                     G( I ) = G( I ) + ABS1( A( K, I ) )*
     $                        ABS1( B( K, J ) )
   40             CONTINUE
   50          CONTINUE
            ELSE
               DO 70 K = 1, KK
                  DO 60 I = 1, M
                     CT( I ) = CT( I ) + A( K, I )*B( K, J )
                     G( I ) = G( I ) + ABS1( A( K, I ) )*
     $                        ABS1( B( K, J ) )
   60             CONTINUE
   70          CONTINUE
            END IF
         ELSE IF( .NOT.TRANA.AND.TRANB )THEN
            IF( CTRANB )THEN
               DO 90 K = 1, KK
                  DO 80 I = 1, M
                     CT( I ) = CT( I ) + A( I, K )*CONJG( B( J, K ) )
                     G( I ) = G( I ) + ABS1( A( I, K ) )*
     $                        ABS1( B( J, K ) )
   80             CONTINUE
   90          CONTINUE
            ELSE
               DO 110 K = 1, KK
                  DO 100 I = 1, M
                     CT( I ) = CT( I ) + A( I, K )*B( J, K )
                     G( I ) = G( I ) + ABS1( A( I, K ) )*
     $                        ABS1( B( J, K ) )
  100             CONTINUE
  110          CONTINUE
            END IF
         ELSE IF( TRANA.AND.TRANB )THEN
            IF( CTRANA )THEN
               IF( CTRANB )THEN
                  DO 130 K = 1, KK
                     DO 120 I = 1, M
                        CT( I ) = CT( I ) + CONJG( A( K, I ) )*
     $                            CONJG( B( J, K ) )
                        G( I ) = G( I ) + ABS1( A( K, I ) )*
     $                           ABS1( B( J, K ) )
  120                CONTINUE
  130             CONTINUE
               ELSE
                  DO 150 K = 1, KK
                     DO 140 I = 1, M
                        CT( I ) = CT( I ) + CONJG( A( K, I ) )*B( J, K )
                        G( I ) = G( I ) + ABS1( A( K, I ) )*
     $                           ABS1( B( J, K ) )
  140                CONTINUE
  150             CONTINUE
               END IF
            ELSE
               IF( CTRANB )THEN
                  DO 170 K = 1, KK
                     DO 160 I = 1, M
                        CT( I ) = CT( I ) + A( K, I )*CONJG( B( J, K ) )
                        G( I ) = G( I ) + ABS1( A( K, I ) )*
     $                           ABS1( B( J, K ) )
  160                CONTINUE
  170             CONTINUE
               ELSE
                  DO 190 K = 1, KK
                     DO 180 I = 1, M
                        CT( I ) = CT( I ) + A( K, I )*B( J, K )
                        G( I ) = G( I ) + ABS1( A( K, I ) )*
     $                           ABS1( B( J, K ) )
  180                CONTINUE
  190             CONTINUE
               END IF
            END IF
         END IF
         DO 200 I = 1, M
            CT( I ) = ALPHA*CT( I ) + BETA*C( I, J )
            G( I ) = ABS1( ALPHA )*G( I ) +
     $               ABS1( BETA )*ABS1( C( I, J ) )
  200    CONTINUE
C
C        Compute the error ratio for this result.
C
         ERR = ZERO
         DO 210 I = 1, M
            ERRI = ABS1( CT( I ) - CC( I, J ) )/EPS
            IF( G( I ).NE.RZERO )
     $         ERRI = ERRI/G( I )
            ERR = MAX( ERR, ERRI )
            IF( ERR*SQRT( EPS ).GE.RONE ) THEN
               FTL = .TRUE.
               IF (KPRINT .GE. 2) THEN
                 WRITE( NOUT, FMT = 9999 )
                 DO 240 K = 1, M
                 IF( MV )THEN
                   WRITE( NOUT, FMT = 9998 )K, CT( K ), CC( K, J )
                 ELSE
                   WRITE( NOUT, FMT = 9998 )K, CC( K, J ), CT( K )
                 END IF
  240            CONTINUE
               ENDIF
          ENDIF
  210   CONTINUE
  220   CONTINUE
      RETURN
C
 9999 FORMAT( ' ******* FATAL ERROR - COMPUTED RESULT IS LESS THAN HAL',
     $      'F ACCURATE *******', /'                       EXPECTED RE',
     $      'SULT                    COMPUTED RESULT' )
 9998 FORMAT( 1X, I7, 2( '  (', G15.6, ',', G15.6, ')' ) )
C
C     End of CMMCH.
C
      END
