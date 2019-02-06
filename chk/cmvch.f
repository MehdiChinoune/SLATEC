*DECK CMVCH
      SUBROUTINE CMVCH (TRANS, M, N, ALPHA, A, NMAX, X, INCX, BETA, Y,
     $   INCY, YT, G, YY, EPS, ERR, FTL, NOUT, MV, KPRINT)
C***BEGIN PROLOGUE  CMVCH
C***SUBSIDIARY
C***PURPOSE  Check the results of the computational tests.
C***LIBRARY   SLATEC (BLAS)
C***AUTHOR  Du Croz, J. J., (NAG)
C           Hanson, R. J., (SNLA)
C***DESCRIPTION
C
C  Checks the results of the computational tests.
C
C  Auxiliary routine for test program for Level 2 Blas.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   870810  DATE WRITTEN
C   910620  Modified to meet SLATEC code and prologue standards.  (BKS)
C***END PROLOGUE  CMVCH
C     .. Parameters ..
      COMPLEX            ZERO
      PARAMETER          ( ZERO = ( 0.0, 0.0 ) )
      REAL               RZERO, RONE
      PARAMETER          ( RZERO = 0.0, RONE = 1.0 )
C     .. Scalar Arguments ..
      COMPLEX            ALPHA, BETA
      REAL               EPS, ERR
      INTEGER            INCX, INCY, KPRINT, M, N, NMAX, NOUT
      LOGICAL            MV, FTL
      CHARACTER*1        TRANS
C     .. Array Arguments ..
      COMPLEX            A( NMAX, * ), X( * ), Y( * ), YT( * ), YY( * )
      REAL               G( * )
C     .. Local Scalars ..
      COMPLEX            C
      REAL               ERRI
      INTEGER            I, INCXL, INCYL, IY, J, JX, K, KX, KY, ML, NL
      LOGICAL            CTRAN, TRAN
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, AIMAG, CONJG, MAX, REAL, SQRT
C     .. Statement Functions ..
      REAL               ABS1
C     .. Statement Function definitions ..
      ABS1( C ) = ABS( REAL( C ) ) + ABS( AIMAG( C ) )
C***FIRST EXECUTABLE STATEMENT  CMVCH
      TRAN = TRANS.EQ.'T'
      CTRAN = TRANS.EQ.'C'
      IF( TRAN.OR.CTRAN )THEN
         ML = N
         NL = M
      ELSE
         ML = M
         NL = N
      END IF
      IF( INCX.LT.0 )THEN
         KX = NL
         INCXL = -1
      ELSE
         KX = 1
         INCXL = 1
      END IF
      IF( INCY.LT.0 )THEN
         KY = ML
         INCYL = -1
      ELSE
         KY = 1
         INCYL = 1
      END IF
C
C     Compute expected result in YT using data in A, X and Y.
C     Compute gauges in G.
C
      IY = KY
      DO 40 I = 1, ML
         YT( IY ) = ZERO
         G( IY ) = RZERO
         JX = KX
         IF( TRAN )THEN
            DO 10 J = 1, NL
               YT( IY ) = YT( IY ) + A( J, I )*X( JX )
               G( IY ) = G( IY ) + ABS1( A( J, I ) )*ABS1( X( JX ) )
               JX = JX + INCXL
   10       CONTINUE
         ELSE IF( CTRAN )THEN
            DO 20 J = 1, NL
               YT( IY ) = YT( IY ) + CONJG( A( J, I ) )*X( JX )
               G( IY ) = G( IY ) + ABS1( A( J, I ) )*ABS1( X( JX ) )
               JX = JX + INCXL
   20       CONTINUE
         ELSE
            DO 30 J = 1, NL
               YT( IY ) = YT( IY ) + A( I, J )*X( JX )
               G( IY ) = G( IY ) + ABS1( A( I, J ) )*ABS1( X( JX ) )
               JX = JX + INCXL
   30       CONTINUE
         END IF
         YT( IY ) = ALPHA*YT( IY ) + BETA*Y( IY )
         G( IY ) = ABS1( ALPHA )*G( IY ) + ABS1( BETA )*ABS1( Y( IY ) )
         IY = IY + INCYL
   40 CONTINUE
C
C     Compute the error ratio for this result.
C
      ERR = ZERO
      DO 50 I = 1, ML
         ERRI = ABS( YT( I ) - YY( 1 + ( I - 1 )*ABS( INCY ) ) )/EPS
         IF( G( I ).NE.RZERO )
     $      ERRI = ERRI/G( I )
         ERR = MAX( ERR, ERRI )
         IF( ERR*SQRT( EPS ).GE.RONE ) THEN
           FTL = .TRUE.
           IF (KPRINT .GE. 2) THEN
              WRITE( NOUT, FMT = 9999 )
              DO 70 K = 1, ML
                 IF( MV )THEN
                    WRITE( NOUT, FMT = 9998 )K, YT( K ),
     $                 YY( 1 + ( K - 1 )*ABS( INCY ) )
                 ELSE
                    WRITE( NOUT, FMT = 9998 )I,
     $                 YY( 1 + ( K - 1 )*ABS( INCY ) ), YT( K )
                 END IF
   70         CONTINUE
            ENDIF
         ENDIF
C
   50 CONTINUE
      RETURN
C
 9999 FORMAT( ' ******* FATAL ERROR - COMPUTED RESULT IS LESS THAN HAL',
     $      'F ACCURATE *******', /'                       EXPECTED RE',
     $      'SULT                    COMPUTED RESULT' )
 9998 FORMAT( 1X, I7, 2( '  (', G15.6, ',', G15.6, ')' ) )
C
C     End of CMVCH.
C
      END
