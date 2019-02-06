*DECK SCHK62
      SUBROUTINE SCHK62 (SNAME, EPS, THRESH, NOUT, KPRINT, FATAL, NIDIM,
     $   IDIM, NALF, ALF, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS,
     $   Y, YY, YS, YT, G, Z)
C***BEGIN PROLOGUE  SCHK62
C***SUBSIDIARY
C***PURPOSE  Quick check for SSYR2 and SSPR2.
C***LIBRARY   SLATEC (BLAS)
C***KEYWORDS  BLAS, QUICK CHECK SERVICE ROUTINE
C***AUTHOR  Du Croz, J. (NAG)
C           Hanson, R. J. (SNLA)
C***DESCRIPTION
C
C  Quick check for SSYR2 and SSPR2.
C
C  Auxiliary routine for test program for Level 2 Blas.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  LSE, LSERES, NUMXER, SMAKE2, SMVCH, SSPR2, SSYR2
C***REVISION HISTORY  (YYMMDD)
C   870810  DATE WRITTEN
C   910619  Modified to meet SLATEC code and prologue standards.  (BKS)
C***END PROLOGUE  SCHK62
C     .. Parameters ..
      REAL               ZERO, HALF, ONE
      PARAMETER          ( ZERO = 0.0, HALF = 0.5, ONE = 1.0 )
C     .. Scalar Arguments ..
      LOGICAL            FATAL
      REAL               EPS, THRESH
      INTEGER            INCMAX, KPRINT, NALF, NIDIM, NINC, NMAX, NOUT
      CHARACTER*6        SNAME
C     .. Array Arguments ..
      REAL               A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ),
     $                   AS( NMAX*NMAX ), G( NMAX ), X( NMAX ),
     $                   XS( NMAX*INCMAX ), XX( NMAX*INCMAX ),
     $                   Y( NMAX ), YS( NMAX*INCMAX ), YT( NMAX ),
     $                   YY( NMAX*INCMAX ), Z( NMAX, 2 )
      INTEGER            IDIM( NIDIM ), INC( NINC )
C     .. Local Scalars ..
      REAL               ALPHA, ALS, ERR, ERRMAX, TRANSL
      INTEGER            I, IA, IC, IN, INCX, INCXS, INCY, INCYS, IX,
     $                   IY, J, JA, JJ, LAA, LDA, LDAS, LJ, LX, LY, N,
     $                   NARGS, NC, NS, NERR
      LOGICAL            FTL, FULL, NULL, PACKED, RESET, UPPER
      CHARACTER*1        UPLO, UPLOS
      CHARACTER*2        ICH
C     .. Local Arrays ..
      REAL               W( 2 )
      LOGICAL            ISAME( 13 )
C     .. External Functions ..
      INTEGER            NUMXER
      LOGICAL            LSE, LSERES
      EXTERNAL           LSE, LSERES, NUMXER
C     .. External Subroutines ..
      EXTERNAL           SMAKE2, SMVCH, SSPR2, SSYR2
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
C     .. Data statements ..
      DATA               ICH/'UL'/
C***FIRST EXECUTABLE STATEMENT  SCHK62
      FULL = SNAME( 3: 3 ).EQ.'Y'
      PACKED = SNAME( 3: 3 ).EQ.'P'
C     Define the number of arguments.
      IF( FULL )THEN
         NARGS = 9
      ELSE IF( PACKED )THEN
         NARGS = 8
      END IF
C
      NC = 0
      RESET = .TRUE.
      ERRMAX = ZERO
C
      DO 140 IN = 1, NIDIM
         N = IDIM( IN )
C        Set LDA to 1 more than minimum value if room.
         LDA = N
         IF( LDA.LT.NMAX )
     $      LDA = LDA + 1
C        Skip tests if not enough room.
         IF( LDA.GT.NMAX )
     $      GO TO 140
         IF( PACKED )THEN
            LAA = ( N*( N + 1 ) )/2
         ELSE
            LAA = LDA*N
         END IF
C
         DO 130 IC = 1, 2
            UPLO = ICH( IC: IC )
            UPPER = UPLO.EQ.'U'
C
            DO 120 IX = 1, NINC
               INCX = INC( IX )
               LX = ABS( INCX )*N
C
C              Generate the vector X.
C
               TRANSL = HALF
               CALL SMAKE2( 'GE', ' ', ' ', 1, N, X, 1, XX, ABS( INCX ),
     $                     0, N - 1, RESET, TRANSL )
               IF( N.GT.1 )THEN
                  X( N/2 ) = ZERO
                  XX( 1 + ABS( INCX )*( N/2 - 1 ) ) = ZERO
               END IF
C
               DO 110 IY = 1, NINC
                  INCY = INC( IY )
                  LY = ABS( INCY )*N
C
C                 Generate the vector Y.
C
                  TRANSL = ZERO
                  CALL SMAKE2( 'GE', ' ', ' ', 1, N, Y, 1, YY,
     $                        ABS( INCY ), 0, N - 1, RESET, TRANSL )
                  IF( N.GT.1 )THEN
                     Y( N/2 ) = ZERO
                     YY( 1 + ABS( INCY )*( N/2 - 1 ) ) = ZERO
                  END IF
C
                  DO 100 IA = 1, NALF
                     ALPHA = ALF( IA )
                     NULL = N.LE.0.OR.ALPHA.EQ.ZERO
C
C                    Generate the matrix A.
C
                     TRANSL = ZERO
                     CALL SMAKE2( SNAME( 2: 3 ), UPLO, ' ', N, N, A,
     $                           NMAX, AA, LDA, N - 1, N - 1, RESET,
     $                           TRANSL )
C
                     NC = NC + 1
C
C                    Save every datum before calling the subroutine.
C
                     UPLOS = UPLO
                     NS = N
                     ALS = ALPHA
                     DO 10 I = 1, LAA
                        AS( I ) = AA( I )
   10                CONTINUE
                     LDAS = LDA
                     DO 20 I = 1, LX
                        XS( I ) = XX( I )
   20                CONTINUE
                     INCXS = INCX
                     DO 30 I = 1, LY
                        YS( I ) = YY( I )
   30                CONTINUE
                     INCYS = INCY
C
C                    Call the subroutine.
C
                     IF( FULL )THEN
                        CALL SSYR2( UPLO, N, ALPHA, XX, INCX, YY, INCY,
     $                              AA, LDA )
                     ELSE IF( PACKED )THEN
                        CALL SSPR2( UPLO, N, ALPHA, XX, INCX, YY, INCY,
     $                              AA )
                     END IF
C
C                    Check if error-exit was taken incorrectly.
C
                     IF (NUMXER(NERR) .NE. 0) THEN
                        IF (KPRINT .GE. 2) THEN
                          WRITE( NOUT, FMT = 9992 )
                        ENDIF
                        FATAL = .TRUE.
                     END IF
C
C                    See what data changed inside subroutines.
C
                     ISAME( 1 ) = UPLO.EQ.UPLOS
                     ISAME( 2 ) = NS.EQ.N
                     ISAME( 3 ) = ALS.EQ.ALPHA
                     ISAME( 4 ) = LSE( XS, XX, LX )
                     ISAME( 5 ) = INCXS.EQ.INCX
                     ISAME( 6 ) = LSE( YS, YY, LY )
                     ISAME( 7 ) = INCYS.EQ.INCY
                     IF( NULL )THEN
                        ISAME( 8 ) = LSE( AS, AA, LAA )
                     ELSE
                        ISAME( 8 ) = LSERES( SNAME( 2: 3 ), UPLO, N, N,
     $                               AS, AA, LDA )
                     END IF
                     IF( .NOT.PACKED )THEN
                        ISAME( 9 ) = LDAS.EQ.LDA
                     END IF
C
C                    If data was incorrectly changed, report and return.
C
                     DO 40 I = 1, NARGS
                       IF (.NOT. ISAME( I )) THEN
                         FATAL = .TRUE.
                         IF (KPRINT .GE. 2) THEN
                           WRITE( NOUT, FMT = 9998 )I
                         ENDIF
                       ENDIF
  40                 CONTINUE
C
                     FTL = .FALSE.
                     IF( .NOT.NULL )THEN
C
C                       Check the result column by column.
C
                        IF( INCX.GT.0 )THEN
                           DO 50 I = 1, N
                              Z( I, 1 ) = X( I )
   50                      CONTINUE
                        ELSE
                           DO 60 I = 1, N
                              Z( I, 1 ) = X( N - I + 1 )
   60                      CONTINUE
                        END IF
                        IF( INCY.GT.0 )THEN
                           DO 70 I = 1, N
                              Z( I, 2 ) = Y( I )
   70                      CONTINUE
                        ELSE
                           DO 80 I = 1, N
                              Z( I, 2 ) = Y( N - I + 1 )
   80                      CONTINUE
                        END IF
                        JA = 1
                        DO 90 J = 1, N
                           W( 1 ) = Z( J, 2 )
                           W( 2 ) = Z( J, 1 )
                           IF( UPPER )THEN
                              JJ = 1
                              LJ = J
                           ELSE
                              JJ = J
                              LJ = N - J + 1
                           END IF
                           CALL SMVCH( 'N', LJ, 2, ALPHA, Z( JJ, 1 ),
     $                                 NMAX, W, 1, ONE, A( JJ, J ), 1,
     $                                 YT, G, AA( JA ), EPS, ERR, FTL,
     $                                 NOUT, .TRUE., KPRINT )
                           IF( FULL )THEN
                              IF( UPPER )THEN
                                 JA = JA + LDA
                              ELSE
                                 JA = JA + LDA + 1
                              END IF
                           ELSE
                              JA = JA + LJ
                           END IF
                           ERRMAX = MAX( ERRMAX, ERR )
                           IF (FTL) THEN
                             FATAL = .TRUE.
                             IF (KPRINT .GE. 3) THEN
                               WRITE( NOUT, FMT = 9995 )J
                               WRITE( NOUT, FMT = 9996 )SNAME
                               IF( FULL )THEN
                                  WRITE( NOUT, FMT = 9993 )NC, SNAME,
     $                             UPLO, N, ALPHA, INCX,
     $                             INCY, LDA
                               ELSE IF( PACKED )THEN
                                  WRITE( NOUT, FMT = 9994 )NC, SNAME,
     $                             UPLO, N, ALPHA, INCX, INCY
                               END IF
                             ENDIF
                           ENDIF
   90                   CONTINUE
                     END IF
C
  100             CONTINUE
C
  110          CONTINUE
C
  120       CONTINUE
C
  130    CONTINUE
C
  140 CONTINUE
C
C     Report result.
C
      IF (.NOT. (FATAL)) THEN
        IF (KPRINT .GE. 3) THEN
          IF( ERRMAX.LT.THRESH )THEN
             WRITE( NOUT, FMT = 9999 )SNAME, NC
          ELSE
             WRITE( NOUT, FMT = 9997 )SNAME, NC, ERRMAX
          END IF
        ENDIF
      ENDIF
      RETURN
C
 9999 FORMAT( ' ', A6, ' PASSED THE COMPUTATIONAL TESTS (', I6, ' CALL',
     $      'S)' )
 9998 FORMAT( ' ******* FATAL ERROR - PARAMETER NUMBER ', I2, ' WAS CH',
     $      'ANGED INCORRECTLY *******' )
 9997 FORMAT( ' ', A6, ' COMPLETED THE COMPUTATIONAL TESTS (', I6, ' C',
     $      'ALLS)', /' ******* BUT WITH MAXIMUM TEST RATIO', F8.2,
     $      ' - SUSPECT *******' )
 9996 FORMAT( ' ******* ', A6, ' FAILED ON CALL NUMBER:' )
 9995 FORMAT( '      THESE ARE THE RESULTS FOR COLUMN ', I3 )
 9994 FORMAT( 1X, I6, ': ', A6, '(''', A1, ''',', I3, ',', F4.1, ', X,',
     $      I2, ', Y,', I2, ', AP)                     .' )
 9993 FORMAT( 1X, I6, ': ', A6, '(''', A1, ''',', I3, ',', F4.1, ', X,',
     $      I2, ', Y,', I2, ', A,', I3, ')                  .' )
 9992 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',
     $      '******' )
C
C     End of SCHK62.
C
      END
