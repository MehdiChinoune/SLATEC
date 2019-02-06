*DECK SCHK12
      SUBROUTINE SCHK12 (SNAME, EPS, THRESH, NOUT, KPRINT, FATAL, NIDIM,
     $   IDIM, NKB, KB, NALF, ALF, NBET, BET, NINC, INC, NMAX, INCMAX,
     $   A, AA, AS, X, XX, XS, Y, YY, YS, YT, G)
C***BEGIN PROLOGUE  SCHK12
C***SUBSIDIARY
C***PURPOSE  Quick check for SGEMV and SGBMV.
C***LIBRARY   SLATEC (BLAS)
C***KEYWORDS  BLAS, QUICK CHECK SERVICE ROUTINE
C***AUTHOR  Du Croz, J. (NAG)
C           Hanson, R. J. (SNLA)
C***DESCRIPTION
C
C  Quick check for SGEMV and SGBMV.
C
C  Auxiliary routine for test program for Level 2 Blas.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  LSE, LSERES, NUMXER, SGBMV, SGEMV, SMAKE2, SMVCH
C***REVISION HISTORY  (YYMMDD)
C   870810  DATE WRITTEN
C   910619  Modified to meet SLATEC code and prologue standards.  (BKS)
C***END PROLOGUE  SCHK12
C     .. Parameters ..
      REAL               ZERO, HALF
      PARAMETER          ( ZERO = 0.0, HALF = 0.5 )
C     .. Scalar Arguments ..
      LOGICAL            FATAL
      REAL               EPS, THRESH
      INTEGER            INCMAX, KPRINT, NALF, NBET, NIDIM, NINC, NKB,
     $                   NMAX, NOUT
      CHARACTER*6        SNAME
C     .. Array Arguments ..
      REAL               A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ),
     $                   AS( NMAX*NMAX ), BET( NBET ), G( NMAX ),
     $                   X( NMAX ), XS( NMAX*INCMAX ),
     $                   XX( NMAX*INCMAX ), Y( NMAX ),
     $                   YS( NMAX*INCMAX ), YT( NMAX ),
     $                   YY( NMAX*INCMAX )
      INTEGER            IDIM( NIDIM ), INC( NINC ), KB( NKB )
C     .. Local Scalars ..
      REAL               ALPHA, ALS, BETA, BLS, ERR, ERRMAX, TRANSL
      INTEGER            I, IA, IB, IC, IKU, IM, IN, INCX, INCXS, INCY,
     $                   INCYS, IX, IY, KL, KLS, KU, KUS, LAA, LDA,
     $                   LDAS, LX, LY, M, ML, MS, N, NARGS, NC, ND, NK,
     $                   NL, NS, NERR
      LOGICAL            BANDED, FTL, FULL, NULL, RESET, TRAN
      CHARACTER*1        TRANS, TRANSS
      CHARACTER*3        ICH
C     .. Local Arrays ..
      LOGICAL            ISAME( 13 )
C     .. External Functions ..
      INTEGER            NUMXER
      LOGICAL            LSE, LSERES
      EXTERNAL           LSE, LSERES, NUMXER
C     .. External Subroutines ..
      EXTERNAL           SGBMV, SGEMV, SMAKE2, SMVCH
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
C     .. Data statements ..
      DATA               ICH/'NTC'/
C***FIRST EXECUTABLE STATEMENT  SCHK12
      FULL = SNAME( 3: 3 ).EQ.'E'
      BANDED = SNAME( 3: 3 ).EQ.'B'
C     Define the number of arguments.
      IF( FULL )THEN
         NARGS = 11
      ELSE IF( BANDED )THEN
         NARGS = 13
      END IF
C
      NC = 0
      RESET = .TRUE.
      ERRMAX = ZERO
C
      DO 120 IN = 1, NIDIM
         N = IDIM( IN )
         ND = N/2 + 1
C
         DO 110 IM = 1, 2
            IF( IM.EQ.1 )
     $         M = MAX( N - ND, 0 )
            IF( IM.EQ.2 )
     $         M = MIN( N + ND, NMAX )
C
            IF( BANDED )THEN
               NK = NKB
            ELSE
               NK = 1
            END IF
            DO 100 IKU = 1, NK
               IF( BANDED )THEN
                  KU = KB( IKU )
                  KL = MAX( KU - 1, 0 )
               ELSE
                  KU = N - 1
                  KL = M - 1
               END IF
C              Set LDA to 1 more than minimum value if room.
               IF( BANDED )THEN
                  LDA = KL + KU + 1
               ELSE
                  LDA = M
               END IF
               IF( LDA.LT.NMAX )
     $            LDA = LDA + 1
C              Skip tests if not enough room.
               IF( LDA.GT.NMAX )
     $            GO TO 100
               LAA = LDA*N
               NULL = N.LE.0.OR.M.LE.0
C
C              Generate the matrix A.
C
               TRANSL = ZERO
               CALL SMAKE2( SNAME( 2: 3 ), ' ', ' ', M, N, A, NMAX, AA,
     $                     LDA, KL, KU, RESET, TRANSL )
C
               DO 90 IC = 1, 3
                  TRANS = ICH( IC: IC )
                  TRAN = TRANS.EQ.'T'.OR.TRANS.EQ.'C'
C
                  IF( TRAN )THEN
                     ML = N
                     NL = M
                  ELSE
                     ML = M
                     NL = N
                  END IF
C
                  DO 80 IX = 1, NINC
                     INCX = INC( IX )
                     LX = ABS( INCX )*NL
C
C                    Generate the vector X.
C
                     TRANSL = HALF
                     CALL SMAKE2( 'GE', ' ', ' ', 1, NL, X, 1, XX,
     $                           ABS( INCX ), 0, NL - 1, RESET, TRANSL )
                     IF( NL.GT.1 )THEN
                        X( NL/2 ) = ZERO
                        XX( 1 + ABS( INCX )*( NL/2 - 1 ) ) = ZERO
                     END IF
C
                     DO 70 IY = 1, NINC
                        INCY = INC( IY )
                        LY = ABS( INCY )*ML
C
                        DO 60 IA = 1, NALF
                           ALPHA = ALF( IA )
C
                           DO 50 IB = 1, NBET
                              BETA = BET( IB )
C
C                             Generate the vector Y.
C
                              TRANSL = ZERO
                              CALL SMAKE2( 'GE', ' ', ' ', 1, ML, Y, 1,
     $                                    YY, ABS( INCY ), 0, ML - 1,
     $                                    RESET, TRANSL )
C
                              NC = NC + 1
C
C                             Save every datum before calling the
C                             subroutine.
C
                              TRANSS = TRANS
                              MS = M
                              NS = N
                              KLS = KL
                              KUS = KU
                              ALS = ALPHA
                              DO 10 I = 1, LAA
                                 AS( I ) = AA( I )
   10                         CONTINUE
                              LDAS = LDA
                              DO 20 I = 1, LX
                                 XS( I ) = XX( I )
   20                         CONTINUE
                              INCXS = INCX
                              BLS = BETA
                              DO 30 I = 1, LY
                                 YS( I ) = YY( I )
   30                         CONTINUE
                              INCYS = INCY
C
C                             Call the subroutine.
C
                              IF( FULL )THEN
                                 CALL SGEMV( TRANS, M, N, ALPHA, AA,
     $                                       LDA, XX, INCX, BETA, YY,
     $                                       INCY )
                              ELSE IF( BANDED )THEN
                                 CALL SGBMV( TRANS, M, N, KL, KU, ALPHA,
     $                                       AA, LDA, XX, INCX, BETA,
     $                                       YY, INCY )
                              END IF
C
C                             Check if error-exit was taken incorrectly.
C
                              IF (NUMXER(NERR) .NE. 0) THEN
                                 IF (KPRINT .GE. 2) THEN
                                   WRITE( NOUT, FMT = 9993 )
                                 ENDIF
                                 FATAL = .TRUE.
                              END IF
C
C                             See what data changed inside subroutines.
C
                              ISAME( 1 ) = TRANS.EQ.TRANSS
                              ISAME( 2 ) = MS.EQ.M
                              ISAME( 3 ) = NS.EQ.N
                              IF( FULL )THEN
                                 ISAME( 4 ) = ALS.EQ.ALPHA
                                 ISAME( 5 ) = LSE( AS, AA, LAA )
                                 ISAME( 6 ) = LDAS.EQ.LDA
                                 ISAME( 7 ) = LSE( XS, XX, LX )
                                 ISAME( 8 ) = INCXS.EQ.INCX
                                 ISAME( 9 ) = BLS.EQ.BETA
                                 IF( NULL )THEN
                                    ISAME( 10 ) = LSE( YS, YY, LY )
                                 ELSE
                                    ISAME( 10 ) = LSERES( 'GE', ' ', 1,
     $                                            ML, YS, YY,
     $                                            ABS( INCY ) )
                                 END IF
                                 ISAME( 11 ) = INCYS.EQ.INCY
                              ELSE IF( BANDED )THEN
                                 ISAME( 4 ) = KLS.EQ.KL
                                 ISAME( 5 ) = KUS.EQ.KU
                                 ISAME( 6 ) = ALS.EQ.ALPHA
                                 ISAME( 7 ) = LSE( AS, AA, LAA )
                                 ISAME( 8 ) = LDAS.EQ.LDA
                                 ISAME( 9 ) = LSE( XS, XX, LX )
                                 ISAME( 10 ) = INCXS.EQ.INCX
                                 ISAME( 11 ) = BLS.EQ.BETA
                                 IF( NULL )THEN
                                    ISAME( 12 ) = LSE( YS, YY, LY )
                                 ELSE
                                    ISAME( 12 ) = LSERES( 'GE', ' ', 1,
     $                                            ML, YS, YY,
     $                                            ABS( INCY ) )
                                 END IF
                                 ISAME( 13 ) = INCYS.EQ.INCY
                              END IF
C
C                             If data was incorrectly changed, report
C                             and return.
C
                              DO 40 I = 1, NARGS
                                IF (.NOT. ISAME( I )) THEN
                                  FATAL = .TRUE.
                                  IF (KPRINT .GE. 2) THEN
                                    WRITE( NOUT, FMT = 9998 )I
                                  ENDIF
                                ENDIF
  40                          CONTINUE
C
                              FTL = .FALSE.
                              IF( .NOT.NULL )THEN
C
C                                Check the result.
C
                                 CALL SMVCH( TRANS, M, N, ALPHA, A,
     $                                       NMAX, X, INCX, BETA, Y,
     $                                       INCY, YT, G, YY, EPS, ERR,
     $                                       FTL, NOUT, .TRUE.,
     $                                       KPRINT )
                                 ERRMAX = MAX( ERRMAX, ERR )
                              END IF
                              IF (FTL) THEN
                                FATAL = .TRUE.
                                IF (KPRINT .GE. 3) THEN
                                  WRITE( NOUT, FMT = 9996 )SNAME
                                  IF( FULL )THEN
                                     WRITE( NOUT, FMT = 9994 )NC, SNAME,
     $                                TRANS, M, N, ALPHA,
     $                                LDA, INCX, BETA, INCY
                                  ELSE IF( BANDED )THEN
                                     WRITE( NOUT, FMT = 9995 )NC, SNAME,
     $                                TRANS, M, N, KL, KU,
     $                                ALPHA, LDA, INCX, BETA, INCY
                                  END IF
                                ENDIF
                              ENDIF
C
   50                      CONTINUE
C
   60                   CONTINUE
C
   70                CONTINUE
C
   80             CONTINUE
C
   90          CONTINUE
C
  100       CONTINUE
C
  110    CONTINUE
C
  120 CONTINUE
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
 9995 FORMAT( 1X, I6, ': ', A6, '(''', A1, ''',', 4( I3, ',' ), F4.1,
     $      ', A,', I3, ', X,', I2, ',', F4.1, ', Y,', I2, ') .' )
 9994 FORMAT( 1X, I6, ': ', A6, '(''', A1, ''',', 2( I3, ',' ), F4.1,
     $      ', A,', I3, ', X,', I2, ',', F4.1, ', Y,', I2,
     $      ')         .' )
 9993 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',
     $      '******' )
C
C     End of SCHK12.
C
      END
