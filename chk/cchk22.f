*DECK CCHK22
      SUBROUTINE CCHK22 (SNAME, EPS, THRESH, NOUT, KPRINT, FATAL, NIDIM,
     $   IDIM, NKB, KB, NALF, ALF, NBET, BET, NINC, INC, NMAX, INCMAX,
     $   A, AA, AS, X, XX, XS, Y, YY, YS, YT, G)
C***BEGIN PROLOGUE  CCHK22
C***SUBSIDIARY
C***PURPOSE  Quick check for CHEMV, CHBMV, CHPMV.
C***LIBRARY   SLATEC (BLAS)
C***KEYWORDS  BLAS, QUICK CHECK SERVICE ROUTINE
C***AUTHOR  Du Croz, J. (NAG)
C           Hanson, R. J. (SNLA)
C***DESCRIPTION
C
C  Quick check for CHEMV, CHBMV and CHPMV.
C
C  Auxiliary routine for test program for Level 2 Blas.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CHBMV, CHEMV, CHPMV, CMAKE2, CMVCH, LCE, LCERES,
C                    NUMXER
C***REVISION HISTORY  (YYMMDD)
C   870810  DATE WRITTEN
C   910619  Modified to meet SLATEC code and prologue standards.  (BKS)
C***END PROLOGUE  CCHK22
C     .. Parameters ..
      COMPLEX            ZERO, HALF
      PARAMETER          ( ZERO = ( 0.0, 0.0 ), HALF = ( 0.5, 0.0 ) )
      REAL               RZERO
      PARAMETER          ( RZERO = 0.0 )
C     .. Scalar Arguments ..
      LOGICAL            FATAL
      REAL               EPS, THRESH
      INTEGER            INCMAX, KPRINT, NALF, NBET, NIDIM, NINC, NKB,
     $                   NMAX, NOUT
      CHARACTER*6        SNAME
C     .. Array Arguments ..
      COMPLEX            A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ),
     $                   AS( NMAX*NMAX ), BET( NBET ), X( NMAX ),
     $                   XS( NMAX*INCMAX ), XX( NMAX*INCMAX ),
     $                   Y( NMAX ), YS( NMAX*INCMAX ), YT( NMAX ),
     $                   YY( NMAX*INCMAX )
      REAL               G( NMAX )
      INTEGER            IDIM( NIDIM ), INC( NINC ), KB( NKB )
C     .. Local Scalars ..
      COMPLEX            ALPHA, ALS, BETA, BLS, TRANSL
      REAL               ERR, ERRMAX
      INTEGER            I, IA, IB, IC, IK, IN, INCX, INCXS, INCY,
     $                   INCYS, IX, IY, K, KS, LAA, LDA, LDAS, LX, LY,
     $                   N, NARGS, NC, NERR, NK, NS
      LOGICAL            BANDED, FTL, FULL, NULL, PACKED, RESET
      CHARACTER*1        UPLO, UPLOS
      CHARACTER*2        ICH
C     .. Local Arrays ..
      LOGICAL            ISAME( 13 )
C     .. External Functions ..
      INTEGER            NUMXER
      LOGICAL            LCE, LCERES
      EXTERNAL           LCE, LCERES, NUMXER
C     .. External Subroutines ..
      EXTERNAL           CHBMV, CHEMV, CHPMV, CMAKE2, CMVCH
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
C     .. Data statements ..
      DATA               ICH/'UL'/
C***FIRST EXECUTABLE STATEMENT  CCHK22
      FULL = SNAME( 3: 3 ).EQ.'E'
      BANDED = SNAME( 3: 3 ).EQ.'B'
      PACKED = SNAME( 3: 3 ).EQ.'P'
C     Define the number of arguments.
      IF( FULL )THEN
         NARGS = 10
      ELSE IF( BANDED )THEN
         NARGS = 11
      ELSE IF( PACKED )THEN
         NARGS = 9
      END IF
C
      NC = 0
      RESET = .TRUE.
      ERRMAX = RZERO
C
      DO 110 IN = 1, NIDIM
         N = IDIM( IN )
C
         IF( BANDED )THEN
            NK = NKB
         ELSE
            NK = 1
         END IF
         DO 100 IK = 1, NK
            IF( BANDED )THEN
               K = KB( IK )
            ELSE
               K = N - 1
            END IF
C           Set LDA to 1 more than minimum value if room.
            IF( BANDED )THEN
               LDA = K + 1
            ELSE
               LDA = N
            END IF
            IF( LDA.LT.NMAX )
     $         LDA = LDA + 1
C           Skip tests if not enough room.
            IF( LDA.GT.NMAX )
     $         GO TO 100
            IF( PACKED )THEN
               LAA = ( N*( N + 1 ) )/2
            ELSE
               LAA = LDA*N
            END IF
            NULL = N.LE.0
C
            DO 90 IC = 1, 2
               UPLO = ICH( IC: IC )
C
C              Generate the matrix A.
C
               TRANSL = ZERO
               CALL CMAKE2( SNAME( 2: 3 ), UPLO, ' ', N, N, A, NMAX, AA,
     $                     LDA, K, K, RESET, TRANSL )
C
               DO 80 IX = 1, NINC
                  INCX = INC( IX )
                  LX = ABS( INCX )*N
C
C                 Generate the vector X.
C
                  TRANSL = HALF
                  CALL CMAKE2( 'GE', ' ', ' ', 1, N, X, 1, XX,
     $                        ABS( INCX ), 0, N - 1, RESET, TRANSL )
                  IF( N.GT.1 )THEN
                     X( N/2 ) = ZERO
                     XX( 1 + ABS( INCX )*( N/2 - 1 ) ) = ZERO
                  END IF
C
                  DO 70 IY = 1, NINC
                     INCY = INC( IY )
                     LY = ABS( INCY )*N
C
                     DO 60 IA = 1, NALF
                        ALPHA = ALF( IA )
C
                        DO 50 IB = 1, NBET
                           BETA = BET( IB )
C
C                          Generate the vector Y.
C
                           TRANSL = ZERO
                           CALL CMAKE2( 'GE', ' ', ' ', 1, N, Y, 1, YY,
     $                                 ABS( INCY ), 0, N - 1, RESET,
     $                                 TRANSL )
C
                           NC = NC + 1
C
C                          Save every datum before calling the
C                          subroutine.
C
                           UPLOS = UPLO
                           NS = N
                           KS = K
                           ALS = ALPHA
                           DO 10 I = 1, LAA
                              AS( I ) = AA( I )
   10                      CONTINUE
                           LDAS = LDA
                           DO 20 I = 1, LX
                              XS( I ) = XX( I )
   20                      CONTINUE
                           INCXS = INCX
                           BLS = BETA
                           DO 30 I = 1, LY
                              YS( I ) = YY( I )
   30                      CONTINUE
                           INCYS = INCY
C
C                          Call the subroutine.
C
                           IF( FULL )THEN
                              CALL CHEMV( UPLO, N, ALPHA, AA, LDA, XX,
     $                                    INCX, BETA, YY, INCY )
                           ELSE IF( BANDED )THEN
                              CALL CHBMV( UPLO, N, K, ALPHA, AA, LDA,
     $                                    XX, INCX, BETA, YY, INCY )
                           ELSE IF( PACKED )THEN
                              CALL CHPMV( UPLO, N, ALPHA, AA, XX, INCX,
     $                                    BETA, YY, INCY )
                           END IF
C
C                          Check if error-exit was taken incorrectly.
C
                           IF( NUMXER(NERR) .NE. 0 )THEN
                             IF (KPRINT .GE. 2) THEN
                              WRITE( NOUT, FMT = 9992 )
                             ENDIF
                             FATAL = .TRUE.
                           END IF
C
C                          See what data changed inside subroutines.
C
                           ISAME( 1 ) = UPLO.EQ.UPLOS
                           ISAME( 2 ) = NS.EQ.N
                           IF( FULL )THEN
                              ISAME( 3 ) = ALS.EQ.ALPHA
                              ISAME( 4 ) = LCE( AS, AA, LAA )
                              ISAME( 5 ) = LDAS.EQ.LDA
                              ISAME( 6 ) = LCE( XS, XX, LX )
                              ISAME( 7 ) = INCXS.EQ.INCX
                              ISAME( 8 ) = BLS.EQ.BETA
                              IF( NULL )THEN
                                 ISAME( 9 ) = LCE( YS, YY, LY )
                              ELSE
                                 ISAME( 9 ) = LCERES( 'GE', ' ', 1, N,
     $                                        YS, YY, ABS( INCY ) )
                              END IF
                              ISAME( 10 ) = INCYS.EQ.INCY
                           ELSE IF( BANDED )THEN
                              ISAME( 3 ) = KS.EQ.K
                              ISAME( 4 ) = ALS.EQ.ALPHA
                              ISAME( 5 ) = LCE( AS, AA, LAA )
                              ISAME( 6 ) = LDAS.EQ.LDA
                              ISAME( 7 ) = LCE( XS, XX, LX )
                              ISAME( 8 ) = INCXS.EQ.INCX
                              ISAME( 9 ) = BLS.EQ.BETA
                              IF( NULL )THEN
                                 ISAME( 10 ) = LCE( YS, YY, LY )
                              ELSE
                                 ISAME( 10 ) = LCERES( 'GE', ' ', 1, N,
     $                                         YS, YY, ABS( INCY ) )
                              END IF
                              ISAME( 11 ) = INCYS.EQ.INCY
                           ELSE IF( PACKED )THEN
                              ISAME( 3 ) = ALS.EQ.ALPHA
                              ISAME( 4 ) = LCE( AS, AA, LAA )
                              ISAME( 5 ) = LCE( XS, XX, LX )
                              ISAME( 6 ) = INCXS.EQ.INCX
                              ISAME( 7 ) = BLS.EQ.BETA
                              IF( NULL )THEN
                                 ISAME( 8 ) = LCE( YS, YY, LY )
                              ELSE
                                 ISAME( 8 ) = LCERES( 'GE', ' ', 1, N,
     $                                        YS, YY, ABS( INCY ) )
                              END IF
                              ISAME( 9 ) = INCYS.EQ.INCY
                           END IF
C
C                          If data was incorrectly changed, report and
C                          return.
C
                           DO 40 I = 1, NARGS
                             IF (.NOT. ISAME( I )) THEN
                               FATAL = .TRUE.
                               IF (KPRINT .GE. 2) THEN
                                 WRITE( NOUT, FMT = 9998 )I
                               ENDIF
                             ENDIF
  40                       CONTINUE
C
                           FTL = .FALSE.
                           IF( .NOT.NULL )THEN
C
C                             Check the result.
C
                              CALL CMVCH( 'N', N, N, ALPHA, A, NMAX, X,
     $                                    INCX, BETA, Y, INCY, YT, G,
     $                                    YY, EPS, ERR, FTL, NOUT,
     $                                    .TRUE., KPRINT )
                              ERRMAX = MAX( ERRMAX, ERR )
                           END IF
                           IF (FTL) THEN
                            FATAL = .TRUE.
                            IF (KPRINT .GE. 3) THEN
                              WRITE (NOUT, FMT = 9996) SNAME
                              IF( FULL )THEN
                                 WRITE( NOUT, FMT = 9994 )NC, SNAME,
     $                              UPLO, N, ALPHA, LDA,
     $                              INCX, BETA, INCY
                              ELSE IF( BANDED )THEN
                                 WRITE( NOUT, FMT = 9995 )NC, SNAME,
     $                              UPLO, N,
     $                              ALPHA, LDA, INCX, BETA, INCY
                              ELSE IF( PACKED )THEN
                                 WRITE( NOUT, FMT = 9995 )NC, SNAME,
     $                              UPLO, N, ALPHA, INCX,
     $                              BETA, INCY
                              END IF
                            ENDIF
                          ENDIF
C
   50                   CONTINUE
C
   60                CONTINUE
C
   70             CONTINUE
C
   80          CONTINUE
C
   90       CONTINUE
C
  100    CONTINUE
C
  110 CONTINUE
C
C     Report result.
C
      IF (.NOT. FATAL) THEN
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
 9995 FORMAT( 1X, I6, ': ', A6, '(''', A1, ''',', I3, ',(', F4.1, ',',
     $      F4.1, '), AP, X,', I2, ',(', F4.1, ',', F4.1, '), Y,', I2,
     $      ')                .' )
 9994 FORMAT( 1X, I6, ': ', A6, '(''', A1, ''',', 2( I3, ',' ), '(',
     $      F4.1, ',', F4.1, '), A,', I3, ', X,', I2, ',(', F4.1, ',',
     $      F4.1, '), Y,', I2, ')         .' )
 9992 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',
     $      '******' )
C
C     End of CCHK22.
C
      END
