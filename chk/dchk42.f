*DECK DCHK42
      SUBROUTINE DCHK42 (SNAME, EPS, THRESH, NOUT, KPRINT, FATAL, NIDIM,
     $   IDIM, NALF, ALF, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS,
     $   Y, YY, YS, YT, G, Z)
C***BEGIN PROLOGUE  DCHK42
C***SUBSIDIARY
C***PURPOSE  Test DGER.
C***LIBRARY   SLATEC (BLAS)
C***KEYWORDS  BLAS, QUICK CHECK SERVICE ROUTINE
C***AUTHOR  Du Croz, J. (NAG)
C           Hanson, R. J. (SNLA)
C***DESCRIPTION
C
C  Quick check for DGER.
C
C  Auxiliary routine for test program for Level 2 Blas.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DGER, DMAKE2, DMVCH, LDE, LDERES, NUMXER
C***REVISION HISTORY  (YYMMDD)
C   870810  DATE WRITTEN
C   910619  Modified to meet SLATEC code and prologue standards. (BKS)
C***END PROLOGUE  DCHK42
C     .. Parameters ..
      DOUBLE PRECISION   ZERO, HALF, ONE
      PARAMETER          ( ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0 )
C     .. Scalar Arguments ..
      LOGICAL            FATAL
      DOUBLE PRECISION   EPS, THRESH
      INTEGER            INCMAX, KPRINT, NALF, NIDIM, NINC, NMAX, NOUT
      CHARACTER*6        SNAME
C     .. Array Arguments ..
      DOUBLE PRECISION   A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ),
     $                   AS( NMAX*NMAX ), G( NMAX ), X( NMAX ),
     $                   XS( NMAX*INCMAX ), XX( NMAX*INCMAX ),
     $                   Y( NMAX ), YS( NMAX*INCMAX ), YT( NMAX ),
     $                   YY( NMAX*INCMAX ), Z( NMAX )
      INTEGER            IDIM( NIDIM ), INC( NINC )
C     .. Local Scalars ..
      DOUBLE PRECISION   ALPHA, ALS, ERR, ERRMAX, TRANSL
      INTEGER            I, IA, IM, IN, INCX, INCXS, INCY, INCYS, IX,
     $                   IY, J, LAA, LDA, LDAS, LX, LY, M, MS, N, NARGS,
     $                   NC, ND, NS, NERR
      LOGICAL            FTL, NULL, RESET
C     .. Local Arrays ..
      DOUBLE PRECISION   W( 1 )
      LOGICAL            ISAME( 13 )
C     .. External Functions ..
      INTEGER            NUMXER
      LOGICAL            LDE, LDERES
      EXTERNAL           LDE, LDERES, NUMXER
C     .. External Subroutines ..
      EXTERNAL           DGER, DMAKE2, DMVCH
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
C***FIRST EXECUTABLE STATEMENT  DCHK42
C     Define the number of arguments.
      NARGS = 9
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
C           Set LDA to 1 more than minimum value if room.
            LDA = M
            IF( LDA.LT.NMAX )
     $         LDA = LDA + 1
C           Skip tests if not enough room.
            IF( LDA.GT.NMAX )
     $         GO TO 110
            LAA = LDA*N
            NULL = N.LE.0.OR.M.LE.0
C
            DO 100 IX = 1, NINC
               INCX = INC( IX )
               LX = ABS( INCX )*M
C
C              Generate the vector X.
C
               TRANSL = HALF
               CALL DMAKE2( 'GE', ' ', ' ', 1, M, X, 1, XX, ABS( INCX ),
     $                     0, M - 1, RESET, TRANSL )
               IF( M.GT.1 )THEN
                  X( M/2 ) = ZERO
                  XX( 1 + ABS( INCX )*( M/2 - 1 ) ) = ZERO
               END IF
C
               DO 90 IY = 1, NINC
                  INCY = INC( IY )
                  LY = ABS( INCY )*N
C
C                 Generate the vector Y.
C
                  TRANSL = ZERO
                  CALL DMAKE2( 'GE', ' ', ' ', 1, N, Y, 1, YY,
     $                        ABS( INCY ), 0, N - 1, RESET, TRANSL )
                  IF( N.GT.1 )THEN
                     Y( N/2 ) = ZERO
                     YY( 1 + ABS( INCY )*( N/2 - 1 ) ) = ZERO
                  END IF
C
                  DO 80 IA = 1, NALF
                     ALPHA = ALF( IA )
C
C                    Generate the matrix A.
C
                     TRANSL = ZERO
                     CALL DMAKE2(SNAME( 2: 3 ), ' ', ' ', M, N, A, NMAX,
     $                           AA, LDA, M - 1, N - 1, RESET, TRANSL )
C
                     NC = NC + 1
C
C                    Save every datum before calling the subroutine.
C
                     MS = M
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
                     CALL DGER( M, N, ALPHA, XX, INCX, YY, INCY, AA,
     $                          LDA )
C
C                    Check if error-exit was taken incorrectly.
C
                     IF(  NUMXER(NERR) .NE. 0 )THEN
                      IF (KPRINT .GE. 2) THEN
                        WRITE( NOUT, FMT = 9993 )
                      ENDIF
                        FATAL = .TRUE.
                     END IF
C
C                    See what data changed inside subroutine.
C
                     ISAME( 1 ) = MS.EQ.M
                     ISAME( 2 ) = NS.EQ.N
                     ISAME( 3 ) = ALS.EQ.ALPHA
                     ISAME( 4 ) = LDE( XS, XX, LX )
                     ISAME( 5 ) = INCXS.EQ.INCX
                     ISAME( 6 ) = LDE( YS, YY, LY )
                     ISAME( 7 ) = INCYS.EQ.INCY
                     IF( NULL )THEN
                        ISAME( 8 ) = LDE( AS, AA, LAA )
                     ELSE
                        ISAME( 8 ) = LDERES( 'GE', ' ', M, N, AS, AA,
     $                               LDA )
                     END IF
                     ISAME( 9 ) = LDAS.EQ.LDA
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
                           DO 50 I = 1, M
                              Z( I ) = X( I )
   50                      CONTINUE
                        ELSE
                           DO 60 I = 1, M
                              Z( I ) = X( M - I + 1 )
   60                      CONTINUE
                        END IF
                        DO 70 J = 1, N
                           IF( INCY.GT.0 )THEN
                              W( 1 ) = Y( J )
                           ELSE
                              W( 1 ) = Y( N - J + 1 )
                           END IF
                           CALL DMVCH( 'N', M, 1, ALPHA, Z, NMAX, W, 1,
     $                                 ONE, A( 1, J ), 1, YT, G,
     $                                 AA( 1 + ( J - 1 )*LDA ), EPS,
     $                                 ERR, FTL, NOUT, .TRUE., KPRINT)
                           ERRMAX = MAX( ERRMAX, ERR )
   70                   CONTINUE
                     END IF
                     IF (FTL) THEN
                       FATAL = .TRUE.
                       IF (KPRINT .GE. 3) THEN
                         WRITE( NOUT, FMT = 9995 )J
                         WRITE( NOUT, FMT = 9996 )SNAME
                         WRITE( NOUT, FMT = 9994 )NC,
     $                     SNAME, M, N, ALPHA, INCX, INCY,
     $                     LDA
                       ENDIF
                     ENDIF
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
 9995 FORMAT( '      THESE ARE THE RESULTS FOR COLUMN ', I3 )
 9994 FORMAT( 1X, I6, ': ', A6, '(', 2( I3, ',' ), F4.1, ', X,', I2,
     $      ', Y,', I2, ', A,', I3, ')                  .' )
 9993 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',
     $      '******' )
C
C     End of DCHK42.
C
      END
