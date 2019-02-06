*DECK SCHK32
      SUBROUTINE SCHK32 (SNAME, EPS, THRESH, NOUT, KPRINT, FATAL, NIDIM,
     $   IDIM, NKB, KB, NINC, INC, NMAX, INCMAX, A, AA, AS, X, XX, XS,
     $   XT, G, Z)
C***BEGIN PROLOGUE  SCHK32
C***SUBSIDIARY
C***PURPOSE  Quick check for STRMV, STBMV, STPMV, STRSV, STBSV and
C            STPSV.
C***LIBRARY   SLATEC (BLAS)
C***KEYWORDS  BLAS, QUICK CHECK SERVICE ROUTINE
C***AUTHOR  Du Croz, J. (NAG)
C           Hanson, R. J. (SNLA)
C***DESCRIPTION
C
C  Quick check for STRMV, STBMV, STPMV, STRSV, STBSV and STPSV.
C
C  Auxiliary routine for test program for Level 2 Blas.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  LSE, LSERES, NUMXER, SMAKE2, SMVCH, STBMV, STBSV,
C                    STPMV, STPSV, STRMV, STRSV
C***REVISION HISTORY  (YYMMDD)
C   870810  DATE WRITTEN
C   910619  Modified to meet SLATEC code and prologue standards.  (BKS)
C***END PROLOGUE  SCHK32
C     .. Parameters ..
      REAL               ZERO, HALF, ONE
      PARAMETER          ( ZERO = 0.0, HALF = 0.5, ONE = 1.0 )
C     .. Scalar Arguments ..
      LOGICAL            FATAL
      REAL               EPS, THRESH
      INTEGER            INCMAX, KPRINT, NIDIM, NINC, NKB, NMAX, NOUT
      CHARACTER*6        SNAME
C     .. Array Arguments ..
      REAL               A( NMAX, NMAX ), AA( NMAX*NMAX ),
     $                   AS( NMAX*NMAX ), G( NMAX ), X( NMAX ),
     $                   XS( NMAX*INCMAX ), XT( NMAX ),
     $                   XX( NMAX*INCMAX ), Z( NMAX )
      INTEGER            IDIM( NIDIM ), INC( NINC ), KB( NKB )
C     .. Local Scalars ..
      REAL               ERR, ERRMAX, TRANSL
      INTEGER            I, ICD, ICT, ICU, IK, IN, INCX, INCXS, IX, K,
     $                   KS, LAA, LDA, LDAS, LX, N, NARGS, NC, NK, NS,
     $                   NERR
      LOGICAL            BANDED, FTL, FULL, NULL, PACKED, RESET
      CHARACTER*1        DIAG, DIAGS, TRANS, TRANSS, UPLO, UPLOS
      CHARACTER*2        ICHD, ICHU
      CHARACTER*3        ICHT
C     .. Local Arrays ..
      LOGICAL            ISAME( 13 )
C     .. External Functions ..
      INTEGER            NUMXER
      LOGICAL            LSE, LSERES
      EXTERNAL           LSE, LSERES, NUMXER
C     .. External Subroutines ..
      EXTERNAL           SMAKE2, SMVCH, STBMV, STBSV, STPMV, STPSV,
     $                   STRMV, STRSV
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
C     .. Data statements ..
      DATA               ICHU/'UL'/, ICHT/'NTC'/, ICHD/'UN'/
C***FIRST EXECUTABLE STATEMENT  SCHK32
      FULL = SNAME( 3: 3 ).EQ.'R'
      BANDED = SNAME( 3: 3 ).EQ.'B'
      PACKED = SNAME( 3: 3 ).EQ.'P'
C     Define the number of arguments.
      IF( FULL )THEN
         NARGS = 8
      ELSE IF( BANDED )THEN
         NARGS = 9
      ELSE IF( PACKED )THEN
         NARGS = 7
      END IF
C
      NC = 0
      RESET = .TRUE.
      ERRMAX = ZERO
C     Set up zero vector for SMVCH.
      DO 10 I = 1, NMAX
         Z( I ) = ZERO
   10 CONTINUE
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
            DO 90 ICU = 1, 2
               UPLO = ICHU( ICU: ICU )
C
               DO 80 ICT = 1, 3
                  TRANS = ICHT( ICT: ICT )
C
                  DO 70 ICD = 1, 2
                     DIAG = ICHD( ICD: ICD )
C
C                    Generate the matrix A.
C
                     TRANSL = ZERO
                     CALL SMAKE2( SNAME( 2: 3 ), UPLO, DIAG, N, N, A,
     $                           NMAX, AA, LDA, K, K, RESET, TRANSL )
C
                     DO 60 IX = 1, NINC
                        INCX = INC( IX )
                        LX = ABS( INCX )*N
C
C                       Generate the vector X.
C
                        TRANSL = HALF
                        CALL SMAKE2( 'GE', ' ', ' ', 1, N, X, 1, XX,
     $                              ABS( INCX ), 0, N - 1, RESET,
     $                              TRANSL )
                        IF( N.GT.1 )THEN
                           X( N/2 ) = ZERO
                           XX( 1 + ABS( INCX )*( N/2 - 1 ) ) = ZERO
                        END IF
C
                        NC = NC + 1
C
C                       Save every datum before calling the subroutine.
C
                        UPLOS = UPLO
                        TRANSS = TRANS
                        DIAGS = DIAG
                        NS = N
                        KS = K
                        DO 20 I = 1, LAA
                           AS( I ) = AA( I )
   20                   CONTINUE
                        LDAS = LDA
                        DO 30 I = 1, LX
                           XS( I ) = XX( I )
   30                   CONTINUE
                        INCXS = INCX
C
C                       Call the subroutine.
C
                        IF( SNAME( 4: 5 ).EQ.'MV' )THEN
                           IF( FULL )THEN
                              CALL STRMV( UPLO, TRANS, DIAG, N, AA, LDA,
     $                                    XX, INCX )
                           ELSE IF( BANDED )THEN
                              CALL STBMV( UPLO, TRANS, DIAG, N, K, AA,
     $                                    LDA, XX, INCX )
                           ELSE IF( PACKED )THEN
                              CALL STPMV( UPLO, TRANS, DIAG, N, AA, XX,
     $                                    INCX )
                           END IF
                        ELSE IF( SNAME( 4: 5 ).EQ.'SV' )THEN
                           IF( FULL )THEN
                              CALL STRSV( UPLO, TRANS, DIAG, N, AA, LDA,
     $                                    XX, INCX )
                           ELSE IF( BANDED )THEN
                              CALL STBSV( UPLO, TRANS, DIAG, N, K, AA,
     $                                    LDA, XX, INCX )
                           ELSE IF( PACKED )THEN
                              CALL STPSV( UPLO, TRANS, DIAG, N, AA, XX,
     $                                    INCX )
                           END IF
                        END IF
C
C                       Check if error-exit was taken incorrectly.
C
                        IF (NUMXER(NERR) .NE. 0) THEN
                           IF (KPRINT .GE. 2) THEN
                             WRITE( NOUT, FMT = 9992 )
                           ENDIF
                           FATAL = .TRUE.
                        END IF
C
C                       See what data changed inside subroutines.
C
                        ISAME( 1 ) = UPLO.EQ.UPLOS
                        ISAME( 2 ) = TRANS.EQ.TRANSS
                        ISAME( 3 ) = DIAG.EQ.DIAGS
                        ISAME( 4 ) = NS.EQ.N
                        IF( FULL )THEN
                           ISAME( 5 ) = LSE( AS, AA, LAA )
                           ISAME( 6 ) = LDAS.EQ.LDA
                           IF( NULL )THEN
                              ISAME( 7 ) = LSE( XS, XX, LX )
                           ELSE
                              ISAME( 7 ) = LSERES( 'GE', ' ', 1, N, XS,
     $                                     XX, ABS( INCX ) )
                           END IF
                           ISAME( 8 ) = INCXS.EQ.INCX
                        ELSE IF( BANDED )THEN
                           ISAME( 5 ) = KS.EQ.K
                           ISAME( 6 ) = LSE( AS, AA, LAA )
                           ISAME( 7 ) = LDAS.EQ.LDA
                           IF( NULL )THEN
                              ISAME( 8 ) = LSE( XS, XX, LX )
                           ELSE
                              ISAME( 8 ) = LSERES( 'GE', ' ', 1, N, XS,
     $                                     XX, ABS( INCX ) )
                           END IF
                           ISAME( 9 ) = INCXS.EQ.INCX
                        ELSE IF( PACKED )THEN
                           ISAME( 5 ) = LSE( AS, AA, LAA )
                           IF( NULL )THEN
                              ISAME( 6 ) = LSE( XS, XX, LX )
                           ELSE
                              ISAME( 6 ) = LSERES( 'GE', ' ', 1, N, XS,
     $                                     XX, ABS( INCX ) )
                           END IF
                           ISAME( 7 ) = INCXS.EQ.INCX
                        END IF
C
C                       If data was incorrectly changed, report and
C                       return.
C
                        DO 40 I = 1, NARGS
                          IF (.NOT. ISAME( I )) THEN
                            FATAL = .TRUE.
                            IF (KPRINT .GE. 2) THEN
                              WRITE( NOUT, FMT = 9998 )I
                            ENDIF
                          ENDIF
  40                    CONTINUE
C
                        FTL = .FALSE.
                        IF( .NOT.NULL )THEN
                           IF( SNAME( 4: 5 ).EQ.'MV' )THEN
C
C                             Check the result.
C
                              CALL SMVCH( TRANS, N, N, ONE, A, NMAX, X,
     $                                    INCX, ZERO, Z, INCX, XT, G,
     $                                    XX, EPS, ERR, FTL, NOUT,
     $                                    .TRUE., KPRINT )
                           ELSE IF( SNAME( 4: 5 ).EQ.'SV' )THEN
C
C                             Compute approximation to original vector.
C
                              DO 50 I = 1, N
                                 Z( I ) = XX( 1 + ( I - 1 )*
     $                                    ABS( INCX ) )
                                 XX( 1 + ( I - 1 )*ABS( INCX ) )
     $                              = X( I )
   50                         CONTINUE
                              CALL SMVCH( TRANS, N, N, ONE, A, NMAX, Z,
     $                                    INCX, ZERO, X, INCX, XT, G,
     $                                    XX, EPS, ERR, FTL, NOUT,
     $                                    .FALSE., KPRINT )
                           END IF
                           ERRMAX = MAX( ERRMAX, ERR )
                        END IF
                        IF (FTL) THEN
                          FATAL = .TRUE.
                          IF (KPRINT .GE. 3) THEN
                            WRITE( NOUT, FMT = 9996 )SNAME
                            IF( FULL )THEN
                               WRITE( NOUT, FMT = 9993 )NC, SNAME,
     $                          UPLO, TRANS, DIAG, N,
     $                          LDA, INCX
                            ELSE IF( BANDED )THEN
                               WRITE( NOUT, FMT = 9994 )NC, SNAME,
     $                          UPLO, TRANS, DIAG, N,
     $                          K, LDA, INCX
                            ELSE IF( PACKED )THEN
                               WRITE( NOUT, FMT = 9995 )NC, SNAME,
     $                          UPLO, TRANS, DIAG, N,
     $                          INCX
                            END IF
                          ENDIF
                        ENDIF
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
 9995 FORMAT( 1X, I6, ': ', A6, '(', 3( '''', A1, ''',' ), I3, ', AP, ',
     $      'X,', I2, ')                        .' )
 9994 FORMAT( 1X, I6, ': ', A6, '(', 3( '''', A1, ''',' ), 2( I3, ',' ),
     $      ' A,', I3, ', X,', I2, ')                 .' )
 9993 FORMAT( 1X, I6, ': ', A6, '(', 3( '''', A1, ''',' ), I3, ', A,',
     $      I3, ', X,', I2, ')                     .' )
 9992 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',
     $      '******' )
C
C     End of SCHK32.
C
      END
