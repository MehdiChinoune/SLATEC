*DECK CCHK33
      SUBROUTINE CCHK33 (SNAME, EPS, THRESH, NOUT, KPRINT, FATAL, NIDIM,
     $   IDIM, NALF, ALF, NMAX, A, AA, AS, B, BB, BS, CT, G, C)
C***BEGIN PROLOGUE  CCHK33
C***SUBSIDIARY
C***PURPOSE  Quick check for CTRMM and CTRSM.
C***LIBRARY   SLATEC (BLAS)
C***KEYWORDS  BLAS, QUICK CHECK SERVICE ROUTINE
C***AUTHOR  Dongarra, J. J., (ANL)
C           Duff, I., (AERE)
C           Du Croz, J., (NAG)
C           Hammarling, S., (NAG)
C***DESCRIPTION
C
C  Quick check for CTRMM and CTRSM.
C
C  Auxiliary routine for test program for Level 3 Blas.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CMAKE3, CMMCH, CTRMM, CTRSM, LCE, LCERES, NUMXER
C***REVISION HISTORY  (YYMMDD)
C   890208  DATE WRITTEN
C   910619  Modified to meet SLATEC code and prologue standards.  (BKS)
C***END PROLOGUE  CCHK33
C     .. Parameters ..
      COMPLEX            ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0, 0.0 ), ONE = ( 1.0, 0.0 ) )
      REAL               RZERO
      PARAMETER          ( RZERO = 0.0 )
C     .. Scalar Arguments ..
      LOGICAL            FATAL
      REAL               EPS, THRESH
      INTEGER            KPRINT, NALF, NIDIM, NMAX, NOUT
      CHARACTER*6        SNAME
C     .. Array Arguments ..
      COMPLEX            A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ),
     $                   AS( NMAX*NMAX ), B( NMAX, NMAX ),
     $                   BB( NMAX*NMAX ), BS( NMAX*NMAX ),
     $                   C( NMAX, NMAX ),
     $                   CT( NMAX )
      REAL               G( NMAX )
      INTEGER            IDIM( NIDIM )
C     .. Local Scalars ..
      COMPLEX            ALPHA, ALS
      REAL               ERR, ERRMAX
      INTEGER            I, IA, ICD, ICS, ICT, ICU, IM, IN, J, LAA, LBB,
     $                   LDA, LDAS, LDB, LDBS, M, MS, N, NA, NARGS, NC,
     $                   NERR, NS
      LOGICAL            FTL, LEFT, NULL, RESET
      CHARACTER*1        DIAG, DIAGS, SIDE, SIDES, TRANAS, TRANSA, UPLO,
     $                   UPLOS
      CHARACTER*2        ICHU, ICHS, ICHD
      CHARACTER*3        ICHT
C     .. Local Arrays ..
      LOGICAL            ISAME( 13 )
C     .. External Functions ..
      INTEGER            NUMXER
      LOGICAL            LCE, LCERES
      EXTERNAL           LCE, LCERES, NUMXER
C     .. External Subroutines ..
      EXTERNAL           CTRMM, CTRSM, CMAKE3, CMMCH
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
C     .. Data statements ..
      DATA               ICHS/'LR'/, ICHU/'UL'/, ICHT/'NTC'/, ICHD/'UN'/
C***FIRST EXECUTABLE STATEMENT  CCHK33
      NARGS = 11
      NC = 0
      RESET = .TRUE.
      ERRMAX = RZERO
C     Set up zero matrix for CMMCH.
      DO 20 J = 1, NMAX
         DO 10 I = 1, NMAX
            C( I, J ) = ZERO
   10    CONTINUE
   20 CONTINUE
C
      DO 140 IM = 1, NIDIM
         M = IDIM( IM )
C
         DO 130 IN = 1, NIDIM
            N = IDIM( IN )
C           Set LDB to 1 more than minimum value if room.
            LDB = M
            IF( LDB.LT.NMAX )
     $         LDB = LDB + 1
C           Skip tests if not enough room.
            IF( LDB.GT.NMAX )
     $         GO TO 130
            LBB = LDB*N
            NULL = M.LE.0.OR.N.LE.0
C
            DO 120 ICS = 1, 2
               SIDE = ICHS( ICS: ICS )
               LEFT = SIDE.EQ.'L'
               IF( LEFT )THEN
                  NA = M
               ELSE
                  NA = N
               END IF
C              Set LDA to 1 more than minimum value if room.
               LDA = NA
               IF( LDA.LT.NMAX )
     $            LDA = LDA + 1
C              Skip tests if not enough room.
               IF( LDA.GT.NMAX )
     $            GO TO 130
               LAA = LDA*NA
C
               DO 110 ICU = 1, 2
                  UPLO = ICHU( ICU: ICU )
C
                  DO 100 ICT = 1, 3
                     TRANSA = ICHT( ICT: ICT )
C
                     DO 90 ICD = 1, 2
                        DIAG = ICHD( ICD: ICD )
C
                        DO 80 IA = 1, NALF
                           ALPHA = ALF( IA )
C
C                          Generate the matrix A.
C
                           CALL CMAKE3( 'TR', UPLO, DIAG, NA, NA, A,
     $                                 NMAX, AA, LDA, RESET, ZERO )
C
C                          Generate the matrix B.
C
                           CALL CMAKE3( 'GE', ' ', ' ', M, N, B, NMAX,
     $                                 BB, LDB, RESET, ZERO )
C
                           NC = NC + 1
C
C                          Save every datum before calling the
C                          subroutine.
C
                           SIDES = SIDE
                           UPLOS = UPLO
                           TRANAS = TRANSA
                           DIAGS = DIAG
                           MS = M
                           NS = N
                           ALS = ALPHA
                           DO 30 I = 1, LAA
                              AS( I ) = AA( I )
   30                      CONTINUE
                           LDAS = LDA
                           DO 40 I = 1, LBB
                              BS( I ) = BB( I )
   40                      CONTINUE
                           LDBS = LDB
C
C                          Call the subroutine.
C
                           IF( SNAME( 4: 5 ).EQ.'MM' )THEN
                              CALL CTRMM( SIDE, UPLO, TRANSA, DIAG, M,
     $                                    N, ALPHA, AA, LDA, BB, LDB )
                           ELSE IF( SNAME( 4: 5 ).EQ.'SM' )THEN
                              CALL CTRSM( SIDE, UPLO, TRANSA, DIAG, M,
     $                                    N, ALPHA, AA, LDA, BB, LDB )
                           END IF
C
C                          Check if error-exit was taken incorrectly.
C
                           IF( NUMXER(NERR) .NE. 0 )THEN
                             IF (KPRINT .GE. 2) THEN
                              WRITE( NOUT, FMT = 9994 )
                              END IF
                             FATAL = .TRUE.
                           END IF
C
C                          See what data changed inside subroutines.
C
                           ISAME( 1 ) = SIDES.EQ.SIDE
                           ISAME( 2 ) = UPLOS.EQ.UPLO
                           ISAME( 3 ) = TRANAS.EQ.TRANSA
                           ISAME( 4 ) = DIAGS.EQ.DIAG
                           ISAME( 5 ) = MS.EQ.M
                           ISAME( 6 ) = NS.EQ.N
                           ISAME( 7 ) = ALS.EQ.ALPHA
                           ISAME( 8 ) = LCE( AS, AA, LAA )
                           ISAME( 9 ) = LDAS.EQ.LDA
                           IF( NULL )THEN
                              ISAME( 10 ) = LCE( BS, BB, LBB )
                           ELSE
                              ISAME( 10 ) = LCERES( 'GE', ' ', M, N, BS,
     $                                      BB, LDB )
                           END IF
                           ISAME( 11 ) = LDBS.EQ.LDB
C
C                          If data was incorrectly changed, report and
C                          return.
C
                           DO 50 I = 1, NARGS
                             IF (.NOT. ISAME( I )) THEN
                               FATAL = .TRUE.
                               IF (KPRINT .GE. 2) THEN
                                 WRITE( NOUT, FMT = 9998 )I
                               ENDIF
                             ENDIF
  50                          CONTINUE
C
                           FTL = .FALSE.
                           IF( .NOT.NULL )THEN
                              IF( SNAME( 4: 5 ).EQ.'MM' )THEN
C
C                                Check the result.
C
                                 IF( LEFT )THEN
                                    CALL CMMCH( TRANSA, 'N', M, N, M,
     $                                          ALPHA, A, NMAX, B, NMAX,
     $                                          ZERO, C, NMAX, CT, G,
     $                                          BB, LDB, EPS, ERR,
     $                                          FTL, NOUT, .TRUE.,
     $                                          KPRINT )
                                 ELSE
                                    CALL CMMCH( 'N', TRANSA, M, N, N,
     $                                          ALPHA, B, NMAX, A, NMAX,
     $                                          ZERO, C, NMAX, CT, G,
     $                                          BB, LDB, EPS, ERR,
     $                                          FTL, NOUT, .TRUE.,
     $                                          KPRINT )
                                 END IF
                              ELSE IF( SNAME( 4: 5 ).EQ.'SM' )THEN
C
C                                Compute approximation to original
C                                matrix.
C
                                 DO 70 J = 1, N
                                    DO 60 I = 1, M
                                       C( I, J ) = BB( I + ( J - 1 )*
     $                                             LDB )
                                       BB( I + ( J - 1 )*LDB ) = ALPHA*
     $                                    B( I, J )
   60                               CONTINUE
   70                            CONTINUE
C
                                 IF( LEFT )THEN
                                    CALL CMMCH( TRANSA, 'N', M, N, M,
     $                                          ONE, A, NMAX, C, NMAX,
     $                                          ZERO, B, NMAX, CT, G,
     $                                          BB, LDB, EPS, ERR,
     $                                          FTL, NOUT, .FALSE.,
     $                                          KPRINT )
                                 ELSE
                                    CALL CMMCH( 'N', TRANSA, M, N, N,
     $                                          ONE, C, NMAX, A, NMAX,
     $                                          ZERO, B, NMAX, CT, G,
     $                                          BB, LDB, EPS, ERR,
     $                                          FTL, NOUT, .FALSE.,
     $                                          KPRINT )
                                 END IF
                              END IF
                              ERRMAX = MAX( ERRMAX, ERR )
                           END IF
                           IF (FTL) THEN
                             FATAL = .TRUE.
                             IF (KPRINT .GE. 3) THEN
                               WRITE( NOUT, FMT = 9996 )SNAME
                               WRITE( NOUT, FMT = 9995 )NC,
     $                              SNAME, SIDE, UPLO, TRANSA,
     $                              DIAG, M, N, ALPHA, LDA, LDB
                             ENDIF
                           ENDIF
C
   80                   CONTINUE
C
   90                CONTINUE
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
      IF (.NOT. FATAL) THEN
        IF (KPRINT .GE. 3 ) THEN
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
 9995 FORMAT( 1X, I6, ': ', A6, '(', 4( '''', A1, ''',' ), 2( I3, ',' ),
     $      '(', F4.1, ',', F4.1, '), A,', I3, ', B,', I3, ')         ',
     $      '      .' )
 9994 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',
     $      '******' )
C
C     End of CCHK33.
C
      END
