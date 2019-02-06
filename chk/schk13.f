*DECK SCHK13
      SUBROUTINE SCHK13 (SNAME, EPS, THRESH, NOUT, KPRINT, FATAL, NIDIM,
     $   IDIM, NALF, ALF, NBET, BET, NMAX, A, AA, AS, B, BB, BS, C, CC,
     $   CS, CT, G)
C***BEGIN PROLOGUE  SCHK13
C***SUBSIDIARY
C***PURPOSE  Quick check for SGEMM.
C***LIBRARY   SLATEC (BLAS)
C***KEYWORDS  BLAS, QUICK CHECK SERVICE ROUTINE
C***AUTHOR  Dongarra, J. J., (ANL)
C           Duff, I., (AERE)
C           Du Croz, J., (NAG)
C           Hammarling, S., (NAG)
C***DESCRIPTION
C
C  Quick check for SGEMM.
C
C  Auxiliary routine for test program for Level 3 Blas.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  LSE, LSERES, NUMXER, SGEMM, SMAKE3, SMMCH
C***REVISION HISTORY  (YYMMDD)
C   890208  DATE WRITTEN
C   910619  Modified to meet SLATEC code and prologue standards.  (BKS)
C***END PROLOGUE  SCHK13
C     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0 )
C     .. Scalar Arguments ..
      LOGICAL            FATAL
      REAL               EPS, THRESH
      INTEGER            KPRINT, NALF, NBET, NIDIM, NMAX, NOUT
      CHARACTER*6        SNAME
C     .. Array Arguments ..
      REAL               A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ),
     $                   AS( NMAX*NMAX ), BET( NBET ), G( NMAX ),
     $                   B( NMAX, NMAX ), BB( NMAX*NMAX ),
     $                   BS( NMAX*NMAX ), C( NMAX, NMAX ),
     $                   CC( NMAX*NMAX ), CS( NMAX*NMAX ),
     $                   CT( NMAX )
      INTEGER            IDIM( NIDIM )
C     .. Local Scalars ..
      REAL               ALPHA, ALS, BETA, BLS, ERR, ERRMAX
      INTEGER            I, IA, IB, ICA, ICB, IK, IM, IN, K, KS, LAA,
     $                   LBB, LCC, LDA, LDAS, LDB, LDBS, LDC, LDCS, M,
     $                   MA, MB, MS, N, NA, NARGS, NB, NC, NERR, NS
      LOGICAL            FTL, NULL, RESET, TRANA, TRANB
      CHARACTER*1        TRANAS, TRANBS, TRANSA, TRANSB
      CHARACTER*3        ICH
C     .. Local Arrays ..
      LOGICAL            ISAME( 13 )
C     .. External Functions ..
      INTEGER            NUMXER
      LOGICAL            LSE, LSERES
      EXTERNAL           LSE, LSERES, NUMXER
C     .. External Subroutines ..
      EXTERNAL           SGEMM, SMAKE3, SMMCH
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     .. Data statements ..
      DATA               ICH/'NTC'/
C***FIRST EXECUTABLE STATEMENT  SCHK13
      NARGS = 13
C
      NC = 0
      RESET = .TRUE.
      ERRMAX = ZERO
C
C
      DO 110 IM = 1, NIDIM
         M = IDIM( IM )
C
         DO 100 IN = 1, NIDIM
            N = IDIM( IN )
C           Set LDC to 1 more than minimum value if room.
            LDC = M
            IF( LDC.LT.NMAX )
     $         LDC = LDC + 1
C           Skip tests if not enough room.
            IF( LDC.GT.NMAX )
     $         GO TO 100
            LCC = LDC*N
            NULL = N.LE.0.OR.M.LE.0
C
            DO 90 IK = 1, NIDIM
               K = IDIM( IK )
C
               DO 80 ICA = 1, 3
                  TRANSA = ICH( ICA: ICA )
                  TRANA = TRANSA.EQ.'T'.OR.TRANSA.EQ.'C'
C
                  IF( TRANA )THEN
                     MA = K
                     NA = M
                  ELSE
                     MA = M
                     NA = K
                  END IF
C                 Set LDA to 1 more than minimum value if room.
                  LDA = MA
                  IF( LDA.LT.NMAX )
     $               LDA = LDA + 1
C                 Skip tests if not enough room.
                  IF( LDA.GT.NMAX )
     $               GO TO 80
                  LAA = LDA*NA
C
C                 Generate the matrix A.
C
                  CALL SMAKE3( 'GE', ' ', ' ', MA, NA, A, NMAX, AA, LDA,
     $                        RESET, ZERO )
C
                  DO 70 ICB = 1, 3
                     TRANSB = ICH( ICB: ICB )
                     TRANB = TRANSB.EQ.'T'.OR.TRANSB.EQ.'C'
C
                     IF( TRANB )THEN
                        MB = N
                        NB = K
                     ELSE
                        MB = K
                        NB = N
                     END IF
C                    Set LDB to 1 more than minimum value if room.
                     LDB = MB
                     IF( LDB.LT.NMAX )
     $                  LDB = LDB + 1
C                    Skip tests if not enough room.
                     IF( LDB.GT.NMAX )
     $                  GO TO 70
                     LBB = LDB*NB
C
C                    Generate the matrix B.
C
                     CALL SMAKE3( 'GE', ' ', ' ', MB, NB, B, NMAX, BB,
     $                           LDB, RESET, ZERO )
C
                     DO 60 IA = 1, NALF
                        ALPHA = ALF( IA )
C
                        DO 50 IB = 1, NBET
                           BETA = BET( IB )
C
C                          Generate the matrix C.
C
                           CALL SMAKE3( 'GE', ' ', ' ', M, N, C, NMAX,
     $                                 CC, LDC, RESET, ZERO )
C
                           NC = NC + 1
C
C                          Save every datum before calling the
C                          subroutine.
C
                           TRANAS = TRANSA
                           TRANBS = TRANSB
                           MS = M
                           NS = N
                           KS = K
                           ALS = ALPHA
                           DO 10 I = 1, LAA
                              AS( I ) = AA( I )
   10                      CONTINUE
                           LDAS = LDA
                           DO 20 I = 1, LBB
                              BS( I ) = BB( I )
   20                      CONTINUE
                           LDBS = LDB
                           BLS = BETA
                           DO 30 I = 1, LCC
                              CS( I ) = CC( I )
   30                      CONTINUE
                           LDCS = LDC
C
C                          Call the subroutine.
C
                           CALL SGEMM( TRANSA, TRANSB, M, N, K, ALPHA,
     $                                 AA, LDA, BB, LDB, BETA, CC, LDC )
C
C                          Check if error-exit was taken incorrectly.
C
                           IF( NUMXER(NERR) .NE. 0) THEN
                            IF (KPRINT .GE. 2) THEN
                              WRITE( NOUT, FMT = 9994 )
                            ENDIF
                            FATAL = .TRUE.
                           END IF
C
C                          See what data changed inside subroutines.
C
                           ISAME( 1 ) = TRANSA.EQ.TRANAS
                           ISAME( 2 ) = TRANSB.EQ.TRANBS
                           ISAME( 3 ) = MS.EQ.M
                           ISAME( 4 ) = NS.EQ.N
                           ISAME( 5 ) = KS.EQ.K
                           ISAME( 6 ) = ALS.EQ.ALPHA
                           ISAME( 7 ) = LSE( AS, AA, LAA )
                           ISAME( 8 ) = LDAS.EQ.LDA
                           ISAME( 9 ) = LSE( BS, BB, LBB )
                           ISAME( 10 ) = LDBS.EQ.LDB
                           ISAME( 11 ) = BLS.EQ.BETA
                           IF( NULL )THEN
                              ISAME( 12 ) = LSE( CS, CC, LCC )
                           ELSE
                              ISAME( 12 ) = LSERES( 'GE', ' ', M, N, CS,
     $                                      CC, LDC )
                           END IF
                           ISAME( 13 ) = LDCS.EQ.LDC
C
C                          If data was incorrectly changed, report
C                          and return.
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
C                             Check the result.
C
                              CALL SMMCH( TRANSA, TRANSB, M, N, K,
     $                                    ALPHA, A, NMAX, B, NMAX, BETA,
     $                                    C, NMAX, CT, G, CC, LDC, EPS,
     $                                    ERR, FTL, NOUT, .TRUE.,
     $                                    KPRINT )
                              ERRMAX = MAX( ERRMAX, ERR )
                           END IF
                           IF (FTL) THEN
                             FATAL = .TRUE.
                             IF (KPRINT .GE. 3) THEN
                               WRITE (NOUT, FMT = 9996) SNAME
                               WRITE( NOUT, FMT = 9995 )NC, SNAME,
     $                           TRANSA, TRANSB, M, N, K,
     $                           ALPHA, LDA, LDB, BETA, LDC
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
 9995 FORMAT( 1X, I6, ': ', A6, '(''', A1, ''',''', A1, ''',',
     $      3( I3, ',' ), F4.1, ', A,', I3, ', B,', I3, ',', F4.1, ', ',
     $      'C,', I3, ').' )
 9994 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',
     $      '******' )
C
C     End of SCHK13.
C
      END
