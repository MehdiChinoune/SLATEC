*DECK DCHK53
      SUBROUTINE DCHK53 (SNAME, EPS, THRESH, NOUT, KPRINT, FATAL, NIDIM,
     $   IDIM, NALF, ALF, NBET, BET, NMAX, AB, AA, AS, BB, BS, C, CC,
     $   CS, CT, G, W)
C***BEGIN PROLOGUE  DCHK53
C***SUBSIDIARY
C***PURPOSE  Test DSYR2K.
C***LIBRARY   SLATEC (BLAS)
C***KEYWORDS  BLAS, QUICK CHECK SERVICE ROUTINE
C***AUTHOR  Dongarra, J. J., (ANL)
C           Duff, I., (AERE)
C           Du Croz, J., (NAG)
C           Hammarling, S., (NAG)
C***DESCRIPTION
C
C  Quick check for DSYR2K.
C
C  Auxiliary routine for test program for Level 3 Blas.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DMAKE3, DMMCH, DSYR2K, LDE, LDERES, NUMXER
C***REVISION HISTORY  (YYMMDD)
C   890208  DATE WRITTEN
C   910619  Modified to meet SLATEC code and prologue standards. (BKS)
C***END PROLOGUE  DCHK53
C     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
C     .. Scalar Arguments ..
      LOGICAL            FATAL
      DOUBLE PRECISION   EPS, THRESH
      INTEGER            KPRINT, NALF, NBET, NIDIM, NMAX, NOUT
      CHARACTER*6        SNAME
C     .. Array Arguments ..
      DOUBLE PRECISION   AA( NMAX*NMAX ), AB( 2*NMAX*NMAX ),
     $                   ALF( NALF ), AS( NMAX*NMAX ), BB( NMAX*NMAX ),
     $                   BET( NBET ), BS( NMAX*NMAX ), C( NMAX, NMAX ),
     $                   CC( NMAX*NMAX ), CS( NMAX*NMAX ), CT( NMAX ),
     $                   G( NMAX ), W( 2*NMAX )
      INTEGER            IDIM( NIDIM )
C     .. Local Scalars ..
      DOUBLE PRECISION   ALPHA, ALS, BETA, BETS, ERR, ERRMAX
      INTEGER            I, IA, IB, ICT, ICU, IK, IN, J, JC, JJ, JJAB,
     $                   K, KS, LAA, LBB, LCC, LDA, LDAS, LDB, LDBS,
     $                   LDC, LDCS, LJ, MA, N, NA, NARGS, NC, NERR, NS
      LOGICAL            FTL, NULL, RESET, TRAN, UPPER
      CHARACTER*1        TRANS, TRANSS, UPLO, UPLOS
      CHARACTER*2        ICHU
      CHARACTER*3        ICHT
C     .. Local Arrays ..
      LOGICAL            ISAME( 13 )
C     .. External Functions ..
      INTEGER            NUMXER
      LOGICAL            LDE, LDERES
      EXTERNAL           LDE, LDERES, NUMXER
C     .. External Subroutines ..
      EXTERNAL           DSYR2K, DMAKE3, DMMCH
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
C     .. Data statements ..
      DATA               ICHU/'UL'/, ICHT/'NTC'/
C***FIRST EXECUTABLE STATEMENT  DCHK53
      NARGS = 12
      NC = 0
      RESET = .TRUE.
      ERRMAX = ZERO
C
      DO 130 IN = 1, NIDIM
         N = IDIM( IN )
C        Set LDC to 1 more than minimum value if room.
         LDC = N
         IF( LDC.LT.NMAX )
     $      LDC = LDC + 1
C        Skip tests if not enough room.
         IF( LDC.GT.NMAX )
     $      GO TO 130
         LCC = LDC*N
         NULL = N.LE.0
C
         DO 120 IK = 1, NIDIM
            K = IDIM( IK )
C
            DO 110 ICT = 1, 3
               TRANS = ICHT( ICT: ICT )
               TRAN = TRANS.EQ.'T'.OR.TRANS.EQ.'C'
               IF( TRAN )THEN
                  MA = K
                  NA = N
               ELSE
                  MA = N
                  NA = K
               END IF
C              Set LDA to 1 more than minimum value if room.
               LDA = MA
               IF( LDA.LT.NMAX )
     $            LDA = LDA + 1
C              Skip tests if not enough room.
               IF( LDA.GT.NMAX )
     $            GO TO 110
               LAA = LDA*NA
C
C              Generate the matrix A.
C
               IF( TRAN )THEN
                  CALL DMAKE3( 'GE', ' ', ' ', MA, NA, AB, 2*NMAX, AA,
     $                        LDA, RESET, ZERO )
               ELSE
                  CALL DMAKE3('GE', ' ', ' ', MA, NA, AB, NMAX, AA, LDA,
     $                        RESET, ZERO )
               END IF
C
C              Generate the matrix B.
C
               LDB = LDA
               LBB = LAA
               IF( TRAN )THEN
                  CALL DMAKE3( 'GE', ' ', ' ', MA, NA, AB( K + 1 ),
     $                        2*NMAX, BB, LDB, RESET, ZERO )
               ELSE
                  CALL DMAKE3( 'GE', ' ', ' ', MA, NA, AB( K*NMAX + 1 ),
     $                        NMAX, BB, LDB, RESET, ZERO )
               END IF
C
               DO 100 ICU = 1, 2
                  UPLO = ICHU( ICU: ICU )
                  UPPER = UPLO.EQ.'U'
C
                  DO 90 IA = 1, NALF
                     ALPHA = ALF( IA )
C
                     DO 80 IB = 1, NBET
                        BETA = BET( IB )
C
C                       Generate the matrix C.
C
                        CALL DMAKE3( 'SY', UPLO, ' ', N, N, C, NMAX, CC,
     $                              LDC, RESET, ZERO )
C
                        NC = NC + 1
C
C                       Save every datum before calling the subroutine.
C
                        UPLOS = UPLO
                        TRANSS = TRANS
                        NS = N
                        KS = K
                        ALS = ALPHA
                        DO 10 I = 1, LAA
                           AS( I ) = AA( I )
   10                   CONTINUE
                        LDAS = LDA
                        DO 20 I = 1, LBB
                           BS( I ) = BB( I )
   20                   CONTINUE
                        LDBS = LDB
                        BETS = BETA
                        DO 30 I = 1, LCC
                           CS( I ) = CC( I )
   30                   CONTINUE
                        LDCS = LDC
C
C                       Call the subroutine.
C
                        CALL DSYR2K( UPLO, TRANS, N, K, ALPHA, AA, LDA,
     $                               BB, LDB, BETA, CC, LDC )
C
C                       Check if error-exit was taken incorrectly.
C
                        IF( NUMXER(NERR) .NE. 0 )THEN
                          IF (KPRINT .GE. 2) THEN
                           WRITE( NOUT, FMT = 9993 )
                           END IF
                          FATAL = .TRUE.
                        END IF
C
C                       See what data changed inside subroutines.
C
                        ISAME( 1 ) = UPLOS.EQ.UPLO
                        ISAME( 2 ) = TRANSS.EQ.TRANS
                        ISAME( 3 ) = NS.EQ.N
                        ISAME( 4 ) = KS.EQ.K
                        ISAME( 5 ) = ALS.EQ.ALPHA
                        ISAME( 6 ) = LDE( AS, AA, LAA )
                        ISAME( 7 ) = LDAS.EQ.LDA
                        ISAME( 8 ) = LDE( BS, BB, LBB )
                        ISAME( 9 ) = LDBS.EQ.LDB
                        ISAME( 10 ) = BETS.EQ.BETA
                        IF( NULL )THEN
                           ISAME( 11 ) = LDE( CS, CC, LCC )
                        ELSE
                           ISAME( 11 ) = LDERES( 'SY', UPLO, N, N, CS,
     $                                   CC, LDC )
                        END IF
                        ISAME( 12 ) = LDCS.EQ.LDC
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
C
C                          Check the result column by column.
C
                           JJAB = 1
                           JC = 1
                           DO 70 J = 1, N
                              IF( UPPER )THEN
                                 JJ = 1
                                 LJ = J
                              ELSE
                                 JJ = J
                                 LJ = N - J + 1
                              END IF
                              IF( TRAN )THEN
                                 DO 50 I = 1, K
                                    W( I ) = AB( ( J - 1 )*2*NMAX + K +
     $                                       I )
                                    W( K + I ) = AB( ( J - 1 )*2*NMAX +
     $                                           I )
   50                            CONTINUE
                                 CALL DMMCH( 'T', 'N', LJ, 1, 2*K,
     $                                       ALPHA, AB( JJAB ), 2*NMAX,
     $                                       W, 2*NMAX, BETA,
     $                                       C( JJ, J ), NMAX, CT, G,
     $                                       CC( JC ), LDC, EPS, ERR,
     $                                       FTL, NOUT, .TRUE.,
     $                                       KPRINT )
                              ELSE
                                 DO 60 I = 1, K
                                    W( I ) = AB( ( K + I - 1 )*NMAX +
     $                                       J )
                                    W( K + I ) = AB( ( I - 1 )*NMAX +
     $                                           J )
   60                            CONTINUE
                                 CALL DMMCH( 'N', 'N', LJ, 1, 2*K,
     $                                       ALPHA, AB( JJ ), NMAX, W,
     $                                       2*NMAX, BETA, C( JJ, J ),
     $                                       NMAX, CT, G, CC( JC ), LDC,
     $                                       EPS, ERR, FTL, NOUT,
     $                                       .TRUE., KPRINT )
                              END IF
                              IF( UPPER )THEN
                                 JC = JC + LDC
                              ELSE
                                 JC = JC + LDC + 1
                                 IF( TRAN )
     $                              JJAB = JJAB + 2*NMAX
                              END IF
                              ERRMAX = MAX( ERRMAX, ERR )
   70                      CONTINUE
                        END IF
                        IF (FTL) THEN
                          FATAL = .TRUE.
                          IF (KPRINT .GE. 3) THEN
                            WRITE( NOUT, FMT = 9996 )SNAME
                            WRITE( NOUT, FMT = 9994 )NC,
     $                           SNAME, UPLO, TRANS,
     $                           N, K, ALPHA, LDA, LDB,
     $                           BETA, LDC
                          ENDIF
                        ENDIF
C
   80                CONTINUE
C
   90             CONTINUE
C
  100          CONTINUE
C
  110       CONTINUE
C
  120    CONTINUE
C
  130 CONTINUE
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
 9994 FORMAT( 1X, I6, ': ', A6, '(', 2( '''', A1, ''',' ), 2( I3, ',' ),
     $      F4.1, ', A,', I3, ', B,', I3, ',', F4.1, ', C,', I3, ')   ',
     $      ' .' )
 9993 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',
     $      '******' )
C
C     End of DCHK53.
C
      END
