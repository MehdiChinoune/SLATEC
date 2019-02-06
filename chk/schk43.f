*DECK SCHK43
      SUBROUTINE SCHK43 (SNAME, EPS, THRESH, NOUT, KPRINT, FATAL, NIDIM,
     $   IDIM, NALF, ALF, NBET, BET, NMAX, A, AA, AS, B, BB, BS, C, CC,
     $   CS, CT, G)
C***BEGIN PROLOGUE  SCHK43
C***SUBSIDIARY
C***PURPOSE  Quick check for SSYRK.
C***LIBRARY   SLATEC (BLAS)
C***KEYWORDS  BLAS, QUICK CHECK SERVICE ROUTINE
C***AUTHOR  Dongarra, J. J., (ANL)
C           Duff, I., (AERE)
C           Du Croz, J., (NAG)
C           Hammarling, S., (NAG)
C***DESCRIPTION
C
C  Quick check for SSYRK.
C
C  Auxiliary routine for test program for Level 3 Blas.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  LSE, LSERES, NUMXER, SMAKE3, SMMCH, SSYRK
C***REVISION HISTORY  (YYMMDD)
C   890208  DATE WRITTEN
C   910619  Modified to meet SLATEC code and prologue standards.  (BKS)
C***END PROLOGUE  SCHK43
C     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0, ONE = 1.0)
C     .. Scalar Arguments ..
      LOGICAL            FATAL
      REAL               EPS, THRESH
      INTEGER            KPRINT, NALF, NBET, NIDIM, NMAX, NOUT
      CHARACTER*6        SNAME
C     .. Array Arguments ..
      REAL               A( NMAX, NMAX), AA( NMAX*NMAX), ALF( NALF),
     $                   AS( NMAX*NMAX ), BET( NBET ), G( NMAX ),
     $                   B( NMAX, NMAX ), BB( NMAX*NMAX ),
     $                   BS( NMAX*NMAX ), C( NMAX, NMAX ),
     $                   CC( NMAX*NMAX ), CS( NMAX*NMAX ),
     $                   CT( NMAX )
      INTEGER            IDIM( NIDIM )
C     .. Local Scalars ..
      REAL               ALPHA, ALS, BETA, BETS, ERR, ERRMAX
      INTEGER            I, IA, IB, ICT, ICU, IK, IN, J, JC, JJ, K, LAA,
     $                   LCC, LDA, LDAS, LDC, LDCS,
     $                   N, NA, NARGS, NC, NERR, NS, KS, LJ, MA
      LOGICAL            FTL, NULL, RESET, TRAN, UPPER
      CHARACTER*1        UPLO, UPLOS, TRANS, TRANSS
      CHARACTER*2        ICHU
      CHARACTER*3        ICHT
C     .. Local Arrays ..
      LOGICAL            ISAME( 13 )
C     .. External Functions ..
      INTEGER            NUMXER
      LOGICAL            LSE, LSERES
      EXTERNAL           LSE, LSERES, NUMXER
C     .. External Subroutines ..
      EXTERNAL           SSYRK, SMAKE3, SMMCH
C     .. Intrinsic Functions ..
      INTRINSIC          MAX
C     .. Data statements ..
      DATA               ICHT/'NTC'/, ICHU/'UL'/
C***FIRST EXECUTABLE STATEMENT  SCHK43
      NARGS = 10
C
      NC = 0
      RESET = .TRUE.
      ERRMAX = ZERO
C
C
      DO 100 IN = 1, NIDIM
         N = IDIM( IN )
C        Set LDC to 1 more than minimum value if room.
         LDC = N
         IF( LDC.LT.NMAX )
     $      LDC = LDC + 1
C        Skip tests if not enough room.
         IF( LDC.GT.NMAX )
     $      GO TO 100
         LCC = LDC*N
         NULL = N.LE.0
C
         DO 90 IK = 1, NIDIM
            K = IDIM( IK )
C
            DO 80 ICT = 1, 3
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
     $            GO TO 80
               LAA = LDA*NA
C
C              Generate the matrix A.
C
               CALL SMAKE3( 'GE', ' ', ' ', MA, NA, A, NMAX, AA, LDA,
     $                     RESET, ZERO )
C
               DO 70 ICU = 1, 2
                  UPLO = ICHU( ICU: ICU )
                  UPPER = UPLO.EQ.'U'
C
                  DO 60 IA = 1, NALF
                     ALPHA = ALF( IA )
C
                     DO 50 IB = 1, NBET
                        BETA = BET( IB )
C
C                       Generate the matrix C.
C
                        CALL SMAKE3( 'SY', UPLO, ' ', N, N, C, NMAX, CC,
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
                        BETS = BETA
                        DO 20 I = 1, LCC
                           CS( I ) = CC( I )
   20                   CONTINUE
                        LDCS = LDC
C
C                       Call the subroutine.
C
                        CALL SSYRK( UPLO, TRANS, N, K, ALPHA, AA, LDA,
     $                              BETA, CC, LDC )
C
C                       Check if error-exit was taken incorrectly.
C
                        IF(NUMXER(NERR) .NE. 0) THEN
                           IF (KPRINT .GE. 2) THEN
                           WRITE( NOUT, FMT = 9993 )
                           ENDIF
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
                        ISAME( 6 ) = LSE( AS, AA, LAA )
                        ISAME( 7 ) = LDAS.EQ.LDA
                        ISAME( 8 ) = BETS.EQ.BETA
                        IF( NULL )THEN
                           ISAME( 9 ) = LSE( CS, CC, LCC )
                        ELSE
                           ISAME( 9 ) = LSERES( 'SY', UPLO, N, N, CS,
     $                                  CC, LDC )
                        END IF
                        ISAME( 10 ) = LDCS.EQ.LDC
C
C                       If data was incorrectly changed, report and
C                       return.
C
                        DO 30 I = 1, NARGS
                          IF (.NOT. ISAME( I )) THEN
                            FATAL = .TRUE.
                            IF (KPRINT .GE. 2) THEN
                              WRITE( NOUT, FMT = 9998 )I
                            ENDIF
                          ENDIF
  30                    CONTINUE
C
                        FTL = .FALSE.
                        IF( .NOT.NULL )THEN
C
C                          Check the result column by column.
C
                           JC = 1
                           DO 40 J = 1, N
                              IF( UPPER )THEN
                                 JJ = 1
                                 LJ = J
                              ELSE
                                 JJ = J
                                 LJ = N - J + 1
                              END IF
                              IF( TRAN )THEN
                                 CALL SMMCH( 'T', 'N', LJ, 1, K, ALPHA,
     $                                       A( 1, JJ ), NMAX,
     $                                       A( 1, J ), NMAX, BETA,
     $                                       C( JJ, J ), NMAX, CT, G,
     $                                       CC( JC ), LDC, EPS, ERR,
     $                                       FTL, NOUT, .TRUE.,KPRINT)
                              ELSE
                                 CALL SMMCH( 'N', 'T', LJ, 1, K, ALPHA,
     $                                       A( JJ, 1 ), NMAX,
     $                                       A( J, 1 ), NMAX, BETA,
     $                                       C( JJ, J ), NMAX, CT, G,
     $                                       CC( JC ), LDC, EPS, ERR,
     $                                       FTL, NOUT, .TRUE.,KPRINT)
                              END IF
                              IF( UPPER )THEN
                                 JC = JC + LDC
                              ELSE
                                 JC = JC + LDC + 1
                              END IF
                              ERRMAX = MAX( ERRMAX, ERR )
   40                      CONTINUE
                        END IF
                        IF (FTL) THEN
                          FATAL = .TRUE.
                          IF (KPRINT .GE. 3) THEN
                            WRITE (NOUT, FMT = 9996) SNAME
                            WRITE( NOUT, FMT = 9994 )NC, SNAME,
     $                        UPLO, TRANS, N, K,
     $                        ALPHA, LDA, BETA, LDC
                          ENDIF
                        ENDIF
C
   50                CONTINUE
C
   60             CONTINUE
C
   70          CONTINUE
C
   80       CONTINUE
C
   90    CONTINUE
C
  100 CONTINUE
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
 9994 FORMAT( 1X, I6, ': ', A6, '(', 2( '''', A1, ''',' ), 2( I3, ',' ),
     $      F4.1, ', A,', I3, ',', F4.1, ', C,', I3, ')           .' )
 9993 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',
     $      '******' )
C
C     End of SCHK43.
C
      END
