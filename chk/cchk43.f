*DECK CCHK43
      SUBROUTINE CCHK43 (SNAME, EPS, THRESH, NOUT, KPRINT, FATAL, NIDIM,
     $   IDIM, NALF, ALF, NBET, BET, NMAX, A, AA, AS, B, BB, BS, C, CC,
     $   CS, CT, G)
C***BEGIN PROLOGUE  CCHK43
C***SUBSIDIARY
C***PURPOSE  Quick check for CHERK and CSYRK.
C***LIBRARY   SLATEC (BLAS)
C***KEYWORDS  BLAS, QUICK CHECK SERVICE ROUTINE
C***AUTHOR  Dongarra, J. J., (ANL)
C           Duff, I., (AERE)
C           Du Croz, J., (NAG)
C           Hammarling, S., (NAG)
C***DESCRIPTION
C
C  Quick check for CHERK and CSYRK.
C
C  Auxiliary routine for test program for Level 3 Blas.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CHERK, CMAKE3, CMMCH, CSYRK, LCE, LCERES, NUMXER
C***REVISION HISTORY  (YYMMDD)
C   890208  DATE WRITTEN
C   910619  Modified to meet SLATEC code and prologue standards.  (BKS)
C***END PROLOGUE  CCHK43
C     .. Parameters ..
      COMPLEX            ZERO
      PARAMETER          ( ZERO = ( 0.0, 0.0 ))
      REAL               RZERO, RONE
      PARAMETER          ( RZERO = 0.0, RONE = 1.0)
C     .. Scalar Arguments ..
      LOGICAL            FATAL
      REAL               EPS, THRESH
      INTEGER            KPRINT, NALF, NBET, NIDIM, NMAX, NOUT
      CHARACTER*6        SNAME
C     .. Array Arguments ..
      COMPLEX            A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ),
     $                   AS( NMAX*NMAX ), B( NMAX, NMAX ),
     $                   BB( NMAX*NMAX ), BS( NMAX*NMAX ), BET(NBET),
     $                   C( NMAX, NMAX ), CC( NMAX*NMAX ),
     $                   CT( NMAX ), CS( NMAX*NMAX )
      REAL               G( NMAX )
      INTEGER            IDIM( NIDIM )
C     .. Local Scalars ..
      COMPLEX            ALPHA, ALS, BETA, BETS
      REAL               ERR, ERRMAX, RALPHA, RALS, RBETA, RBETS
      INTEGER            I, IA, IB, ICT, ICU, IK, IN, J, JC, JJ, K, KS,
     $                   LAA, LCC, LDA, LDAS, LDC, LDCS, LJ, MA, N, NA,
     $                   NARGS, NC, NERR, NS
      LOGICAL            CONJ, FTL, NULL, RESET, TRAN, UPPER
      CHARACTER*1        TRANS, TRANSS, UPLO, TRANST, UPLOS
      CHARACTER*2        ICHU, ICHT
C     .. Local Arrays ..
      LOGICAL            ISAME( 13 )
C     .. External Functions ..
      INTEGER            NUMXER
      LOGICAL            LCE, LCERES
      EXTERNAL           LCE, LCERES, NUMXER
C     .. External Subroutines ..
      EXTERNAL           CHERK, CSYRK, CMAKE3, CMMCH
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
C     .. Data statements ..
      DATA               ICHU/'UL'/, ICHT/'NC'/
C***FIRST EXECUTABLE STATEMENT  CCHK43
      CONJ = SNAME( 2: 3 ).EQ.'HE'
C
      NARGS = 10
      NC = 0
      RESET = .TRUE.
      ERRMAX = RZERO
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
C
         DO 90 IK = 1, NIDIM
            K = IDIM( IK )
C
            DO 80 ICT = 1, 2
               TRANS = ICHT( ICT: ICT )
               TRAN = TRANS.EQ.'C'
               IF( TRAN.AND..NOT.CONJ )
     $            TRANS = 'T'
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
               CALL CMAKE3( 'GE', ' ', ' ', MA, NA, A, NMAX, AA, LDA,
     $                     RESET, ZERO )
C
               DO 70 ICU = 1, 2
                  UPLO = ICHU( ICU: ICU )
                  UPPER = UPLO.EQ.'U'
C
                  DO 60 IA = 1, NALF
                     ALPHA = ALF( IA )
                     IF( CONJ )THEN
                        RALPHA = REAL( ALPHA )
                        ALPHA = CMPLX( RALPHA, RZERO )
                     END IF
C
                     DO 50 IB = 1, NBET
                        BETA = BET( IB )
                        IF( CONJ )THEN
                           RBETA = REAL( BETA )
                           BETA = CMPLX( RBETA, RZERO )
                        END IF
                        NULL = N.LE.0
                        IF( CONJ )
     $                     NULL = NULL.OR.( ( K.LE.0.OR.RALPHA.EQ.
     $                            RZERO ).AND.RBETA.EQ.RONE )
C
C                       Generate the matrix C.
C
                        CALL CMAKE3( SNAME( 2: 3 ), UPLO, ' ', N, N, C,
     $                              NMAX, CC, LDC, RESET, ZERO )
C
                        NC = NC + 1
C
C                       Save every datum before calling the subroutine.
C
                        UPLOS = UPLO
                        TRANSS = TRANS
                        NS = N
                        KS = K
                        IF( CONJ )THEN
                           RALS = RALPHA
                        ELSE
                           ALS = ALPHA
                        END IF
                        DO 10 I = 1, LAA
                           AS( I ) = AA( I )
   10                   CONTINUE
                        LDAS = LDA
                        IF( CONJ )THEN
                           RBETS = RBETA
                        ELSE
                           BETS = BETA
                        END IF
                        DO 20 I = 1, LCC
                           CS( I ) = CC( I )
   20                   CONTINUE
                        LDCS = LDC
C
C                       Call the subroutine.
C
                        IF( CONJ )THEN
                           CALL CHERK( UPLO, TRANS, N, K, RALPHA, AA,
     $                                 LDA, RBETA, CC, LDC )
                        ELSE
                           CALL CSYRK( UPLO, TRANS, N, K, ALPHA, AA,
     $                                 LDA, BETA, CC, LDC )
                        END IF
C
C                       Check if error-exit was taken incorrectly.
C
                        IF( NUMXER(NERR) .NE. 0 )THEN
                          IF (KPRINT .GE. 2) THEN
                           WRITE( NOUT, FMT = 9992 )
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
                        IF( CONJ )THEN
                           ISAME( 5 ) = RALS.EQ.RALPHA
                        ELSE
                           ISAME( 5 ) = ALS.EQ.ALPHA
                        END IF
                        ISAME( 6 ) = LCE( AS, AA, LAA )
                        ISAME( 7 ) = LDAS.EQ.LDA
                        IF( CONJ )THEN
                           ISAME( 8 ) = RBETS.EQ.RBETA
                        ELSE
                           ISAME( 8 ) = BETS.EQ.BETA
                        END IF
                        IF( NULL )THEN
                           ISAME( 9 ) = LCE( CS, CC, LCC )
                        ELSE
                           ISAME( 9 ) = LCERES( SNAME( 2: 3 ), UPLO, N,
     $                                  N, CS, CC, LDC )
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
                        IF( .NOT.NULL )THEN
C
C                          Check the result column by column.
C
                           IF( CONJ )THEN
                              TRANST = 'C'
                           ELSE
                              TRANST = 'T'
                           END IF
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
                                FTL = .FALSE.
                                 CALL CMMCH( TRANST, 'N', LJ, 1, K,
     $                                       ALPHA, A( 1, JJ ), NMAX,
     $                                       A( 1, J ), NMAX, BETA,
     $                                       C( JJ, J ), NMAX, CT, G,
     $                                       CC( JC ), LDC, EPS, ERR,
     $                                       FTL, NOUT, .TRUE.,
     $                                       KPRINT )
                              ELSE
                                FTL = .FALSE.
                                 CALL CMMCH( 'N', TRANST, LJ, 1, K,
     $                                       ALPHA, A( JJ, 1 ), NMAX,
     $                                       A( J, 1 ), NMAX, BETA,
     $                                       C( JJ, J ), NMAX, CT, G,
     $                                       CC( JC ), LDC, EPS, ERR,
     $                                       FTL, NOUT, .TRUE.,
     $                                       KPRINT )
                              END IF
                              IF( UPPER )THEN
                                 JC = JC + LDC
                              ELSE
                                 JC = JC + LDC + 1
                              END IF
                              ERRMAX = MAX( ERRMAX, ERR )
                           IF (FTL) THEN
                           FATAL = .TRUE.
                           IF (KPRINT .GE. 3) THEN
                             WRITE( NOUT, FMT = 9996 )SNAME
                             IF( CONJ )THEN
                                WRITE( NOUT, FMT = 9994 )NC, SNAME,
     $                             UPLO, TRANS, N, K, RALPHA,
     $                             LDA, RBETA, LDC
                             ELSE
                                WRITE( NOUT, FMT = 9993 )NC, SNAME,
     $                             UPLO, TRANS, N, K, ALPHA,
     $                             LDA, BETA, LDC
                             ENDIF
                           ENDIF
                           ENDIF
   40                      CONTINUE
                        END IF
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
 9994 FORMAT( 1X, I6, ': ', A6, '(', 2( '''', A1, ''',' ), 2( I3, ',' ),
     $      F4.1, ', A,', I3, ',', F4.1, ', C,', I3, ')               ',
     $      '          .' )
 9993 FORMAT( 1X, I6, ': ', A6, '(', 2( '''', A1, ''',' ), 2( I3, ',' ),
     $      '(', F4.1, ',', F4.1, ') , A,', I3, ',(', F4.1, ',', F4.1,
     $      '), C,', I3, ')          .' )
 9992 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',
     $      '******' )
C
C     End of CCHK43.
C
      END
