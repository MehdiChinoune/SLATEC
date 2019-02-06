*DECK CCHK53
      SUBROUTINE CCHK53 (SNAME, EPS, THRESH, NOUT, KPRINT, FATAL, NIDIM,
     $   IDIM, NALF, ALF, NBET, BET, NMAX, AB, AA, AS, BB, BS, C, CC,
     $   CS, CT, G, W)
C***BEGIN PROLOGUE  CCHK53
C***SUBSIDIARY
C***PURPOSE  Quick check for CHER2K and CSYR2K.
C***LIBRARY   SLATEC (BLAS)
C***KEYWORDS  BLAS, QUICK CHECK SERVICE ROUTINE
C***AUTHOR  Dongarra, J. J., (ANL)
C           Duff, I., (AERE)
C           Du Croz, J., (NAG)
C           Hammarling, S., (NAG)
C***DESCRIPTION
C
C  Quick check for CHER2K and CSYR2K.
C
C  Auxiliary routine for test program for Level 3 Blas.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CHER2K, CMAKE3, CMMCH, CSYR2K, LCE, LCERES, NUMXER
C***REVISION HISTORY  (YYMMDD)
C   890208  DATE WRITTEN
C   910619  Modified to meet SLATEC code and prologue standards.  (BKS)
C***END PROLOGUE  CCHK53
C     .. Parameters ..
      COMPLEX            ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0, 0.0 ), ONE = ( 1.0, 0.0 ) )
      REAL               RZERO, RONE
      PARAMETER          ( RZERO = 0.0, RONE = 1.0 )
C     .. Scalar Arguments ..
      LOGICAL            FATAL
      REAL               EPS, THRESH
      INTEGER            KPRINT, NALF, NBET, NIDIM, NMAX, NOUT
      CHARACTER*6        SNAME
C     .. Array Arguments ..
      COMPLEX            AA( NMAX*NMAX ), AB( 2*NMAX*NMAX ),
     $                   ALF( NALF ), AS( NMAX*NMAX ), BB( NMAX*NMAX ),
     $                   BET( NBET ), BS( NMAX*NMAX ), C( NMAX, NMAX ),
     $                   CC( NMAX*NMAX ), CS( NMAX*NMAX ), CT( NMAX ),
     $                   W( 2*NMAX )
      REAL               G( NMAX )
      INTEGER            IDIM( NIDIM )
C     .. Local Scalars ..
      COMPLEX            ALPHA, ALS, BETA, BETS
      REAL               ERR, ERRMAX, RBETA, RBETS
      INTEGER            I, IA, IB, ICT, ICU, IK, IN, J, JC, JJ, JJAB,
     $                   K, KS, LAA, LBB, LCC, LDA, LDAS, LDB, LDBS,
     $                   LDC, LDCS, LJ, MA, N, NA, NARGS, NC, NERR, NS
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
      EXTERNAL           CHER2K, CSYR2K, CMAKE3, CMMCH
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
C     .. Data statements ..
      DATA               ICHU/'UL'/, ICHT/'NC'/
C***FIRST EXECUTABLE STATEMENT  CCHK53
      CONJ = SNAME( 2: 3 ).EQ.'HE'
C
      NARGS = 12
      NC = 0
      RESET = .TRUE.
      ERRMAX = RZERO
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
C
         DO 120 IK = 1, NIDIM
            K = IDIM( IK )
C
            DO 110 ICT = 1, 2
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
     $            GO TO 110
               LAA = LDA*NA
C
C              Generate the matrix A.
C
               IF( TRAN )THEN
                  CALL CMAKE3( 'GE', ' ', ' ', MA, NA, AB, 2*NMAX, AA,
     $                        LDA, RESET, ZERO )
               ELSE
                  CALL CMAKE3( 'GE', ' ', ' ', MA, NA, AB, NMAX, AA,
     $                        LDA,  RESET, ZERO )
               END IF
C
C              Generate the matrix B.
C
               LDB = LDA
               LBB = LAA
               IF( TRAN )THEN
                  CALL CMAKE3( 'GE', ' ', ' ', MA, NA, AB( K + 1 ),
     $                        2*NMAX, BB, LDB, RESET, ZERO )
               ELSE
                  CALL CMAKE3( 'GE', ' ', ' ', MA, NA, AB( K*NMAX + 1 ),
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
                        IF( CONJ )THEN
                           RBETA = REAL( BETA )
                           BETA = CMPLX( RBETA, RZERO )
                        END IF
                        NULL = N.LE.0
                        IF( CONJ )
     $                     NULL = NULL.OR.( ( K.LE.0.OR.ALPHA.EQ.
     $                            ZERO ).AND.RBETA.EQ.RONE )
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
                        ALS = ALPHA
                        DO 10 I = 1, LAA
                           AS( I ) = AA( I )
   10                   CONTINUE
                        LDAS = LDA
                        DO 20 I = 1, LBB
                           BS( I ) = BB( I )
   20                   CONTINUE
                        LDBS = LDB
                        IF( CONJ )THEN
                           RBETS = RBETA
                        ELSE
                           BETS = BETA
                        END IF
                        DO 30 I = 1, LCC
                           CS( I ) = CC( I )
   30                   CONTINUE
                        LDCS = LDC
C
C                       Call the subroutine.
C
                        IF( CONJ )THEN
                           CALL CHER2K( UPLO, TRANS, N, K, ALPHA, AA,
     $                                  LDA, BB, LDB, RBETA, CC, LDC )
                        ELSE
                           CALL CSYR2K( UPLO, TRANS, N, K, ALPHA, AA,
     $                                  LDA, BB, LDB, BETA, CC, LDC )
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
                        ISAME( 5 ) = ALS.EQ.ALPHA
                        ISAME( 6 ) = LCE( AS, AA, LAA )
                        ISAME( 7 ) = LDAS.EQ.LDA
                        ISAME( 8 ) = LCE( BS, BB, LBB )
                        ISAME( 9 ) = LDBS.EQ.LDB
                        IF( CONJ )THEN
                           ISAME( 10 ) = RBETS.EQ.RBETA
                        ELSE
                           ISAME( 10 ) = BETS.EQ.BETA
                        END IF
                        IF( NULL )THEN
                           ISAME( 11 ) = LCE( CS, CC, LCC )
                        ELSE
                           ISAME( 11 ) = LCERES( 'HE', UPLO, N, N, CS,
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
                        IF( .NOT.NULL )THEN
C
C                          Check the result column by column.
C
                           IF( CONJ )THEN
                              TRANST = 'C'
                           ELSE
                              TRANST = 'T'
                           END IF
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
                                    W( I ) = ALPHA*AB( ( J - 1 )*2*
     $                                       NMAX + K + I )
                                    IF( CONJ )THEN
                                       W( K + I ) = CONJG( ALPHA )*
     $                                              AB( ( J - 1 )*2*
     $                                              NMAX + I )
                                    ELSE
                                       W( K + I ) = ALPHA*
     $                                              AB( ( J - 1 )*2*
     $                                              NMAX + I )
                                    END IF
   50                            CONTINUE
                                 FTL = .FALSE.
                                 CALL CMMCH( TRANST, 'N', LJ, 1, 2*K,
     $                                       ONE, AB( JJAB ), 2*NMAX, W,
     $                                       2*NMAX, BETA, C( JJ, J ),
     $                                       NMAX, CT, G, CC( JC ), LDC,
     $                                       EPS, ERR, FTL, NOUT,
     $                                       .TRUE., KPRINT )
                              ELSE
                                 DO 60 I = 1, K
                                    IF( CONJ )THEN
                                       W( I ) = ALPHA*CONJG( AB( ( K +
     $                                          I - 1 )*NMAX + J ) )
                                       W( K + I ) = CONJG( ALPHA*
     $                                              AB( ( I - 1 )*NMAX +
     $                                              J ) )
                                    ELSE
                                       W( I ) = ALPHA*AB( ( K + I - 1 )*
     $                                          NMAX + J )
                                       W( K + I ) = ALPHA*
     $                                              AB( ( I - 1 )*NMAX +
     $                                              J )
                                    END IF
   60                            CONTINUE
                                 FTL = .FALSE.
                                 CALL CMMCH( 'N', 'N', LJ, 1, 2*K, ONE,
     $                                       AB( JJ ), NMAX, W, 2*NMAX,
     $                                       BETA, C( JJ, J ), NMAX, CT,
     $                                       G, CC( JC ), LDC, EPS, ERR,
     $                                       FTL, NOUT, .TRUE.,KPRINT)
                              END IF
                              IF( UPPER )THEN
                                 JC = JC + LDC
                              ELSE
                                 JC = JC + LDC + 1
                                 IF( TRAN )
     $                              JJAB = JJAB + 2*NMAX
                              END IF
                              ERRMAX = MAX( ERRMAX, ERR )
                              IF (FTL) THEN
                              FATAL = .TRUE.
                              IF (KPRINT .GE. 3) THEN
                                WRITE( NOUT, FMT = 9996 )SNAME
                                IF( CONJ )THEN
                                   WRITE( NOUT, FMT = 9994 )NC, SNAME,
     $                                UPLO, TRANS, N, K, ALPHA,
     $                                LDA, LDB, RBETA, LDC
                                ELSE
                                   WRITE( NOUT, FMT = 9993 )NC, SNAME,
     $                                UPLO, TRANS, N, K, ALPHA,
     $                                LDA, LDB, BETA, LDC
                                END IF
                              ENDIF
                              ENDIF
   70                      CONTINUE
                        END IF
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
     $      '(', F4.1, ',', F4.1, '), A,', I3, ', B,', I3, ',', F4.1,
     $      ', C,', I3, ')           .' )
 9993 FORMAT( 1X, I6, ': ', A6, '(', 2( '''', A1, ''',' ), 2( I3, ',' ),
     $      '(', F4.1, ',', F4.1, '), A,', I3, ', B,', I3, ',(', F4.1,
     $      ',', F4.1, '), C,', I3, ')    .' )
 9992 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',
     $      '******' )
C
C     End of CCHK53.
C
      END
