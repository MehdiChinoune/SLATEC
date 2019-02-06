*DECK CCHK23
      SUBROUTINE CCHK23 (SNAME, EPS, THRESH, NOUT, KPRINT, FATAL, NIDIM,
     $   IDIM, NALF, ALF, NBET, BET, NMAX, A, AA, AS, B, BB, BS, C, CC,
     $   CS, CT, G)
C***BEGIN PROLOGUE  CCHK23
C***SUBSIDIARY
C***PURPOSE  Quick check for CHEMM and CSYMM.
C***LIBRARY   SLATEC (BLAS)
C***KEYWORDS  BLAS, QUICK CHECK SERVICE ROUTINE
C***AUTHOR  Dongarra, J. J., (ANL)
C           Duff, I., (AERE)
C           Du Croz, J., (NAG)
C           Hammarling, S., (NAG)
C***DESCRIPTION
C
C  Quick check for CHEMM and CSYMM.
C
C  Auxiliary routine for test program for Level 3 Blas.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CHEMM, CMAKE3, CMMCH, CSYMM, LCE, LCERES, NUMXER
C***REVISION HISTORY  (YYMMDD)
C   890208  DATE WRITTEN
C   910619  Modified to meet SLATEC code and prologue standards.  (BKS)
C***END PROLOGUE  CCHK23
C     .. Parameters ..
      COMPLEX            ZERO
      PARAMETER          ( ZERO = ( 0.0, 0.0 ) )
      REAL               RZERO
      PARAMETER          ( RZERO = 0.0 )
C     .. Scalar Arguments ..
      LOGICAL            FATAL
      REAL               EPS, THRESH
      INTEGER            KPRINT, NALF, NBET, NIDIM, NMAX, NOUT
      CHARACTER*6        SNAME
C     .. Array Arguments ..
      COMPLEX            A( NMAX, NMAX ), AA( NMAX*NMAX ), ALF( NALF ),
     $                   AS( NMAX*NMAX ), B( NMAX, NMAX ),
     $                   BB( NMAX*NMAX ), BET( NBET ), BS( NMAX*NMAX ),
     $                   C( NMAX, NMAX ), CC( NMAX*NMAX ),
     $                   CS( NMAX*NMAX ), CT( NMAX )
      REAL               G( NMAX )
      INTEGER            IDIM( NIDIM )
C     .. Local Scalars ..
      COMPLEX            ALPHA, ALS, BETA, BLS
      REAL               ERR, ERRMAX
      INTEGER            I, IA, IB, ICS, ICU, IM, IN, LAA, LBB, LCC,
     $                   LDA, LDAS, LDB, LDBS, LDC, LDCS, M, MS, N, NA,
     $                   NARGS, NC, NERR, NS
      LOGICAL            CONJ, FTL, LEFT, NULL, RESET
      CHARACTER*1        SIDE, SIDES, UPLO, UPLOS
      CHARACTER*2        ICHU, ICHS
C     .. Local Arrays ..
      LOGICAL            ISAME( 13 )
C     .. External Functions ..
      INTEGER            NUMXER
      LOGICAL            LCE, LCERES
      EXTERNAL           LCE, LCERES, NUMXER
C     .. External Subroutines ..
      EXTERNAL           CHEMM, CSYMM, CMAKE3, CMMCH
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
C     .. Data statements ..
      DATA               ICHS/'LR'/, ICHU/'UL'/
C***FIRST EXECUTABLE STATEMENT  CCHK23
      CONJ = SNAME( 2: 3 ).EQ.'HE'
C
      NARGS = 12
      NC = 0
      RESET = .TRUE.
      ERRMAX = RZERO
C
      DO 100 IM = 1, NIDIM
         M = IDIM( IM )
C
         DO 90 IN = 1, NIDIM
            N = IDIM( IN )
C           Set LDC to 1 more than minimum value if room.
            LDC = M
            IF( LDC.LT.NMAX )
     $         LDC = LDC + 1
C           Skip tests if not enough room.
            IF( LDC.GT.NMAX )
     $         GO TO 90
            LCC = LDC*N
            NULL = N.LE.0.OR.M.LE.0
C           Set LDB to 1 more than minimum value if room.
            LDB = M
            IF( LDB.LT.NMAX )
     $         LDB = LDB + 1
C           Skip tests if not enough room.
            IF( LDB.GT.NMAX )
     $         GO TO 90
            LBB = LDB*N
C
C           Generate the matrix B.
C
            CALL CMAKE3( 'GE', ' ', ' ', M, N, B, NMAX, BB, LDB, RESET,
     $                  ZERO )
C
            DO 80 ICS = 1, 2
               SIDE = ICHS( ICS: ICS )
               LEFT = SIDE.EQ.'L'
C
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
     $            GO TO 80
               LAA = LDA*NA
C
               DO 70 ICU = 1, 2
                  UPLO = ICHU( ICU: ICU )
C
C                 Generate the hermitian or symmetric matrix A.
C
                  CALL CMAKE3(SNAME( 2: 3 ), UPLO, ' ', NA, NA, A, NMAX,
     $                        AA, LDA, RESET, ZERO )
C
                  DO 60 IA = 1, NALF
                     ALPHA = ALF( IA )
C
                     DO 50 IB = 1, NBET
                        BETA = BET( IB )
C
C                       Generate the matrix C.
C
                        CALL CMAKE3( 'GE', ' ', ' ', M, N, C, NMAX, CC,
     $                              LDC, RESET, ZERO )
C
                        NC = NC + 1
C
C                       Save every datum before calling the
C                       subroutine.
C
                        SIDES = SIDE
                        UPLOS = UPLO
                        MS = M
                        NS = N
                        ALS = ALPHA
                        DO 10 I = 1, LAA
                           AS( I ) = AA( I )
   10                   CONTINUE
                        LDAS = LDA
                        DO 20 I = 1, LBB
                           BS( I ) = BB( I )
   20                   CONTINUE
                        LDBS = LDB
                        BLS = BETA
                        DO 30 I = 1, LCC
                           CS( I ) = CC( I )
   30                   CONTINUE
                        LDCS = LDC
C
C                       Call the subroutine.
C
                        IF( CONJ )THEN
                           CALL CHEMM( SIDE, UPLO, M, N, ALPHA, AA, LDA,
     $                                 BB, LDB, BETA, CC, LDC )
                        ELSE
                           CALL CSYMM( SIDE, UPLO, M, N, ALPHA, AA, LDA,
     $                                 BB, LDB, BETA, CC, LDC )
                        END IF
C
C                       Check if error-exit was taken incorrectly.
C
                        IF( NUMXER(NERR) .NE. 0 )THEN
                          IF (KPRINT .GE. 2) THEN
                           WRITE( NOUT, FMT = 9994 )
                           END IF
                          FATAL = .TRUE.
                        END IF
C
C                       See what data changed inside subroutines.
C
                        ISAME( 1 ) = SIDES.EQ.SIDE
                        ISAME( 2 ) = UPLOS.EQ.UPLO
                        ISAME( 3 ) = MS.EQ.M
                        ISAME( 4 ) = NS.EQ.N
                        ISAME( 5 ) = ALS.EQ.ALPHA
                        ISAME( 6 ) = LCE( AS, AA, LAA )
                        ISAME( 7 ) = LDAS.EQ.LDA
                        ISAME( 8 ) = LCE( BS, BB, LBB )
                        ISAME( 9 ) = LDBS.EQ.LDB
                        ISAME( 10 ) = BLS.EQ.BETA
                        IF( NULL )THEN
                           ISAME( 11 ) = LCE( CS, CC, LCC )
                        ELSE
                           ISAME( 11 ) = LCERES( 'GE', ' ', M, N, CS,
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
C                          Check the result.
C
                           IF( LEFT )THEN
                              CALL CMMCH( 'N', 'N', M, N, M, ALPHA, A,
     $                                    NMAX, B, NMAX, BETA, C, NMAX,
     $                                    CT, G, CC, LDC, EPS, ERR,
     $                                    FTL, NOUT, .TRUE.,
     $                                    KPRINT )
                           ELSE
                              CALL CMMCH( 'N', 'N', M, N, N, ALPHA, B,
     $                                    NMAX, A, NMAX, BETA, C, NMAX,
     $                                    CT, G, CC, LDC, EPS, ERR,
     $                                    FTL, NOUT, .TRUE., KPRINT )
                           END IF
                           ERRMAX = MAX( ERRMAX, ERR )
                        END IF
C
                        IF (FTL) THEN
                          FATAL = .TRUE.
                          IF (KPRINT .GE. 3) THEN
                            WRITE( NOUT, FMT = 9996 )SNAME
                            WRITE( NOUT, FMT = 9995 )NC, SNAME, SIDE,
     $                         UPLO, M, N, ALPHA, LDA, LDB, BETA,
     $                         LDC
                          ENDIF
                        ENDIF
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
 9995 FORMAT( 1X, I6, ': ', A6, '(', 2( '''', A1, ''',' ), 2( I3, ',' ),
     $      '(', F4.1, ',', F4.1, '), A,', I3, ', B,', I3, ',(', F4.1,
     $      ',', F4.1, '), C,', I3, ')    .' )
 9994 FORMAT( ' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',
     $      '******' )
C
C     End of CCHK23.
C
      END
