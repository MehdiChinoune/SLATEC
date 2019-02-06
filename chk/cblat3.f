*DECK CBLAT3
      SUBROUTINE CBLAT3 (NOUT, KPRINT, IPASS)
C***BEGIN PROLOGUE  CBLAT3
C***PURPOSE  Driver for testing Level 3 BLAS complex subroutines.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  A4
C***TYPE      COMPLEX (SBLAT3-S, DBLAT3-D, CBLAT3-C)
C***KEYWORDS  BLAS, QUICK CHECK DRIVER
C***AUTHOR  Dongarra, J. J., (ANL)
C           Duff, I., (AERE)
C           Du Croz, J., (NAG)
C           Hammarling, S., (NAG)
C***DESCRIPTION
C
C  Test program for the COMPLEX              Level 3 Blas.
C
C***REFERENCES  Dongarra, J., Du Croz, J., Duff, I., and Hammarling, S.
C                 A set of level 3 basic linear algebra subprograms.
C                 ACM TOMS, Vol. 16, No. 1, pp. 1-17, March 1990.
C***ROUTINES CALLED  CCHK13, CCHK23, CCHK33, CCHK43, CCHK53, CCHKE3,
C                    CMMCH, LCE, R1MACH, XERCLR
C***REVISION HISTORY  (YYMMDD)
C   890208  DATE WRITTEN
C   910619  Modified to meet SLATEC code and prologue standards.  (BKS)
C   930315  Removed unused variables.  (WRB)
C   930618  Code modified to improve PASS/FAIL reporting.  (BKS, WRB)
C***END PROLOGUE  CBLAT3
C     .. Parameters ..
      INTEGER            NSUBS
      PARAMETER          ( NSUBS = 9)
      COMPLEX            ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0, 0.0 ), ONE = ( 1.0, 0.0 ) )
      REAL               RZERO
      PARAMETER          ( RZERO = 0.0 )
      INTEGER            NMAX
      PARAMETER          ( NMAX = 65)
C     .. Scalar Arguments ..
      INTEGER            IPASS, KPRINT
C     .. Local Scalars ..
      REAL               EPS, ERR, THRESH
      INTEGER            I, ISNUM, J, N, NALF, NBET, NIDIM, NOUT
      PARAMETER          (NIDIM=6, NALF=3, NBET=3)
      LOGICAL            SAME, TSTERR, FTL, FTL1, FTL2
      CHARACTER*1        TRANSA, TRANSB
C     .. Local Arrays ..
      COMPLEX            AA( NMAX*NMAX ), AB( NMAX, 2*NMAX ),
     $                   ALF( NALF ), AS( NMAX*NMAX ),
     $                   BB( NMAX*NMAX ), BET( NBET ),
     $                   BS( NMAX*NMAX ), C( NMAX, NMAX ),
     $                   CC( NMAX*NMAX ), CS( NMAX*NMAX ), CT( NMAX ),
     $                   W( 2*NMAX )
      REAL               G( NMAX )
      INTEGER            IDIM( NIDIM )
      LOGICAL            LTEST( NSUBS )
      CHARACTER*6        SNAMES( NSUBS )
C     .. External Functions ..
      REAL               R1MACH
      LOGICAL            LCE
      EXTERNAL           LCE, R1MACH
C     .. External Subroutines ..
      EXTERNAL           CCHK13, CCHK23, CCHK33, CCHK43, CCHK53,
     $                   CCHKE3, CMMCH
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
C     .. Data statements ..
      DATA               SNAMES/'CGEMM ', 'CHEMM ', 'CSYMM ', 'CTRMM ',
     $                   'CTRSM ', 'CHERK ', 'CSYRK ', 'CHER2K',
     $                   'CSYR2K'/
      DATA               IDIM/0,1,2,3,5,9/
      DATA               ALF/(0.0,0.0),(1.0,0.0),(0.7,-0.9)/
      DATA               BET/(0.0,0.0),(1.0,0.0),(1.3,-1.1)/
C***FIRST EXECUTABLE STATEMENT  CBLAT3
C
C     Set the flag that indicates whether error exits are to be tested.
      TSTERR = .TRUE.
C     Set the threshold value of the test ratio
      THRESH = 16.0
C
C     Set IPASS = 1 assuming all tests will pass.
C
      IPASS = 1
C
C     Report values of parameters.
C
      IF (KPRINT .GE. 3) THEN
        WRITE( NOUT, FMT = 9995 )
        WRITE( NOUT, FMT = 9994 )( IDIM( I ), I = 1, NIDIM )
        WRITE( NOUT, FMT = 9993 )( ALF( I ), I = 1, NALF )
        WRITE( NOUT, FMT = 9992 )( BET( I ), I = 1, NBET )
        IF( .NOT.TSTERR )THEN
           WRITE( NOUT, FMT = 9984 )
        END IF
        WRITE( NOUT, FMT = 9999 )THRESH
      ENDIF
C
C     Set names of subroutines and flags which indicate
C     whether they are to be tested.
C
      DO 40 I = 1, NSUBS
         LTEST( I ) = .TRUE.
   40 CONTINUE
C
C     Set EPS (the machine precision).
C
      EPS = R1MACH (4)
C
C     Check the reliability of CMMCH using exact data.
C
      N = MIN( 32, NMAX )
      DO 100 J = 1, N
         DO 90 I = 1, N
            AB( I, J ) = MAX( I - J + 1, 0 )
   90    CONTINUE
         AB( J, NMAX + 1 ) = J
         AB( 1, NMAX + J ) = J
         C( J, 1 ) = ZERO
  100 CONTINUE
      DO 110 J = 1, N
         CC( J ) = J*( ( J + 1 )*J )/2 - ( ( J + 1 )*J*( J - 1 ) )/3
  110 CONTINUE
C     CC holds the exact result. On exit from CMMCH CT holds
C     the result computed by CMMCH.
      TRANSA = 'N'
      TRANSB = 'N'
      FTL = .FALSE.
      CALL CMMCH( TRANSA, TRANSB, N, 1, N, ONE, AB, NMAX,
     $            AB( 1, NMAX + 1 ), NMAX, ZERO, C, NMAX, CT, G, CC,
     $            NMAX, EPS, ERR, FTL, NOUT, .TRUE., KPRINT )
      SAME = LCE( CC, CT, N )
      IF( .NOT.SAME.OR.ERR.NE.RZERO )THEN
        IPASS = 0
        IF (KPRINT .GE. 2) THEN
         WRITE( NOUT, FMT = 9989 )TRANSA, TRANSB, SAME, ERR
        ENDIF
      END IF
      TRANSB = 'C'
      FTL = .FALSE.
      CALL CMMCH( TRANSA, TRANSB, N, 1, N, ONE, AB, NMAX,
     $            AB( 1, NMAX + 1 ), NMAX, ZERO, C, NMAX, CT, G, CC,
     $            NMAX, EPS, ERR, FTL, NOUT, .TRUE., KPRINT )
      SAME = LCE( CC, CT, N )
      IF( .NOT.SAME.OR.ERR.NE.RZERO )THEN
        IPASS = 0
        IF (KPRINT .GE. 2) THEN
         WRITE( NOUT, FMT = 9989 )TRANSA, TRANSB, SAME, ERR
        ENDIF
      END IF
      DO 120 J = 1, N
         AB( J, NMAX + 1 ) = N - J + 1
         AB( 1, NMAX + J ) = N - J + 1
  120 CONTINUE
      DO 130 J = 1, N
         CC( N - J + 1 ) = J*( ( J + 1 )*J )/2 -
     $                     ( ( J + 1 )*J*( J - 1 ) )/3
  130 CONTINUE
      TRANSA = 'C'
      TRANSB = 'N'
      FTL = .FALSE.
      CALL CMMCH( TRANSA, TRANSB, N, 1, N, ONE, AB, NMAX,
     $            AB( 1, NMAX + 1 ), NMAX, ZERO, C, NMAX, CT, G, CC,
     $            NMAX, EPS, ERR, FTL, NOUT, .TRUE., KPRINT )
      SAME = LCE( CC, CT, N )
      IF( .NOT.SAME.OR.ERR.NE.RZERO )THEN
        IPASS = 0
        IF (KPRINT .GE. 2 ) THEN
         WRITE( NOUT, FMT = 9989 )TRANSA, TRANSB, SAME, ERR
        ENDIF
      END IF
      TRANSB = 'C'
      FTL = .FALSE.
      CALL CMMCH( TRANSA, TRANSB, N, 1, N, ONE, AB, NMAX,
     $            AB( 1, NMAX + 1 ), NMAX, ZERO, C, NMAX, CT, G, CC,
     $            NMAX, EPS, ERR, FTL, NOUT, .TRUE., KPRINT )
      SAME = LCE( CC, CT, N )
      IF( .NOT.SAME.OR.ERR.NE.RZERO )THEN
        IPASS = 0
        IF (KPRINT .GE. 2 ) THEN
         WRITE( NOUT, FMT = 9989 )TRANSA, TRANSB, SAME, ERR
        ENDIF
      END IF
C
C     Test each subroutine in turn.
C
      DO 200 ISNUM = 1, NSUBS
         IF( .NOT.LTEST( ISNUM ) )THEN
C           Subprogram is not to be tested.
            WRITE( NOUT, FMT = 9987 )SNAMES( ISNUM )
         ELSE
C           Test error exits.
            FTL1 = .FALSE.
            IF( TSTERR )THEN
               CALL CCHKE3(ISNUM, SNAMES( ISNUM ), NOUT, KPRINT, FTL1)
            END IF
C           Test computations.
            FTL2 = .FALSE.
            CALL XERCLR
            GO TO ( 140, 150, 150, 160, 160, 170, 170,
     $              180, 180 )ISNUM
C           Test CGEMM, 01.
  140       CALL CCHK13( SNAMES( ISNUM ), EPS, THRESH, NOUT, KPRINT,
     $                  FTL2, NIDIM, IDIM, NALF, ALF, NBET, BET,
     $                  NMAX, AB, AA, AS, AB( 1, NMAX + 1 ), BB, BS, C,
     $                  CC, CS, CT, G )
            GO TO 190
C           Test CHEMM, 02, CSYMM, 03.
  150       CALL CCHK23( SNAMES( ISNUM ), EPS, THRESH, NOUT, KPRINT,
     $                  FTL2, NIDIM, IDIM, NALF, ALF, NBET, BET,
     $                  NMAX, AB, AA, AS, AB( 1, NMAX + 1 ), BB, BS, C,
     $                  CC, CS, CT, G )
            GO TO 190
C           Test CTRMM, 04, CTRSM, 05.
  160       CALL CCHK33( SNAMES( ISNUM ), EPS, THRESH, NOUT, KPRINT,
     $                  FTL2, NIDIM, IDIM, NALF, ALF, NMAX, AB,
     $                  AA, AS, AB( 1, NMAX + 1 ), BB, BS, CT, G, C )
            GO TO 190
C           Test CHERK, 06, CSYRK, 07.
  170       CALL CCHK43( SNAMES( ISNUM ), EPS, THRESH, NOUT, KPRINT,
     $                  FTL2, NIDIM, IDIM, NALF, ALF, NBET, BET,
     $                  NMAX, AB, AA, AS, AB( 1, NMAX + 1 ), BB, BS, C,
     $                  CC, CS, CT, G )
            GO TO 190
C           Test CHER2K, 08, CSYR2K, 09.
  180       CALL CCHK53( SNAMES( ISNUM ), EPS, THRESH, NOUT, KPRINT,
     $                  FTL2, NIDIM, IDIM, NALF, ALF, NBET, BET,
     $                  NMAX, AB, AA, AS, BB, BS, C, CC, CS, CT, G, W )
            GO TO 190
C
  190       IF( FTL1 .OR. FTL2 ) THEN
              IPASS = 0
            ENDIF
         END IF
  200 CONTINUE
      RETURN
C
 9999 FORMAT( ' ROUTINES PASS COMPUTATIONAL TESTS IF TEST RATIO IS LES',
     $      'S THAN', F8.2 )
 9995 FORMAT( ' TESTS OF THE COMPLEX          LEVEL 3 BLAS', //' THE F',
     $      'OLLOWING PARAMETER VALUES WILL BE USED:' )
 9994 FORMAT( '   FOR N              ', 9I6 )
 9993 FORMAT( '   FOR ALPHA          ',
     $      7( '(', F4.1, ',', F4.1, ')  ', : ) )
 9992 FORMAT( '   FOR BETA           ',
     $      7( '(', F4.1, ',', F4.1, ')  ', : ) )
 9989 FORMAT( ' ERROR IN CMMCH -  IN-LINE DOT PRODUCTS ARE BEING EVALU',
     $      'ATED WRONGLY.', /' CMMCH WAS CALLED WITH TRANSA = ', A1,
     $      ' AND TRANSB = ', A1, /' AND RETURNED SAME = ', L1, ' AND ',
     $      'ERR = ', F12.3, '.', /' THIS MAY BE DUE TO FAULTS IN THE ',
     $      'ARITHMETIC OR THE COMPILER.')
 9987 FORMAT( 1X, A6, ' WAS NOT TESTED' )
 9984 FORMAT( ' ERROR-EXITS WILL NOT BE TESTED' )
C
C     End of CBLAT3.
C
      END
