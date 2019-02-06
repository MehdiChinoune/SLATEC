*DECK DBLAT3
      SUBROUTINE DBLAT3 (NOUT, KPRINT, IPASS)
C***BEGIN PROLOGUE  DBLAT3
C***PURPOSE  Driver for testing Level 3 BLAS double precision
C            subroutines.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  A3B
C***TYPE      DOUBLE PRECISION (SBLAT3-S, DBLAT3-D, CBLAT3-C)
C***KEYWORDS  BLAS, QUICK CHECK DRIVER
C***AUTHOR  Dongarra, J. J., (ANL)
C           Duff, I., (AERE)
C           Du Croz, J., (NAG)
C           Hammarling, S., (NAG)
C***DESCRIPTION
C
C  Test program for the DOUBLE           Level 3 Blas.
C
C***REFERENCES  Dongarra, J., Du Croz, J., Duff, I., and Hammarling, S.
C                 A set of level 3 basic linear algebra subprograms.
C                 ACM TOMS, Vol. 16, No. 1, pp. 1-17, March 1990.
C***ROUTINES CALLED  DCHK13, DCHK23, DCHK33, DCHK43, DCHK53, DCHKE3,
C                    DMMCH, LDE, R1MACH, XERCLR
C***REVISION HISTORY  (YYMMDD)
C   890208  DATE WRITTEN
C   910619  Modified to meet SLATEC code and prologue standards. (BKS)
C   930315  Removed unused variables.  (WRB)
C   930618  Code modified to improve PASS/FAIL reporting.  (BKS, WRB)
C***END PROLOGUE  DBLAT3
C     .. Parameters ..
      INTEGER            NSUBS
      PARAMETER          ( NSUBS = 6)
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
      INTEGER            NMAX
      PARAMETER          ( NMAX = 65)
C     .. Scalar Arguments ..
      INTEGER            IPASS, KPRINT
C     .. Local Scalars ..
      DOUBLE PRECISION   EPS, ERR, THRESH
      INTEGER            I, ISNUM, J, N, NALF, NBET, NIDIM, NOUT
      PARAMETER          (NIDIM=6, NALF=3, NBET=3)
      LOGICAL            SAME, TSTERR, FTL, FTL1, FTL2
      CHARACTER*1        TRANSA, TRANSB
C     .. Local Arrays ..
      DOUBLE PRECISION   AB( NMAX, 2*NMAX ), AA( NMAX*NMAX ),
     $                   ALF( NALF ), AS( NMAX*NMAX ), BET( NBET ),
     $                   G( NMAX ), BB( NMAX*NMAX ),
     $                   BS( NMAX*NMAX ), C( NMAX, NMAX ),
     $                   CC( NMAX*NMAX ), CS( NMAX*NMAX ),  CT( NMAX ),
     $                   W( 2*NMAX )
      INTEGER            IDIM( NIDIM )
      LOGICAL            LTEST( NSUBS )
      CHARACTER*6        SNAMES( NSUBS )
C     .. External Functions ..
      REAL               R1MACH
      LOGICAL            LDE
      EXTERNAL           LDE, R1MACH
C     .. External Subroutines ..
      EXTERNAL           DCHK13, DCHK23, DCHK33, DCHK43, DCHK53,
     $                   DCHKE3, DMMCH
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
C     .. Data statements ..
      DATA               SNAMES/'DGEMM ', 'DSYMM ', 'DTRMM ', 'DTRSM ',
     $                   'DSYRK ', 'DSYR2K'/
      DATA               IDIM/0,1,2,3,5,9/
      DATA               ALF/0.0,1.0,0.7/
      DATA               BET/0.0,1.0,1.3/
C***FIRST EXECUTABLE STATEMENT  DBLAT3
C     Set the flag that indicates whether error exits are to be tested.
      TSTERR = .TRUE.
C     Set the threshold value of the test ratio
      THRESH = 16.0
C
C     Set IPASS to 1 assuming it will pass.
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
      DO 20 I = 1, NSUBS
         LTEST( I ) = .TRUE.
   20 CONTINUE
C
C     Set EPS (the machine precision).
C
      EPS = R1MACH (4)
C
C     Check the reliability of DMMCH using exact data.
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
C     CC holds the exact result. On exit from DMMCH CT holds
C     the result computed by DMMCH.
      TRANSA = 'N'
      TRANSB = 'N'
      FTL = .FALSE.
      CALL DMMCH( TRANSA, TRANSB, N, 1, N, ONE, AB, NMAX,
     $            AB( 1, NMAX + 1 ), NMAX, ZERO, C, NMAX, CT, G, CC,
     $            NMAX, EPS, ERR, FTL, NOUT, .TRUE., KPRINT )
      SAME = LDE( CC, CT, N )
      IF( .NOT.SAME.OR.ERR.NE.ZERO )THEN
        IPASS = 0
        IF (KPRINT .GE. 2) THEN
         WRITE( NOUT, FMT = 9989 )TRANSA, TRANSB, SAME, ERR
        ENDIF
      END IF
      TRANSB = 'T'
      FTL = .FALSE.
      CALL DMMCH( TRANSA, TRANSB, N, 1, N, ONE, AB, NMAX,
     $            AB( 1, NMAX + 1 ), NMAX, ZERO, C, NMAX, CT, G, CC,
     $            NMAX, EPS, ERR, FTL, NOUT, .TRUE., KPRINT )
      SAME = LDE( CC, CT, N )
      IF( .NOT.SAME.OR.ERR.NE.ZERO )THEN
        IPASS = 0
        IF (KPRINT .GE. 2) THEN
         WRITE( NOUT, FMT = 9989 )TRANSA, TRANSB, SAME, ERR
        END IF
      ENDIF
      DO 120 J = 1, N
         AB( J, NMAX + 1 ) = N - J + 1
         AB( 1, NMAX + J ) = N - J + 1
  120 CONTINUE
      DO 130 J = 1, N
         CC( N - J + 1 ) = J*( ( J + 1 )*J )/2 -
     $                     ( ( J + 1 )*J*( J - 1 ) )/3
  130 CONTINUE
      TRANSA = 'T'
      TRANSB = 'N'
      FTL = .FALSE.
      CALL DMMCH( TRANSA, TRANSB, N, 1, N, ONE, AB, NMAX,
     $            AB( 1, NMAX + 1 ), NMAX, ZERO, C, NMAX, CT, G, CC,
     $            NMAX, EPS, ERR, FTL, NOUT, .TRUE., KPRINT )
      SAME = LDE( CC, CT, N )
      IF( .NOT.SAME.OR.ERR.NE.ZERO )THEN
        IPASS = 0
        IF (KPRINT .GE. 2) THEN
         WRITE( NOUT, FMT = 9989 )TRANSA, TRANSB, SAME, ERR
        ENDIF
      END IF
      TRANSB = 'T'
      FTL = .FALSE.
      CALL DMMCH( TRANSA, TRANSB, N, 1, N, ONE, AB, NMAX,
     $            AB( 1, NMAX + 1 ), NMAX, ZERO, C, NMAX, CT, G, CC,
     $            NMAX, EPS, ERR, FTL, NOUT, .TRUE., KPRINT )
      SAME = LDE( CC, CT, N )
      IF( .NOT.SAME.OR.ERR.NE.ZERO )THEN
        IPASS = 0
        IF (KPRINT .GE. 2) THEN
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
               CALL DCHKE3(ISNUM, SNAMES( ISNUM ), NOUT, KPRINT, FTL1)
            END IF
C           Test computations.
            FTL2 = .FALSE.
            CALL XERCLR
            GO TO ( 140, 150, 160, 160, 170, 180 )ISNUM
C           Test DGEMM, 01.
  140       CALL DCHK13( SNAMES( ISNUM ), EPS, THRESH, NOUT, KPRINT,
     $                  FTL2, NIDIM, IDIM, NALF, ALF, NBET, BET,
     $                  NMAX, AB, AA, AS, AB( 1, NMAX + 1 ), BB, BS, C,
     $                  CC, CS, CT, G )
            GO TO 190
C           Test DSYMM, 02.
  150       CALL DCHK23( SNAMES( ISNUM ), EPS, THRESH, NOUT, KPRINT,
     $                  FTL2, NIDIM, IDIM, NALF, ALF, NBET, BET,
     $                  NMAX, AB, AA, AS, AB( 1, NMAX + 1 ), BB, BS, C,
     $                  CC, CS, CT, G )
            GO TO 190
C           Test DTRMM, 03, DTRSM, 04.
  160       CALL DCHK33( SNAMES( ISNUM ), EPS, THRESH, NOUT, KPRINT,
     $                  FTL2, NIDIM, IDIM, NALF, ALF, NMAX, AB,
     $                  AA, AS, AB( 1, NMAX + 1 ), BB, BS, CT, G, C )
            GO TO 190
C           Test DSYRK, 05.
  170       CALL DCHK43( SNAMES( ISNUM ), EPS, THRESH, NOUT, KPRINT,
     $                  FTL2, NIDIM, IDIM, NALF, ALF, NBET, BET,
     $                  NMAX, AB, AA, AS, AB( 1, NMAX + 1 ), BB, BS, C,
     $                  CC, CS, CT, G )
            GO TO 190
C           Test DSYR2K, 06.
  180       CALL DCHK53( SNAMES( ISNUM ), EPS, THRESH, NOUT, KPRINT,
     $                  FTL2, NIDIM, IDIM, NALF, ALF, NBET, BET,
     $                  NMAX, AB, AA, AS, BB, BS, C, CC, CS, CT, G, W )
            GO TO 190
C
  190       IF( FTL1 .OR. FTL2) THEN
              IPASS = 0
            ENDIF
         END IF
  200 CONTINUE
      RETURN
C
 9999 FORMAT( ' ROUTINES PASS COMPUTATIONAL TESTS IF TEST RATIO IS LES',
     $      'S THAN', F8.2 )
 9995 FORMAT( ' TESTS OF THE DOUBLE PRECISION LEVEL 3 BLAS', //' THE F',
     $      'OLLOWING PARAMETER VALUES WILL BE USED:' )
 9994 FORMAT( '   FOR N              ', 9I6 )
 9993 FORMAT( '   FOR ALPHA          ', 7F6.1 )
 9992 FORMAT( '   FOR BETA           ', 7F6.1 )
 9989 FORMAT( ' ERROR IN DMMCH -  IN-LINE DOT PRODUCTS ARE BEING EVALU',
     $      'ATED WRONGLY.', /' DMMCH WAS CALLED WITH TRANSA = ', A1,
     $      ' AND TRANSB = ', A1, /' AND RETURNED SAME = ', L1, ' AND ',
     $      'ERR = ', F12.3, '.', /' THIS MAY BE DUE TO FAULTS IN THE ',
     $      'ARITHMETIC OR THE COMPILER.')
 9987 FORMAT( 1X, A6, ' WAS NOT TESTED' )
 9984 FORMAT( ' ERROR-EXITS WILL NOT BE TESTED' )
C
C     End of DBLAT3.
C
      END
