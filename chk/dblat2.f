*DECK DBLAT2
      SUBROUTINE DBLAT2 (NOUT, KPRINT, IPASS)
C***BEGIN PROLOGUE  DBLAT2
C***PURPOSE  Driver for testing Level 2 BLAS double precision
C            subroutines.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  A3B
C***TYPE      DOUBLE PRECISION (SBLAT2-S, DBLAT2-D, CBLAT2-C)
C***KEYWORDS  BLAS, QUICK CHECK DRIVER
C***AUTHOR  Du Croz, J. J., (NAG)
C           Hanson, R. J., (SNLA)
C***DESCRIPTION
C
C  Test program for the DOUBLE           Level 2 Blas.
C
C***REFERENCES  Dongarra, J. J., Du Croz, J. J., Hammarling, S. and
C                 Hanson, R. J.  An  extended  set of Fortran Basic
C                 Linear Algebra Subprograms. ACM TOMS, Vol. 14, No. 1,
C                 pp. 1-17, March 1988.
C***ROUTINES CALLED  DCHK12, DCHK22, DCHK32, DCHK42, DCHK52, DCHK62,
C                    DCHKE2, DMVCH, LDE, R1MACH, XERCLR
C***REVISION HISTORY  (YYMMDD)
C   870810  DATE WRITTEN
C   910619  Modified to meet SLATEC code and prologue standards. (BKS)
C   930315  Removed unused variables.  (WRB)
C   930618  Code modified to improve PASS/FAIL reporting.  (BKS, WRB)
C***END PROLOGUE  DBLAT2
C     .. Parameters ..
      INTEGER            NSUBS
      PARAMETER          ( NSUBS = 16)
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
      INTEGER            NMAX, INCMAX
      PARAMETER          ( NMAX = 65, INCMAX = 2 )
C     .. Scalar Arguments ..
      INTEGER            IPASS, KPRINT
C     .. Local Scalars ..
      DOUBLE PRECISION   EPS, ERR, THRESH
      INTEGER            I, ISNUM, J, N, NALF, NBET, NIDIM, NINC, NKB,
     $                   NOUT
      PARAMETER          (NIDIM=6, NKB=4, NINC=4, NALF=3, NBET=3)
      LOGICAL            SAME, TSTERR, FTL, FTL1, FTL2
      CHARACTER*1        TRANS
C     .. Local Arrays ..
      DOUBLE PRECISION   A( NMAX, NMAX ), AA( NMAX*NMAX ),
     $                   ALF( NALF ), AS( NMAX*NMAX ), BET( NBET ),
     $                   G( NMAX ), X( NMAX ), XS( NMAX*INCMAX ),
     $                   XX( NMAX*INCMAX ), Y( NMAX ),
     $                   YS( NMAX*INCMAX ), YT( NMAX ),
     $                   YY( NMAX*INCMAX ), Z( 2*NMAX )
      INTEGER            IDIM( NIDIM ), INC( NINC ), KB( NKB )
      LOGICAL            LTEST( NSUBS )
      CHARACTER*6        SNAMES( NSUBS )
C     .. External Functions ..
      REAL               R1MACH
      LOGICAL            LDE
      EXTERNAL           LDE, R1MACH
C     .. External Subroutines ..
      EXTERNAL           DCHK12, DCHK22, DCHK32, DCHK42, DCHK52, DCHK62,
     $                   DCHKE2, DMVCH
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
C     .. Data statements ..
      DATA               SNAMES/'DGEMV ', 'DGBMV ', 'DSYMV ', 'DSBMV ',
     $                   'DSPMV ', 'DTRMV ', 'DTBMV ', 'DTPMV ',
     $                   'DTRSV ', 'DTBSV ', 'DTPSV ', 'DGER  ',
     $                   'DSYR  ', 'DSPR  ', 'DSYR2 ', 'DSPR2 '/
      DATA               IDIM/0,1,2,3,5,9/
      DATA               KB/0,1,2,4/
      DATA               INC/1,2,-1,-2/
      DATA               ALF/0.0,1.0,0.7/
      DATA               BET/0.0,1.0,0.9/
C***FIRST EXECUTABLE STATEMENT  DBLAT2
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
      WRITE( NOUT, FMT = 9993 )
      WRITE( NOUT, FMT = 9992 )( IDIM( I ), I = 1, NIDIM )
      WRITE( NOUT, FMT = 9991 )( KB( I ), I = 1, NKB )
      WRITE( NOUT, FMT = 9990 )( INC( I ), I = 1, NINC )
      WRITE( NOUT, FMT = 9989 )( ALF( I ), I = 1, NALF )
      WRITE( NOUT, FMT = 9988 )( BET( I ), I = 1, NBET )
      IF( .NOT.TSTERR )THEN
         WRITE( NOUT, FMT = 9980 )
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
C     Check the reliability of DMVCH using exact data.
C
      N = MIN( 32, NMAX )
      DO 120 J = 1, N
         DO 110 I = 1, N
            A( I, J ) = MAX( I - J + 1, 0 )
  110    CONTINUE
         X( J ) = J
         Y( J ) = ZERO
  120 CONTINUE
      DO 130 J = 1, N
         YY( J ) = J*( ( J + 1 )*J )/2 - ( ( J + 1 )*J*( J - 1 ) )/3
  130 CONTINUE
C     YY holds the exact result. On exit from DMVCH YT holds
C     the result computed by DMVCH.
      TRANS = 'N'
      FTL = .FALSE.
      CALL DMVCH( TRANS, N, N, ONE, A, NMAX, X, 1, ZERO, Y, 1, YT, G,
     $            YY, EPS, ERR, FTL, NOUT, .TRUE., KPRINT )
      SAME = LDE( YY, YT, N )
      IF( .NOT.SAME.OR.ERR.NE.ZERO )THEN
        IPASS = 0
        IF (KPRINT .GE. 2) THEN
         WRITE( NOUT, FMT = 9985 )TRANS, SAME, ERR
        END IF
      ENDIF
      TRANS = 'T'
      FTL = .FALSE.
      CALL DMVCH( TRANS, N, N, ONE, A, NMAX, X, -1, ZERO, Y, -1, YT, G,
     $            YY, EPS, ERR, FTL, NOUT, .TRUE., KPRINT )
      SAME = LDE( YY, YT, N )
      IF( .NOT.SAME.OR.ERR.NE.ZERO )THEN
        IPASS = 0
        IF (KPRINT .GE. 2) THEN
         WRITE( NOUT, FMT = 9985 )TRANS, SAME, ERR
        END IF
      ENDIF
C
C     Test each subroutine in turn.
C
      DO 210 ISNUM = 1, NSUBS
         IF( .NOT.LTEST( ISNUM ) )THEN
C           Subprogram is not to be tested.
            WRITE( NOUT, FMT = 9983 )SNAMES( ISNUM )
         ELSE
C           Test error exits.
            FTL1 = .FALSE.
            IF( TSTERR )THEN
               CALL DCHKE2(ISNUM, SNAMES( ISNUM ), NOUT, KPRINT, FTL1)
            END IF
C           Test computations.
            CALL XERCLR
            FTL2 = .FALSE.
            GO TO ( 140, 140, 150, 150, 150, 160, 160,
     $              160, 160, 160, 160, 170, 180, 180,
     $              190, 190 )ISNUM
C           Test DGEMV, 01, and DGBMV, 02.
  140       CALL DCHK12( SNAMES( ISNUM ), EPS, THRESH, NOUT, KPRINT,
     $                  FTL2, NIDIM, IDIM, NKB, KB, NALF, ALF,
     $                  NBET, BET, NINC, INC, NMAX, INCMAX, A, AA, AS,
     $                  X, XX, XS, Y, YY, YS, YT, G )
            GO TO 200
C           Test DSYMV, 03, DSBMV, 04, and DSPMV, 05.
  150       CALL DCHK22( SNAMES( ISNUM ), EPS, THRESH, NOUT, KPRINT,
     $                  FTL2, NIDIM, IDIM, NKB, KB, NALF, ALF,
     $                  NBET, BET, NINC, INC, NMAX, INCMAX, A, AA, AS,
     $                  X, XX, XS, Y, YY, YS, YT, G )
            GO TO 200
C           Test DTRMV, 06, DTBMV, 07, DTPMV, 08,
C           DTRSV, 09, DTBSV, 10, and DTPSV, 11.
  160       CALL DCHK32( SNAMES( ISNUM ), EPS, THRESH, NOUT, KPRINT,
     $                  FTL2, NIDIM, IDIM, NKB, KB, NINC, INC,
     $                  NMAX, INCMAX, A, AA, AS, Y, YY, YS, YT, G, Z )
            GO TO 200
C           Test DGER, 12.
  170       CALL DCHK42( SNAMES( ISNUM ), EPS, THRESH, NOUT, KPRINT,
     $                  FTL2, NIDIM, IDIM, NALF, ALF, NINC, INC,
     $                  NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS,
     $                  YT, G, Z )
            GO TO 200
C           Test DSYR, 13, and DSPR, 14.
  180       CALL DCHK52( SNAMES( ISNUM ), EPS, THRESH, NOUT, KPRINT,
     $                  FTL2, NIDIM, IDIM, NALF, ALF, NINC, INC,
     $                  NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS,
     $                  YT, G, Z )
            GO TO 200
C           Test DSYR2, 15, and DSPR2, 16.
  190       CALL DCHK62( SNAMES( ISNUM ), EPS, THRESH, NOUT, KPRINT,
     $                  FTL2, NIDIM, IDIM, NALF, ALF, NINC, INC,
     $                  NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS,
     $                  YT, G, Z )
C
  200      IF (FTL1 .OR. FTL2) THEN
             IPASS = 0
            ENDIF
         END IF
  210 CONTINUE
      RETURN
C
 9999 FORMAT( ' ROUTINES PASS COMPUTATIONAL TESTS IF TEST RATIO IS LES',
     $      'S THAN', F8.2 )
 9993 FORMAT( ' TESTS OF THE DOUBLE PRECISION LEVEL 2 BLAS', //' THE F',
     $      'OLLOWING PARAMETER VALUES WILL BE USED:' )
 9992 FORMAT( '   FOR N              ', 9I6 )
 9991 FORMAT( '   FOR K              ', 7I6 )
 9990 FORMAT( '   FOR INCX AND INCY  ', 7I6 )
 9989 FORMAT( '   FOR ALPHA          ', 7F6.1 )
 9988 FORMAT( '   FOR BETA           ', 7F6.1 )
 9985 FORMAT( ' ERROR IN DMVCH -  IN-LINE DOT PRODUCTS ARE BEING EVALU',
     $      'ATED WRONGLY.', /' DMVCH WAS CALLED WITH TRANS = ', A1,
     $      ' AND RETURNED SAME = ', L1, ' AND ERR = ', F12.3, '.', /
     $      ' THIS MAY BE DUE TO FAULTS IN THE ARITHMETIC OR THE',
     $      ' COMPILER.')
 9983 FORMAT( 1X, A6, ' WAS NOT TESTED' )
 9980 FORMAT( ' ERROR-EXITS WILL NOT BE TESTED' )
C
C     End of DBLAT2.
C
      END
