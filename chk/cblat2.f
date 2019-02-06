*DECK CBLAT2
      SUBROUTINE CBLAT2 (NOUT, KPRINT, IPASS)
C***BEGIN PROLOGUE  CBLAT2
C***PURPOSE  Driver for testing Level 2 BLAS complex subroutines.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  A4
C***TYPE      COMPLEX (SBLAT2-S, DBLAT2-D, CBLAT2-C)
C***KEYWORDS  BLAS, QUICK CHECK DRIVER
C***AUTHOR  Du Croz, J. J., (NAG)
C           Hanson, R. J., (SNLA)
C***DESCRIPTION
C
C  Test program for the COMPLEX              Level 2 Blas.
C
C***REFERENCES  Dongarra, J. J., Du Croz, J. J., Hammarling, S. and
C                 Hanson, R. J.  An  extended  set of Fortran Basic
C                 Linear Algebra Subprograms. ACM TOMS, Vol. 14, No. 1,
C                 pp. 1-17, March 1988.
C***ROUTINES CALLED  CCHK12, CCHK22, CCHK32, CCHK42, CCHK52, CCHK62,
C                    CCHKE2, CMVCH, LCE, R1MACH, XERCLR
C***REVISION HISTORY  (YYMMDD)
C   870810  DATE WRITTEN
C   910619  Modified to meet SLATEC code and prologue standards.  (BKS)
C   930315  Removed unused variables.  (WRB)
C   930618  Code modified to improve PASS/FAIL reporting.  (BKS, WRB)
C***END PROLOGUE  CBLAT2
C     .. Parameters ..
      INTEGER            NSUBS
      PARAMETER          ( NSUBS = 17)
      COMPLEX            ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0, 0.0 ), ONE = ( 1.0, 0.0 ) )
      REAL               RZERO
      PARAMETER          ( RZERO = 0.0 )
      INTEGER            NMAX, INCMAX
      PARAMETER          ( NMAX = 65, INCMAX = 2 )
C     .. Scalar Arguments ..
      INTEGER            IPASS, KPRINT
C     .. Local Scalars ..
      REAL               EPS, ERR, THRESH
      INTEGER            I, ISNUM, J, N, NALF, NBET, NIDIM, NINC,
     $                   NKB, NOUT
      PARAMETER          (NIDIM=6, NKB=4, NINC=4, NALF=3, NBET=3)
      LOGICAL            SAME, TSTERR, FTL, FTL1, FTL2
      CHARACTER*1        TRANS
C     .. Local Arrays ..
      COMPLEX            A( NMAX, NMAX ), AA( NMAX*NMAX ),
     $                   ALF( NALF ), AS( NMAX*NMAX ), BET( NBET ),
     $                   X( NMAX ), XS( NMAX*INCMAX ),
     $                   XX( NMAX*INCMAX ), Y( NMAX ),
     $                   YS( NMAX*INCMAX ), YT( NMAX ),
     $                   YY( NMAX*INCMAX ), Z( 2*NMAX )
      REAL               G( NMAX )
      INTEGER            IDIM( NIDIM ), INC( NINC ), KB( NKB )
      LOGICAL            LTEST( NSUBS )
      CHARACTER*6        SNAMES( NSUBS )
C     .. External Functions ..
      REAL               R1MACH
      LOGICAL            LCE
      EXTERNAL           LCE, R1MACH
C     .. External Subroutines ..
      EXTERNAL           CCHK12, CCHK22, CCHK32, CCHK42, CCHK52, CCHK62,
     $                   CCHKE2, CMVCH
C     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
C     .. Data statements ..
      DATA               SNAMES/'CGEMV ', 'CGBMV ', 'CHEMV ', 'CHBMV ',
     $                   'CHPMV ', 'CTRMV ', 'CTBMV ', 'CTPMV ',
     $                   'CTRSV ', 'CTBSV ', 'CTPSV ', 'CGERC ',
     $                   'CGERU ', 'CHER  ', 'CHPR  ', 'CHER2 ',
     $                   'CHPR2 '/
      DATA               IDIM/0,1,2,3,5,9/
      DATA               KB/0,1,2,4/
      DATA               INC/1,2,-1,-2/
      DATA               ALF/(0.0,0.0),(1.0,0.0),(0.7,-0.9)/
      DATA               BET/(0.0,0.0),(1.0,0.0),(1.3,-1.1)/
C***FIRST EXECUTABLE STATEMENT  CBLAT2
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
C     Check the reliability of CMVCH using exact data.
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
C     YY holds the exact result. On exit from CMVCH YT holds
C     the result computed by CMVCH.
      TRANS = 'N'
      FTL = .FALSE.
      CALL CMVCH( TRANS, N, N, ONE, A, NMAX, X, 1, ZERO, Y, 1, YT, G,
     $            YY, EPS, ERR, FTL, NOUT, .TRUE., KPRINT )
      SAME = LCE( YY, YT, N )
      IF( .NOT.SAME.OR.ERR.NE.RZERO )THEN
        IPASS = 0
        IF (KPRINT .GE. 2) THEN
           WRITE( NOUT, FMT = 9985 )TRANS, SAME, ERR
        END IF
      ENDIF
      TRANS = 'T'
      FTL = .FALSE.
      CALL CMVCH( TRANS, N, N, ONE, A, NMAX, X, -1, ZERO, Y, -1, YT, G,
     $            YY, EPS, ERR, FTL, NOUT, .TRUE., KPRINT )
      SAME = LCE( YY, YT, N )
      IF( .NOT.SAME.OR.ERR.NE.RZERO )THEN
        IPASS = 0
        IF (KPRINT .GE. 2 ) THEN
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
               CALL CCHKE2(ISNUM, SNAMES( ISNUM ), NOUT, KPRINT, FTL1)
            END IF
C           Test computations.
            FTL2 = .FALSE.
            CALL XERCLR
            GO TO ( 140, 140, 150, 150, 150, 160, 160,
     $              160, 160, 160, 160, 170, 170, 180,
     $              180, 190, 190 )ISNUM
C           Test CGEMV, 01, and CGBMV, 02.
  140       CALL CCHK12( SNAMES( ISNUM ), EPS, THRESH, NOUT, KPRINT,
     $                  FTL2, NIDIM, IDIM, NKB, KB, NALF, ALF,
     $                  NBET, BET, NINC, INC, NMAX, INCMAX, A, AA, AS,
     $                  X, XX, XS, Y, YY, YS, YT, G )
            GO TO 200
C           Test CHEMV, 03, CHBMV, 04, and CHPMV, 05.
  150       CALL CCHK22( SNAMES( ISNUM ), EPS, THRESH, NOUT, KPRINT,
     $                  FTL2, NIDIM, IDIM, NKB, KB, NALF, ALF,
     $                  NBET, BET, NINC, INC, NMAX, INCMAX, A, AA, AS,
     $                  X, XX, XS, Y, YY, YS, YT, G )
            GO TO 200
C           Test CTRMV, 06, CTBMV, 07, CTPMV, 08,
C           CTRSV, 09, CTBSV, 10, and CTPSV, 11.
  160       CALL CCHK32( SNAMES( ISNUM ), EPS, THRESH, NOUT, KPRINT,
     $                  FTL2, NIDIM, IDIM, NKB, KB, NINC, INC,
     $                  NMAX, INCMAX, A, AA, AS, Y, YY, YS, YT, G, Z )
            GO TO 200
C           Test CGERC, 12, CGERU, 13.
  170       CALL CCHK42( SNAMES( ISNUM ), EPS, THRESH, NOUT, KPRINT,
     $                  FTL2, NIDIM, IDIM, NALF, ALF, NINC, INC,
     $                  NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS,
     $                  YT, G, Z )
            GO TO 200
C           Test CHER, 14, and CHPR, 15.
  180       CALL CCHK52( SNAMES( ISNUM ), EPS, THRESH, NOUT, KPRINT,
     $                  FTL2, NIDIM, IDIM, NALF, ALF, NINC, INC,
     $                  NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS,
     $                  YT, G, Z )
            GO TO 200
C           Test CHER2, 16, and CHPR2, 17.
  190       CALL CCHK62( SNAMES( ISNUM ), EPS, THRESH, NOUT, KPRINT,
     $                  FTL2, NIDIM, IDIM, NALF, ALF, NINC, INC,
     $                  NMAX, INCMAX, A, AA, AS, X, XX, XS, Y, YY, YS,
     $                  YT, G, Z )
C
  200       IF (FTL1 .OR. FTL2) THEN
              IPASS = 0
            ENDIF
         END IF
  210 CONTINUE
      RETURN
C
 9999 FORMAT( ' ROUTINES PASS COMPUTATIONAL TESTS IF TEST RATIO IS LES',
     $      'S THAN', F8.2 )
 9993 FORMAT( ' TESTS OF THE COMPLEX          LEVEL 2 BLAS', //' THE F',
     $      'OLLOWING PARAMETER VALUES WILL BE USED:' )
 9992 FORMAT( '   FOR N              ', 9I6 )
 9991 FORMAT( '   FOR K              ', 7I6 )
 9990 FORMAT( '   FOR INCX AND INCY  ', 7I6 )
 9989 FORMAT( '   FOR ALPHA          ',
     $      7( '(', F4.1, ',', F4.1, ')  ', : ) )
 9988 FORMAT( '   FOR BETA           ',
     $      7( '(', F4.1, ',', F4.1, ')  ', : ) )
 9985 FORMAT( ' ERROR IN CMVCH -  IN-LINE DOT PRODUCTS ARE BEING EVALU',
     $      'ATED WRONGLY.', /' CMVCH WAS CALLED WITH TRANS = ', A1,
     $      ' AND RETURNED SAME = ', L1, ' AND ERR = ', F12.3, '. ', /
     $   'THIS MAY BE DUE TO FAULTS IN THE ARITHMETIC OR THE COMPILER.')
 9983 FORMAT( 1X, A6, ' WAS NOT TESTED' )
 9980 FORMAT( ' ERROR-EXITS WILL NOT BE TESTED' )
C
C     End of CBLAT2.
C
      END
