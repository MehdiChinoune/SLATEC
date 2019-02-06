*DECK TESTI
      PROGRAM TESTI
C***BEGIN PROLOGUE  TESTI
C***PURPOSE  Driver for testing SLATEC subprogram
C            Fullerton intrinsics.
C***LIBRARY   FNLIB
C***CATEGORY  Z
C***TYPE      ALL (TESTI-A)
C***KEYWORDS  FULLERTON INTRINSIC FUNCTIONS, QUICK CHECK DRIVER
C***AUTHOR  SLATEC Common Mathematical Library Committee
C***DESCRIPTION
C
C *Usage:
C     One input data record is required
C         READ (UNIT=LIN, FMT='(I1)') KPRINT
C
C *Arguments:
C     KPRINT = 0  Quick checks - No printing.
C                 Driver       - Short pass or fail message printed.
C              1  Quick checks - No message printed for passed tests,
C                                short message printed for failed tests.
C                 Driver       - Short pass or fail message printed.
C              2  Quick checks - Print short message for passed tests,
C                                fuller information for failed tests.
C                 Driver       - Pass or fail message printed.
C              3  Quick checks - Print complete quick check results.
C                 Driver       - Pass or fail message printed.
C
C *Description:
C     Driver for testing SLATEC subprogram
C
C***REFERENCES  Fong, Kirby W., Jefferson, Thomas H., Suyehiro,
C                 Tokihiko, Walton, Lee, Guidelines to the SLATEC Common
C                 Mathematical Library, March 21, 1989.
C***ROUTINES CALLED  I1MACH, QCINTC, QCINTD, QCINTS, XERMAX, XSETF,
C                    XSETUN
C***REVISION HISTORY  (YYMMDD)
C   900709  DATE WRITTEN
C***END PROLOGUE  TESTI
C     .. Local Scalars ..
      INTEGER IPASS, KPRINT, LIN, LUN, NFAIL
C     .. External Functions ..
      INTEGER I1MACH
      EXTERNAL I1MACH
C     .. External Subroutines ..
      EXTERNAL QCINTC, QCINTD, QCINTS, XERMAX, XSETF, XSETUN
C***FIRST EXECUTABLE STATEMENT  TESTI
      LUN = I1MACH(2)
      LIN = I1MACH(1)
      NFAIL = 0
C
C     Read KPRINT parameter
C
      READ (UNIT=LIN, FMT='(I1)') KPRINT
      CALL XSETUN (LUN)
      CALL XERMAX (1000)
      IF (KPRINT .LE. 1) THEN
         CALL XSETF (0)
      ELSE
         CALL XSETF (1)
      ENDIF
C
C     Test single precision Fullerton intrinsics.
C
      CALL QCINTS (LUN, KPRINT, IPASS)
      IF (IPASS .EQ. 0) NFAIL = NFAIL + 1
C
C     Test double precision Fullerton intrinsics.
C
      CALL QCINTD (LUN, KPRINT, IPASS)
      IF (IPASS .EQ. 0) NFAIL = NFAIL + 1
C
C     Test complex Fullerton intrinsics.
C
      CALL QCINTC (LUN, KPRINT, IPASS)
      IF (IPASS .EQ. 0) NFAIL = NFAIL + 1
C
C     Write PASS or FAIL message
C
      IF (NFAIL .EQ. 0) THEN
         WRITE (UNIT=LUN, FMT=9000)
      ELSE
         WRITE (UNIT=LUN, FMT=9010) NFAIL
      ENDIF
      STOP
 9000 FORMAT (/' --------------TESTI PASSED ALL TESTS----------------')
 9010 FORMAT (/' ************* WARNING -- ', I5,
     1        ' TEST(S) FAILED IN PROGRAM TESTI *************')
      END
*DECK QCINTC
      SUBROUTINE QCINTC (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  QCINTC
C***PURPOSE  Quick check for the complex Fullerton elementary
C            intrinsic functions.
C***LIBRARY   FNLIB
C***CATEGORY  C
C***TYPE      COMPLEX (QCINTS-S, QCINTD-D, QCINTC-C)
C***KEYWORDS  QUICK CHECK
C***AUTHOR  Boland, W. Robert, (LANL)
C           Rivera, Shawn M., (LANL)
C***DESCRIPTION
C
C   This subroutine does a quick check for the complex
C   Fullerton elementary intrinsic functions.
C
C   Parameter list-
C
C   LUN      input INTEGER value to designate the external device unit
C            for message output
C   KPRINT   input INTEGER value to specify amount of printing to be
C            done by quick check
C   IPASS    output INTEGER value indicating whether tests passed or
C            failed
C
C***ROUTINES CALLED  CABS, CCOS, CEXP, CLOG, CSIN, CSQRT, R1MACH, SQRT
C***REVISION HISTORY  (YYMMDD)
C   900717  DATE WRITTEN
C***END PROLOGUE  QCINTC
C     .. Scalar Arguments ..
      INTEGER IPASS, KPRINT, LUN
C     .. Local Scalars ..
      REAL ERRTOL
      INTEGER I
C     .. Local Arrays ..
      COMPLEX C(20), W(20)
C     .. External Functions ..
      COMPLEX CCOS, CEXP, CLOG, CSIN, CSQRT
      REAL CABS, R1MACH, SQRT
      EXTERNAL CCOS, CEXP, CLOG, CSIN, CSQRT, CABS, R1MACH, SQRT
C     .. Intrinsic Functions ..
      INTRINSIC CMPLX
C
C     Complex values through different calculations are stored in C(*)
C
C     .. Data statements ..
      DATA C( 1) /(  1.0000000000000, 0.0000000000000) /
      DATA C( 2) /(  89.00280929194, .0078649202825041) /
      DATA C( 3) /(  30.00001041666, .024999991319455) /
      DATA C( 4) /(  6324555.320337, .0000001897366596101) /
      DATA C( 5) /( -0.8414709848079, 0.0000000000000) /
      DATA C( 6) /(  27.23982534694, 1.930412376268) /
      DATA C( 7) /(  0.000000000000000, 1.175201193644) /
      DATA C( 8) /(  1.127805246806, 1.868618519183) /
      DATA C( 9) /(  0.5403023058681, 0.0000000000000) /
      DATA C(10) /(  23.96522893293, 13.0834832507) /
      DATA C(11) /(  1.543080634815, 0.00000000000000) /
      DATA C(12) /(  2.064433656761, -1.020830949598) /
      DATA C(13) /( -2.929427471521, -3.391753471626) /
      DATA C(14) /( -0.7373937155412, 0.6754631805511) /
      DATA C(15) /(  .1699671429002, .9854497299884) /
      DATA C(16) /(  0.7055457557766, 9.949196994152) /
      DATA C(17) /(  3.738352258649, 0.3119690755436) /
      DATA C(18) /(  4.605747852161, .033986907746255) /
      DATA C(19) /(  2.313710397461, 0.1488899476095) /
      DATA C(20) /(  6.907755278982, 0.00000000000000) /
C
C***FIRST EXECUTABLE STATEMENT  QCINTC
C
      IF (KPRINT .GE. 2) WRITE (UNIT=LUN, FMT=9000)
C
C     Exercise routines in Category C2.
C
      W( 1) = CSQRT(CMPLX(1.0, 0.0))
      W( 2) = CSQRT(CMPLX(7921.5, 1.4))
      W( 3) = CSQRT(CMPLX(900.0, 1.5))
      W( 4) = CSQRT(CMPLX(0.4E+14, 2.4))
C
C     Exercise routines in Category C4A.
C
      W( 5) = CSIN(CMPLX(-1.0, 0.0))
      W( 6) = CSIN(CMPLX(1.5, 4.0))
      W( 7) = CSIN(CMPLX(0.0, 1.0))
      W( 8) = CSIN(CMPLX(0.5, 1.5))
      W( 9) = CCOS(CMPLX(-1.0, 0.0))
      W(10) = CCOS(CMPLX(-0.5, 4.0))
      W(11) = CCOS(CMPLX(0.0, 1.0))
      W(12) = CCOS(CMPLX(0.5, 1.5))
C
C     Exercise routines in Category C4B.
C
      W(13) = CEXP(CMPLX(1.5, 4.0))
      W(14) = CEXP(CMPLX(0.0, 2.4))
      W(15) = CEXP(CMPLX(0.0, 1.4))
      W(16) = CEXP(CMPLX(2.3, 1.5))
      W(17) = CLOG(CMPLX(40.0, 12.9))
      W(18) = CLOG(CMPLX(100.0, 3.4))
      W(19) = CLOG(CMPLX(10.0, 1.5))
      W(20) = CLOG(CMPLX(1000.0, 0.0))
C
C     Check for possible errors.
C
      IPASS = 1
      ERRTOL = SQRT(R1MACH(4))
      DO 10 I = 1,20
        IF (CABS(C(I)-W(I)) .GE. ERRTOL*CABS(C(I))+ERRTOL) THEN
          IPASS = 0
          IF (KPRINT .GE. 2) WRITE (UNIT=LUN, FMT=9020) I, W(I), C(I)
        ENDIF
   10 CONTINUE
      IF (IPASS.NE.0 .AND. KPRINT.GE.2) WRITE (UNIT=LUN, FMT=9010)
      RETURN
 9000 FORMAT (// ' Test of complex Fullerton intrinsic routines')
 9010 FORMAT (' Complex Fullerton intrinsic function routines o.k.')
 9020 FORMAT (' For I  = ', I3, '  test fails with  ', /
     +        ' computed result = (', 1P, E22.14, ', ', E22.14,'  ) '/
     +        ' and true result = (', E22.14, ', ', E22.14, '  )')
      END
*DECK QCINTD
      SUBROUTINE QCINTD (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  QCINTD
C***PURPOSE  Quick check for the double precision Fullerton
C            elementary intrinsic functions.
C***LIBRARY   FNLIB
C***CATEGORY  C
C***TYPE      DOUBLE PRECISION (QCINTS-S, QCINTD-D, QCINTC-C)
C***KEYWORDS  QUICK CHECK
C***AUTHOR  Boland, W. Robert, (LANL)
C           Rivera, Shawn M., (LANL)
C***DESCRIPTION
C
C   This subroutine does a quick check for the double precision
C   Fullerton intrinsic functions.
C
C   Parameter list-
C
C   LUN      input INTEGER value to designate the external device unit
C            for message output
C   KPRINT   input INTEGER value to specify amount of printing to be
C            done by quick check
C   IPASS    output INTEGER value indicating whether tests passed or
C            failed
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DACOS, DASIN, DATAN, DATAN2, DCOS, DCOSH,
C                    DEXP, DINT, DLOG, DLOG10, DSIN, DSINH, DSQRT, DTAN,
C                    DTANH
C***REVISION HISTORY  (YYMMDD)
C   900717  DATE WRITTEN
C***END PROLOGUE  QCINTD
C     .. Scalar Arguments ..
      INTEGER IPASS, KPRINT, LUN
C     .. Local Scalars ..
      DOUBLE PRECISION ERRTOL
      INTEGER I
C     .. Local Arrays ..
      DOUBLE PRECISION V(60), Y(60)
C     .. External Functions ..
      DOUBLE PRECISION D1MACH, DACOS, DASIN, DATAN, DATAN2, DCOS, DCOSH,
     +                 DEXP, DINT, DLOG, DLOG10, DSIN, DSINH, DSQRT,
     +                 DTAN, DTANH
      EXTERNAL D1MACH, DACOS, DASIN, DATAN, DATAN2, DCOS, DCOSH, DEXP,
     +                 DINT, DLOG, DLOG10, DSIN, DSINH, DSQRT, DTAN,
     +                 DTANH
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C
C     Correct values through different calculations are stored in V(*)
C
C     .. Data statements ..
      DATA V( 1) /  10.0D0 /
      DATA V( 2) /  79.0D0 /
      DATA V( 3) /  900.0D0 /
      DATA V( 4) /  4.0D0 /
      DATA V( 5) /  1.0D0 /
      DATA V( 6) /  89.0D0 /
      DATA V( 7) /  30.0D0 /
      DATA V( 8) /  6.32455532033675866399778708D06 /
      DATA V( 9) /  3.1415926535897932846264338D0 /
      DATA V(10) /  2.09439510239319549230842892D0 /
      DATA V(11) /  1.57079632679489661923132169D0 /
      DATA V(12) /  1.04719755119659774615421446D0 /
      DATA V(13) / -1.57079632679489661923132169D0 /
      DATA V(14) / -0.52359877559829887307710723D0 /
      DATA V(15) /  0.0D0 /
      DATA V(16) /  0.52359877559829887307710723D0 /
      DATA V(17) /  -0.785398163397448309615660845D0 /
      DATA V(18) / -0.463647609000806116214256231D0 /
      DATA V(19) /  0.0D0 /
      DATA V(20) /  0.463647609000806116214256231D0 /
      DATA V(21) / -0.58800260354756755124561108D0 /
      DATA V(22) / -0.463647609000806116214256231D0 /
      DATA V(23) /  2.034443935795702707025611744029D0 /
      DATA V(24) /  2.158798930342464394982471276307D0 /
      DATA V(25) /  0.540302305868139717400936607D0 /
      DATA V(26) /  0.877582561890372716116281582D0 /
      DATA V(27) /  1.0D0 /
      DATA V(28) /  0. 877582561890372716116281582D0 /
      DATA V(29) / -0.841470984807896506652502321D0 /
      DATA V(30) / -0.479425538604203000273287935D0 /
      DATA V(31) /  0.0D0 /
      DATA V(32) /  0.479425538604203000273287935D0 /
      DATA V(33) /  -1.55740772465490223050697485D0 /
      DATA V(34) / -0.546302489843790513255179465D0 /
      DATA V(35) /  0.0D0 /
      DATA V(36) /  0.546302489843790513255179465D0 /
      DATA V(37) /  2.30258509299404568401799145D0 /
      DATA V(38) /  2.99573227355399099343522357D0 /
      DATA V(39) /  3.40119738166215537541323669D0 /
      DATA V(40) /  3.68887945411393630285245569D0 /
      DATA V(41) /  1.0D0 /
      DATA V(42) /  1.30102999566398119521373889D0 /
      DATA V(43) /  1.4771212547196624372950279D0 /
      DATA V(44) /  1.60205999132796239042747778D0 /
      DATA V(45) /  1.00000100530050531421637777D0 /
      DATA V(46) /  0.999843012323855043126609044D0 /
      DATA V(47) /  1.00003876575137232151808428D0 /
      DATA V(48) /  0.992002154326025434343372944D0 /
      DATA V(49) /  1.54308063481524377847790562D0 /
      DATA V(50) /  1.12762596520638078522622516D0 /
      DATA V(51) /  1.0D0 /
      DATA V(52) /  1.12762596520638078522622516D0 /
      DATA V(53) / -1.175201193643801456882381851D0 /
      DATA V(54) / -0.521095305493747361622425626D0 /
      DATA V(55) /  0.0D0 /
      DATA V(56) /  0.521095305493747361622425626D0 /
      DATA V(57) / -0.761594155955764888119458282D0 /
      DATA V(58) / -0.462117157260009758502318483D0 /
      DATA V(59) /  0.0D0 /
      DATA V(60) /  0.462117157260009758592318483D0 /
C
C***FIRST EXECUTABLE STATEMENT  QCINTD
C
      IF (KPRINT .GE. 2) WRITE (UNIT=LUN, FMT=9000)
C
C     Exercise routines in Category C1.
C
      Y( 1) = DINT(10.465890D0)
      Y( 2) = DINT(79.32178D0)
      Y( 3) = DINT(900.0D0)
      Y( 4) = DINT(4.0D0)
C
C     Exercise routines in Category C2.
C
      Y( 5) = DSQRT(1.0D0)
      Y( 6) = DSQRT(7921.0D0)
      Y( 7) = DSQRT(900.0D0)
      Y( 8) = DSQRT(4000D+10)
C
C     Exercise routines in Category C4A.
C
      Y( 9) = DACOS(-1.0D0)
      Y(10) = DACOS(-0.5D0)
      Y(11) = DACOS(0.0D0)
      Y(12) = DACOS(0.5D0)
      Y(13) = DASIN(-1.0D0)
      Y(14) = DASIN(-0.5D0)
      Y(15) = DASIN(0.0D0)
      Y(16) = DASIN(0.5D0)
      Y(17) = DATAN(-1.0D0)
      Y(18) = DATAN(-0.5D0)
      Y(19) = DATAN(0.0D0)
      Y(20) = DATAN(0.5D0)
      Y(21) = DATAN2(-1.0D0,1.5D0)
      Y(22) = DATAN2(-0.5D0,1.0D0)
      Y(23) = DATAN2(1.0D0,-0.5D0)
      Y(24) = DATAN2(1.5D0,-1.0D0)
      Y(25) = DCOS(-1.0D0)
      Y(26) = DCOS(-0.5D0)
      Y(27) = DCOS(0.0D0)
      Y(28) = DCOS(0.5D0)
      Y(29) = DSIN(-1.0D0)
      Y(30) = DSIN(-0.5D0)
      Y(31) = DSIN(0.0D0)
      Y(32) = DSIN(0.5D0)
      Y(33) = DTAN(-1.0D0)
      Y(34) = DTAN(-0.5D0)
      Y(35) = DTAN(0.0D0)
      Y(36) = DTAN(0.5D0)
C
C     Exercise routines in Category C4B.
C
      Y(37) = DLOG(10.0D0)
      Y(38) = DLOG(20.0D0)
      Y(39) = DLOG(30.0D0)
      Y(40) = DLOG(40.0D0)
      Y(41) = DLOG10(10.0D0)
      Y(42) = DLOG10(20.0D0)
      Y(43) = DLOG10(30.0D0)
      Y(44) = DLOG10(40.0D0)
      Y(45) = DEXP(1.0053D-06)
      Y(46) = DEXP(-1.57D-04)
      Y(47) = DEXP(3.8765D-05)
      Y(48) = DEXP(-8.03D-03)
C
C     Exercise routines in Category C4C.
C
      Y(49) = DCOSH(-1.0D0)
      Y(50) = DCOSH(-0.5D0)
      Y(51) = DCOSH(0.0D0)
      Y(52) = DCOSH(0.5D0)
      Y(53) = DSINH(-1.0D0)
      Y(54) = DSINH(-0.5D0)
      Y(55) = DSINH(0.0D0)
      Y(56) = DSINH(0.5D0)
      Y(57) = DTANH(-1.0D0)
      Y(58) = DTANH(-0.5D0)
      Y(59) = DTANH(0.0D0)
      Y(60) = DTANH(0.5D0)
C
C     Check for possible errors.
C
      IPASS = 1
      ERRTOL = DSQRT(D1MACH(4))
      DO 10 I = 1,60
        IF (ABS(V(I)-Y(I)) .GE. ERRTOL*ABS(V(I))+ERRTOL) THEN
          IPASS = 0
          IF (KPRINT .GE. 2) WRITE (UNIT=LUN, FMT=9020) I, Y(I), V(I)
        ENDIF
   10 CONTINUE
      IF (IPASS.NE.0 .AND. KPRINT.GE.2) WRITE (UNIT=LUN, FMT=9010)
      RETURN
 9000 FORMAT (// ' Test of double precision Fullerton intrinsic ',
     +        'routines')
 9010 FORMAT (' Double precision Fullerton intrinsic function ',
     +        'routines o.k.')
 9020 FORMAT (' For I  = ', I3, '  test fails with ', /
     +        ' computed result = ', 1P, E38.30, /
     +        ' and true result = ', E38.30)
      END
*DECK QCINTS
      SUBROUTINE QCINTS (LUN, KPRINT, IPASS)
C***BEGIN PROLOGUE  QCINTS
C***PURPOSE  Quick check for the single precision Fullerton
C            elementary intrinsic functions.
C***LIBRARY   FNLIB
C***CATEGORY  C
C***TYPE      SINGLE PRECISION (QCINTS-S, QCINTD-D, QCTINC-C)
C***KEYWORDS  QUICK CHECK
C***AUTHOR  Boland, W. Robert, (LANL)
C           Rivera, Shawn M., (LANL)
C***DESCRIPTION
C
C   This subroutine does a quick check for the single precision
C   Fullerton intrinsic functions.
C
C   Parameter list-
C
C   LUN      input INTEGER value to designate the external device unit
C            for message output
C   KPRINT   input INTEGER value to specify amount of printing to be
C            done by quick check
C   IPASS    output INTEGER value indicating whether tests passed or
C            failed
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  ACOS, ALOG, ALOG10, ASIN, ATAN, ATAN2, CABS, COS,
C                    COSH, EXP, R1MACH, SIN, SINH, SQRT, TAN, TANH
C***REVISION HISTORY  (YYMMDD)
C   900711  DATE WRITTEN
C***END PROLOGUE  QCINTS
C     .. Scalar Arguments ..
      INTEGER IPASS, KPRINT, LUN
C     .. Local Scalars ..
      REAL ERRTOL
      INTEGER I
C     .. Local Arrays ..
      REAL V(60), Y(60)
C     .. External Functions ..
      REAL ACOS, ALOG, ALOG10, ASIN, ATAN, ATAN2, CABS, COS, COSH, EXP,
     +     R1MACH, SIN, SINH, SQRT, TAN, TANH
      EXTERNAL ACOS, ALOG, ALOG10, ASIN, ATAN, ATAN2, CABS, COS, COSH,
     +         EXP, R1MACH, SIN, SINH, SQRT, TAN, TANH
C     .. Intrinsic Functions ..
      INTRINSIC ABS, CMPLX
C
C     Correct values through different calculations are stored in V(*)
C
C     .. Data statements ..
      DATA V( 1) /  1.0 /
      DATA V( 2) /  89.0 /
      DATA V( 3) /  30.0 /
      DATA V( 4) /  6.324555320337E+06 /
      DATA V( 5) /  10.55327437339 /
      DATA V( 6) /  79.32157587945 /
      DATA V( 7) /  901.0429556913 /
      DATA V( 8) /  4.00000E+13 /
      DATA V( 9) /  3.14159265359 /
      DATA V(10) /  2.094395102393 /
      DATA V(11) /  1.570796326795 /
      DATA V(12) /  1.047197551197 /
      DATA V(13) / -1.570796326795 /
      DATA V(14) / -0.5235987755983 /
      DATA V(15) /  0.0 /
      DATA V(16) /  0.5235987755983 /
      DATA V(17) / -0.7853981633974 /
      DATA V(18) / -0.4636476090008 /
      DATA V(19) /  0.0 /
      DATA V(20) /  0.4636476090008 /
      DATA V(21) / -0.5880026035475 /
      DATA V(22) / -0.4636476090008 /
      DATA V(23) /  2.0344438552856 /
      DATA V(24) /  2.158798930342 /
      DATA V(25) /  0.5403023058681 /
      DATA V(26) /  0.8775825618903 /
      DATA V(27) /  1.0 /
      DATA V(28) /  0.8775825618903 /
      DATA V(29) / -0.8414709848079 /
      DATA V(30) / -0.4794255386042 /
      DATA V(31) /  0.0 /
      DATA V(32) /  0.4794255386042 /
      DATA V(33) / -1.557407724655 /
      DATA V(34) / -0.5463024898437 /
      DATA V(35) /  0.0 /
      DATA V(36) /  0.5463024898437 /
      DATA V(37) /  2.302585092994 /
      DATA V(38) /  2.995732273554 /
      DATA V(39) /  3.401197381662 /
      DATA V(40) /  3.688879454114 /
      DATA V(41) /  1.0 /
      DATA V(42) /  1.301029995664 /
      DATA V(43) /  1.47712125472 /
      DATA V(44) /  1.602059991328 /
      DATA V(45) /  1.000001005301 /
      DATA V(46) /  0.9998430123238 /
      DATA V(47) /  1.000038765751 /
      DATA V(48) /  0.992002154326 /
      DATA V(49) /  1.543080634815 /
      DATA V(50) /  1.127625965206 /
      DATA V(51) /  1.0 /
      DATA V(52) /  1.127625965206 /
      DATA V(53) / -1.175201193644 /
      DATA V(54) / -0.5210953054937 /
      DATA V(55) /  0.0 /
      DATA V(56) /  0.5210953054937 /
      DATA V(57) / -0.7615941559557 /
      DATA V(58) / -0.46211715726 /
      DATA V(59) /  0.0 /
      DATA V(60) /  0.46211715726 /
C
C***FIRST EXECUTABLE STATEMENT  QCINTS
C
      IF (KPRINT .GE. 2) WRITE (UNIT=LUN, FMT=9000)
C
C     Exercise routines in Category C2.
C
      Y( 1) = SQRT(1.0)
      Y( 2) = SQRT(7921.0)
      Y( 3) = SQRT(900.0)
      Y( 4) = SQRT(4.00000E+13)
C
C     Exercise routines in Category C4.
C
      Y( 5) = CABS(CMPLX(10.46,1.4))
      Y( 6) = CABS(CMPLX(79.32,0.5))
      Y( 7) = CABS(CMPLX(900.999,8.9))
      Y( 8) = CABS(CMPLX(4.00000E+13,1.5))
C
C     Exercise routines in Category C4A.
C
      Y( 9) = ACOS(-1.0)
      Y(10) = ACOS(-0.5)
      Y(11) = ACOS(0.0)
      Y(12) = ACOS(0.5)
      Y(13) = ASIN(-1.0)
      Y(14) = ASIN(-0.5)
      Y(15) = ASIN(0.0)
      Y(16) = ASIN(0.5)
      Y(17) = ATAN(-1.0)
      Y(18) = ATAN(-0.5)
      Y(19) = ATAN(0.0)
      Y(20) = ATAN(0.5)
      Y(21) = ATAN2(-1.0,1.5)
      Y(22) = ATAN2(-0.5,1.0)
      Y(23) = ATAN2(1.0,-0.5)
      Y(24) = ATAN2(1.5,-1.0)
      Y(25) = COS(-1.0)
      Y(26) = COS(-0.5)
      Y(27) = COS(0.0)
      Y(28) = COS(0.5)
      Y(29) = SIN(-1.0)
      Y(30) = SIN(-0.5)
      Y(31) = SIN(0.0)
      Y(32) = SIN(0.5)
      Y(33) = TAN(-1.0)
      Y(34) = TAN(-0.5)
      Y(35) = TAN(0.0)
      Y(36) = TAN(0.5)
C
C     Exercise routines in Category C4B.
C
      Y(37) = ALOG(10.0)
      Y(38) = ALOG(20.0)
      Y(39) = ALOG(30.0)
      Y(40) = ALOG(40.0)
      Y(41) = ALOG10(10.0)
      Y(42) = ALOG10(20.0)
      Y(43) = ALOG10(30.0)
      Y(44) = ALOG10(40.0)
      Y(45) = EXP(1.0053E-06)
      Y(46) = EXP(-1.57000E-04)
      Y(47) = EXP(3.87650E-05)
      Y(48) = EXP(-8.03000E-03)
C
C     Exercise routines in Category C4C.
C
      Y(49) = COSH(-1.0)
      Y(50) = COSH(-0.5)
      Y(51) = COSH(0.0)
      Y(52) = COSH(0.5)
      Y(53) = SINH(-1.00000)
      Y(54) = SINH(-0.50000)
      Y(55) = SINH(0.000000)
      Y(56) = SINH(0.500000)
      Y(57) = TANH(-1.00000)
      Y(58) = TANH(-0.50000)
      Y(59) = TANH(0.000000)
      Y(60) = TANH(0.500000)
C
C     Check for possible errors.
C
      IPASS = 1
      ERRTOL = SQRT(R1MACH(4))
      DO 10 I = 1,60
        IF (ABS(V(I)-Y(I)) .GE. ERRTOL*ABS(V(I))+ERRTOL) THEN
          IPASS = 0
          IF (KPRINT .GE. 2) WRITE (UNIT=LUN, FMT=9020) I, Y(I), V(I)
        ENDIF
   10 CONTINUE
      IF (IPASS.NE.0 .AND. KPRINT.GE.2) WRITE (UNIT=LUN, FMT=9010)
      RETURN
 9000 FORMAT (// ' Test of single precision Fullerton intrinsic ',
     +        'routines')
 9010 FORMAT (' Single precision Fullerton intrinsic function ',
     +        'routines o.k.')
 9020 FORMAT (' For I  = ', I3, '  test fails with ', /
     +        ' computed result = ', 1P, E22.14, /
     +        ' and true result = ', E22.14)
      END
