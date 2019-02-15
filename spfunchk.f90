!DECK TESTI
PROGRAM TESTI
  IMPLICIT NONE
  !***BEGIN PROLOGUE  TESTI
  !***PURPOSE  Driver for testing SLATEC subprogram
  !            Fullerton intrinsics.
  !***LIBRARY   FNLIB
  !***CATEGORY  Z
  !***TYPE      ALL (TESTI-A)
  !***KEYWORDS  FULLERTON INTRINSIC FUNCTIONS, QUICK CHECK DRIVER
  !***AUTHOR  SLATEC Common Mathematical Library Committee
  !***DESCRIPTION
  !
  ! *Usage:
  !     One input data record is required
  !         READ (UNIT=LIN, FMT='(I1)') KPRINT
  !
  ! *Arguments:
  !     KPRINT = 0  Quick checks - No printing.
  !                 Driver       - Short pass or fail message printed.
  !              1  Quick checks - No message printed for passed tests,
  !                                short message printed for failed tests.
  !                 Driver       - Short pass or fail message printed.
  !              2  Quick checks - Print short message for passed tests,
  !                                fuller information for failed tests.
  !                 Driver       - Pass or fail message printed.
  !              3  Quick checks - Print complete quick check results.
  !                 Driver       - Pass or fail message printed.
  !
  ! *Description:
  !     Driver for testing SLATEC subprogram
  !
  !***REFERENCES  Fong, Kirby W., Jefferson, Thomas H., Suyehiro,
  !                 Tokihiko, Walton, Lee, Guidelines to the SLATEC Common
  !                 Mathematical Library, March 21, 1989.
  !***ROUTINES CALLED  I1MACH, QCINTC, QCINTD, QCINTS, XERMAX, XSETF,
  !                    XSETUN
  !***REVISION HISTORY  (YYMMDD)
  !   900709  DATE WRITTEN
  !***END PROLOGUE  TESTI
  !     .. Local Scalars ..
  INTEGER ipass, kprint, lin, lun, nfail
  !     .. External Functions ..
  INTEGER I1MACH
  EXTERNAL I1MACH
  !     .. External Subroutines ..
  EXTERNAL QCINTC, QCINTD, QCINTS, XERMAX, XSETF, XSETUN
  !***FIRST EXECUTABLE STATEMENT  TESTI
  lun = I1MACH(2)
  lin = I1MACH(1)
  nfail = 0
  !
  !     Read KPRINT parameter
  !
  READ (UNIT=lin,FMT='(I1)') kprint
  CALL XSETUN(lun)
  CALL XERMAX(1000)
  IF ( kprint<=1 ) THEN
    CALL XSETF(0)
  ELSE
    CALL XSETF(1)
  ENDIF
  !
  !     Test single precision Fullerton intrinsics.
  !
  CALL QCINTS(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test double precision Fullerton intrinsics.
  !
  CALL QCINTD(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Test complex Fullerton intrinsics.
  !
  CALL QCINTC(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF ( nfail==0 ) THEN
    WRITE (UNIT=lun,FMT=99001)
  ELSE
    WRITE (UNIT=lun,FMT=99002) nfail
  ENDIF
  STOP
  99001 FORMAT (/' --------------TESTI PASSED ALL TESTS----------------')
  99002 FORMAT (/' ************* WARNING -- ',I5,&
    ' TEST(S) FAILED IN PROGRAM TESTI *************')
END PROGRAM TESTI
!DECK QCINTC
SUBROUTINE QCINTC(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  QCINTC
  !***PURPOSE  Quick check for the complex Fullerton elementary
  !            intrinsic functions.
  !***LIBRARY   FNLIB
  !***CATEGORY  C
  !***TYPE      COMPLEX (QCINTS-S, QCINTD-D, QCINTC-C)
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  Boland, W. Robert, (LANL)
  !           Rivera, Shawn M., (LANL)
  !***DESCRIPTION
  !
  !   This subroutine does a quick check for the complex
  !   Fullerton elementary intrinsic functions.
  !
  !   Parameter list-
  !
  !   LUN      input INTEGER value to designate the external device unit
  !            for message output
  !   KPRINT   input INTEGER value to specify amount of printing to be
  !            done by quick check
  !   IPASS    output INTEGER value indicating whether tests passed or
  !            failed
  !
  !***ROUTINES CALLED  CABS, CCOS, CEXP, CLOG, CSIN, CSQRT, R1MACH, SQRT
  !***REVISION HISTORY  (YYMMDD)
  !   900717  DATE WRITTEN
  !***END PROLOGUE  QCINTC
  !     .. Scalar Arguments ..
  INTEGER Ipass, Kprint, Lun
  !     .. Local Scalars ..
  REAL errtol
  INTEGER i
  !     .. Local Arrays ..
  COMPLEX c(20), w(20)
  !     .. External Functions ..
  COMPLEX CCOS, CEXP, CLOG, CSIN, CSQRT
  REAL CABS, R1MACH, SQRT
  EXTERNAL CCOS, CEXP, CLOG, CSIN, CSQRT, CABS, R1MACH, SQRT
  !     .. Intrinsic Functions ..
  INTRINSIC CMPLX
  !
  !     Complex values through different calculations are stored in C(*)
  !
  !     .. Data statements ..
  DATA c(1)/(1.0000000000000,0.0000000000000)/
  DATA c(2)/(89.00280929194,.0078649202825041)/
  DATA c(3)/(30.00001041666,.024999991319455)/
  DATA c(4)/(6324555.320337,.0000001897366596101)/
  DATA c(5)/(-0.8414709848079,0.0000000000000)/
  DATA c(6)/(27.23982534694,1.930412376268)/
  DATA c(7)/(0.000000000000000,1.175201193644)/
  DATA c(8)/(1.127805246806,1.868618519183)/
  DATA c(9)/(0.5403023058681,0.0000000000000)/
  DATA c(10)/(23.96522893293,13.0834832507)/
  DATA c(11)/(1.543080634815,0.00000000000000)/
  DATA c(12)/(2.064433656761,-1.020830949598)/
  DATA c(13)/(-2.929427471521,-3.391753471626)/
  DATA c(14)/(-0.7373937155412,0.6754631805511)/
  DATA c(15)/(.1699671429002,.9854497299884)/
  DATA c(16)/(0.7055457557766,9.949196994152)/
  DATA c(17)/(3.738352258649,0.3119690755436)/
  DATA c(18)/(4.605747852161,.033986907746255)/
  DATA c(19)/(2.313710397461,0.1488899476095)/
  DATA c(20)/(6.907755278982,0.00000000000000)/
  !
  !***FIRST EXECUTABLE STATEMENT  QCINTC
  !
  IF ( Kprint>=2 ) WRITE (UNIT=Lun,FMT=99001)
  !
  !     Exercise routines in Category C2.
  !
  w(1) = CSQRT(CMPLX(1.0,0.0))
  w(2) = CSQRT(CMPLX(7921.5,1.4))
  w(3) = CSQRT(CMPLX(900.0,1.5))
  w(4) = CSQRT(CMPLX(0.4E+14,2.4))
  !
  !     Exercise routines in Category C4A.
  !
  w(5) = CSIN(CMPLX(-1.0,0.0))
  w(6) = CSIN(CMPLX(1.5,4.0))
  w(7) = CSIN(CMPLX(0.0,1.0))
  w(8) = CSIN(CMPLX(0.5,1.5))
  w(9) = CCOS(CMPLX(-1.0,0.0))
  w(10) = CCOS(CMPLX(-0.5,4.0))
  w(11) = CCOS(CMPLX(0.0,1.0))
  w(12) = CCOS(CMPLX(0.5,1.5))
  !
  !     Exercise routines in Category C4B.
  !
  w(13) = CEXP(CMPLX(1.5,4.0))
  w(14) = CEXP(CMPLX(0.0,2.4))
  w(15) = CEXP(CMPLX(0.0,1.4))
  w(16) = CEXP(CMPLX(2.3,1.5))
  w(17) = CLOG(CMPLX(40.0,12.9))
  w(18) = CLOG(CMPLX(100.0,3.4))
  w(19) = CLOG(CMPLX(10.0,1.5))
  w(20) = CLOG(CMPLX(1000.0,0.0))
  !
  !     Check for possible errors.
  !
  Ipass = 1
  errtol = SQRT(R1MACH(4))
  DO i = 1, 20
    IF ( CABS(c(i)-w(i))>=errtol*CABS(c(i))+errtol ) THEN
      Ipass = 0
      IF ( Kprint>=2 ) WRITE (UNIT=Lun,FMT=99003) i, w(i), c(i)
    ENDIF
  ENDDO
  IF ( Ipass/=0.AND.Kprint>=2 ) WRITE (UNIT=Lun,FMT=99002)
  RETURN
  99001 FORMAT (//' Test of complex Fullerton intrinsic routines')
  99002 FORMAT (' Complex Fullerton intrinsic function routines o.k.')
  99003 FORMAT (' For I  = ',I3,'  test fails with  ',/' computed result = (',1P,&
    E22.14,', ',E22.14,'  ) '/' and true result = (',E22.14,', ',&
    E22.14,'  )')
END SUBROUTINE QCINTC
!DECK QCINTD
SUBROUTINE QCINTD(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  QCINTD
  !***PURPOSE  Quick check for the double precision Fullerton
  !            elementary intrinsic functions.
  !***LIBRARY   FNLIB
  !***CATEGORY  C
  !***TYPE      DOUBLE PRECISION (QCINTS-S, QCINTD-D, QCINTC-C)
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  Boland, W. Robert, (LANL)
  !           Rivera, Shawn M., (LANL)
  !***DESCRIPTION
  !
  !   This subroutine does a quick check for the double precision
  !   Fullerton intrinsic functions.
  !
  !   Parameter list-
  !
  !   LUN      input INTEGER value to designate the external device unit
  !            for message output
  !   KPRINT   input INTEGER value to specify amount of printing to be
  !            done by quick check
  !   IPASS    output INTEGER value indicating whether tests passed or
  !            failed
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  D1MACH, DACOS, DASIN, DATAN, DATAN2, DCOS, DCOSH,
  !                    DEXP, DINT, DLOG, DLOG10, DSIN, DSINH, DSQRT, DTAN,
  !                    DTANH
  !***REVISION HISTORY  (YYMMDD)
  !   900717  DATE WRITTEN
  !***END PROLOGUE  QCINTD
  !     .. Scalar Arguments ..
  INTEGER Ipass, Kprint, Lun
  !     .. Local Scalars ..
  REAL(8) :: errtol
  INTEGER i
  !     .. Local Arrays ..
  REAL(8) :: v(60), y(60)
  !     .. External Functions ..
  REAL(8) :: D1MACH, DACOS, DASIN, DATAN, DATAN2, DCOS, DCOSH, &
    DEXP, DINT, DLOG, DLOG10, DSIN, DSINH, DSQRT, &
    DTAN, DTANH
  EXTERNAL D1MACH, DACOS, DASIN, DATAN, DATAN2, DCOS, DCOSH, DEXP, &
    DINT, DLOG, DLOG10, DSIN, DSINH, DSQRT, DTAN, DTANH
  !     .. Intrinsic Functions ..
  INTRINSIC ABS
  !
  !     Correct values through different calculations are stored in V(*)
  !
  !     .. Data statements ..
  DATA v(1)/10.0D0/
  DATA v(2)/79.0D0/
  DATA v(3)/900.0D0/
  DATA v(4)/4.0D0/
  DATA v(5)/1.0D0/
  DATA v(6)/89.0D0/
  DATA v(7)/30.0D0/
  DATA v(8)/6.32455532033675866399778708D06/
  DATA v(9)/3.1415926535897932846264338D0/
  DATA v(10)/2.09439510239319549230842892D0/
  DATA v(11)/1.57079632679489661923132169D0/
  DATA v(12)/1.04719755119659774615421446D0/
  DATA v(13)/ - 1.57079632679489661923132169D0/
  DATA v(14)/ - 0.52359877559829887307710723D0/
  DATA v(15)/0.0D0/
  DATA v(16)/0.52359877559829887307710723D0/
  DATA v(17)/ - 0.785398163397448309615660845D0/
  DATA v(18)/ - 0.463647609000806116214256231D0/
  DATA v(19)/0.0D0/
  DATA v(20)/0.463647609000806116214256231D0/
  DATA v(21)/ - 0.58800260354756755124561108D0/
  DATA v(22)/ - 0.463647609000806116214256231D0/
  DATA v(23)/2.034443935795702707025611744029D0/
  DATA v(24)/2.158798930342464394982471276307D0/
  DATA v(25)/0.540302305868139717400936607D0/
  DATA v(26)/0.877582561890372716116281582D0/
  DATA v(27)/1.0D0/
  DATA v(28)/0.877582561890372716116281582D0/
  DATA v(29)/ - 0.841470984807896506652502321D0/
  DATA v(30)/ - 0.479425538604203000273287935D0/
  DATA v(31)/0.0D0/
  DATA v(32)/0.479425538604203000273287935D0/
  DATA v(33)/ - 1.55740772465490223050697485D0/
  DATA v(34)/ - 0.546302489843790513255179465D0/
  DATA v(35)/0.0D0/
  DATA v(36)/0.546302489843790513255179465D0/
  DATA v(37)/2.30258509299404568401799145D0/
  DATA v(38)/2.99573227355399099343522357D0/
  DATA v(39)/3.40119738166215537541323669D0/
  DATA v(40)/3.68887945411393630285245569D0/
  DATA v(41)/1.0D0/
  DATA v(42)/1.30102999566398119521373889D0/
  DATA v(43)/1.4771212547196624372950279D0/
  DATA v(44)/1.60205999132796239042747778D0/
  DATA v(45)/1.00000100530050531421637777D0/
  DATA v(46)/0.999843012323855043126609044D0/
  DATA v(47)/1.00003876575137232151808428D0/
  DATA v(48)/0.992002154326025434343372944D0/
  DATA v(49)/1.54308063481524377847790562D0/
  DATA v(50)/1.12762596520638078522622516D0/
  DATA v(51)/1.0D0/
  DATA v(52)/1.12762596520638078522622516D0/
  DATA v(53)/ - 1.175201193643801456882381851D0/
  DATA v(54)/ - 0.521095305493747361622425626D0/
  DATA v(55)/0.0D0/
  DATA v(56)/0.521095305493747361622425626D0/
  DATA v(57)/ - 0.761594155955764888119458282D0/
  DATA v(58)/ - 0.462117157260009758502318483D0/
  DATA v(59)/0.0D0/
  DATA v(60)/0.462117157260009758592318483D0/
  !
  !***FIRST EXECUTABLE STATEMENT  QCINTD
  !
  IF ( Kprint>=2 ) WRITE (UNIT=Lun,FMT=99001)
  !
  !     Exercise routines in Category C1.
  !
  y(1) = DINT(10.465890D0)
  y(2) = DINT(79.32178D0)
  y(3) = DINT(900.0D0)
  y(4) = DINT(4.0D0)
  !
  !     Exercise routines in Category C2.
  !
  y(5) = DSQRT(1.0D0)
  y(6) = DSQRT(7921.0D0)
  y(7) = DSQRT(900.0D0)
  y(8) = DSQRT(4000D+10)
  !
  !     Exercise routines in Category C4A.
  !
  y(9) = DACOS(-1.0D0)
  y(10) = DACOS(-0.5D0)
  y(11) = DACOS(0.0D0)
  y(12) = DACOS(0.5D0)
  y(13) = DASIN(-1.0D0)
  y(14) = DASIN(-0.5D0)
  y(15) = DASIN(0.0D0)
  y(16) = DASIN(0.5D0)
  y(17) = DATAN(-1.0D0)
  y(18) = DATAN(-0.5D0)
  y(19) = DATAN(0.0D0)
  y(20) = DATAN(0.5D0)
  y(21) = DATAN2(-1.0D0,1.5D0)
  y(22) = DATAN2(-0.5D0,1.0D0)
  y(23) = DATAN2(1.0D0,-0.5D0)
  y(24) = DATAN2(1.5D0,-1.0D0)
  y(25) = DCOS(-1.0D0)
  y(26) = DCOS(-0.5D0)
  y(27) = DCOS(0.0D0)
  y(28) = DCOS(0.5D0)
  y(29) = DSIN(-1.0D0)
  y(30) = DSIN(-0.5D0)
  y(31) = DSIN(0.0D0)
  y(32) = DSIN(0.5D0)
  y(33) = DTAN(-1.0D0)
  y(34) = DTAN(-0.5D0)
  y(35) = DTAN(0.0D0)
  y(36) = DTAN(0.5D0)
  !
  !     Exercise routines in Category C4B.
  !
  y(37) = DLOG(10.0D0)
  y(38) = DLOG(20.0D0)
  y(39) = DLOG(30.0D0)
  y(40) = DLOG(40.0D0)
  y(41) = DLOG10(10.0D0)
  y(42) = DLOG10(20.0D0)
  y(43) = DLOG10(30.0D0)
  y(44) = DLOG10(40.0D0)
  y(45) = DEXP(1.0053D-06)
  y(46) = DEXP(-1.57D-04)
  y(47) = DEXP(3.8765D-05)
  y(48) = DEXP(-8.03D-03)
  !
  !     Exercise routines in Category C4C.
  !
  y(49) = DCOSH(-1.0D0)
  y(50) = DCOSH(-0.5D0)
  y(51) = DCOSH(0.0D0)
  y(52) = DCOSH(0.5D0)
  y(53) = DSINH(-1.0D0)
  y(54) = DSINH(-0.5D0)
  y(55) = DSINH(0.0D0)
  y(56) = DSINH(0.5D0)
  y(57) = DTANH(-1.0D0)
  y(58) = DTANH(-0.5D0)
  y(59) = DTANH(0.0D0)
  y(60) = DTANH(0.5D0)
  !
  !     Check for possible errors.
  !
  Ipass = 1
  errtol = DSQRT(D1MACH(4))
  DO i = 1, 60
    IF ( ABS(v(i)-y(i))>=errtol*ABS(v(i))+errtol ) THEN
      Ipass = 0
      IF ( Kprint>=2 ) WRITE (UNIT=Lun,FMT=99003) i, y(i), v(i)
    ENDIF
  ENDDO
  IF ( Ipass/=0.AND.Kprint>=2 ) WRITE (UNIT=Lun,FMT=99002)
  RETURN
  99001 FORMAT (//' Test of double precision Fullerton intrinsic ','routines')
  99002 FORMAT (' Double precision Fullerton intrinsic function ','routines o.k.')
  99003 FORMAT (' For I  = ',I3,'  test fails with ',/' computed result = ',1P,&
    E38.30,/' and true result = ',E38.30)
END SUBROUTINE QCINTD
!DECK QCINTS
SUBROUTINE QCINTS(Lun,Kprint,Ipass)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  QCINTS
  !***PURPOSE  Quick check for the single precision Fullerton
  !            elementary intrinsic functions.
  !***LIBRARY   FNLIB
  !***CATEGORY  C
  !***TYPE      SINGLE PRECISION (QCINTS-S, QCINTD-D, QCTINC-C)
  !***KEYWORDS  QUICK CHECK
  !***AUTHOR  Boland, W. Robert, (LANL)
  !           Rivera, Shawn M., (LANL)
  !***DESCRIPTION
  !
  !   This subroutine does a quick check for the single precision
  !   Fullerton intrinsic functions.
  !
  !   Parameter list-
  !
  !   LUN      input INTEGER value to designate the external device unit
  !            for message output
  !   KPRINT   input INTEGER value to specify amount of printing to be
  !            done by quick check
  !   IPASS    output INTEGER value indicating whether tests passed or
  !            failed
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  ACOS, ALOG, ALOG10, ASIN, ATAN, ATAN2, CABS, COS,
  !                    COSH, EXP, R1MACH, SIN, SINH, SQRT, TAN, TANH
  !***REVISION HISTORY  (YYMMDD)
  !   900711  DATE WRITTEN
  !***END PROLOGUE  QCINTS
  !     .. Scalar Arguments ..
  INTEGER Ipass, Kprint, Lun
  !     .. Local Scalars ..
  REAL errtol
  INTEGER i
  !     .. Local Arrays ..
  REAL v(60), y(60)
  !     .. External Functions ..
  REAL ACOS, ALOG, ALOG10, ASIN, ATAN, ATAN2, CABS, COS, COSH, &
    EXP, R1MACH, SIN, SINH, SQRT, TAN, TANH
  EXTERNAL ACOS, ALOG, ALOG10, ASIN, ATAN, ATAN2, CABS, COS, COSH, &
    EXP, R1MACH, SIN, SINH, SQRT, TAN, TANH
  !     .. Intrinsic Functions ..
  INTRINSIC ABS, CMPLX
  !
  !     Correct values through different calculations are stored in V(*)
  !
  !     .. Data statements ..
  DATA v(1)/1.0/
  DATA v(2)/89.0/
  DATA v(3)/30.0/
  DATA v(4)/6.324555320337E+06/
  DATA v(5)/10.55327437339/
  DATA v(6)/79.32157587945/
  DATA v(7)/901.0429556913/
  DATA v(8)/4.00000E+13/
  DATA v(9)/3.14159265359/
  DATA v(10)/2.094395102393/
  DATA v(11)/1.570796326795/
  DATA v(12)/1.047197551197/
  DATA v(13)/ - 1.570796326795/
  DATA v(14)/ - 0.5235987755983/
  DATA v(15)/0.0/
  DATA v(16)/0.5235987755983/
  DATA v(17)/ - 0.7853981633974/
  DATA v(18)/ - 0.4636476090008/
  DATA v(19)/0.0/
  DATA v(20)/0.4636476090008/
  DATA v(21)/ - 0.5880026035475/
  DATA v(22)/ - 0.4636476090008/
  DATA v(23)/2.0344438552856/
  DATA v(24)/2.158798930342/
  DATA v(25)/0.5403023058681/
  DATA v(26)/0.8775825618903/
  DATA v(27)/1.0/
  DATA v(28)/0.8775825618903/
  DATA v(29)/ - 0.8414709848079/
  DATA v(30)/ - 0.4794255386042/
  DATA v(31)/0.0/
  DATA v(32)/0.4794255386042/
  DATA v(33)/ - 1.557407724655/
  DATA v(34)/ - 0.5463024898437/
  DATA v(35)/0.0/
  DATA v(36)/0.5463024898437/
  DATA v(37)/2.302585092994/
  DATA v(38)/2.995732273554/
  DATA v(39)/3.401197381662/
  DATA v(40)/3.688879454114/
  DATA v(41)/1.0/
  DATA v(42)/1.301029995664/
  DATA v(43)/1.47712125472/
  DATA v(44)/1.602059991328/
  DATA v(45)/1.000001005301/
  DATA v(46)/0.9998430123238/
  DATA v(47)/1.000038765751/
  DATA v(48)/0.992002154326/
  DATA v(49)/1.543080634815/
  DATA v(50)/1.127625965206/
  DATA v(51)/1.0/
  DATA v(52)/1.127625965206/
  DATA v(53)/ - 1.175201193644/
  DATA v(54)/ - 0.5210953054937/
  DATA v(55)/0.0/
  DATA v(56)/0.5210953054937/
  DATA v(57)/ - 0.7615941559557/
  DATA v(58)/ - 0.46211715726/
  DATA v(59)/0.0/
  DATA v(60)/0.46211715726/
  !
  !***FIRST EXECUTABLE STATEMENT  QCINTS
  !
  IF ( Kprint>=2 ) WRITE (UNIT=Lun,FMT=99001)
  !
  !     Exercise routines in Category C2.
  !
  y(1) = SQRT(1.0)
  y(2) = SQRT(7921.0)
  y(3) = SQRT(900.0)
  y(4) = SQRT(4.00000E+13)
  !
  !     Exercise routines in Category C4.
  !
  y(5) = CABS(CMPLX(10.46,1.4))
  y(6) = CABS(CMPLX(79.32,0.5))
  y(7) = CABS(CMPLX(900.999,8.9))
  y(8) = CABS(CMPLX(4.00000E+13,1.5))
  !
  !     Exercise routines in Category C4A.
  !
  y(9) = ACOS(-1.0)
  y(10) = ACOS(-0.5)
  y(11) = ACOS(0.0)
  y(12) = ACOS(0.5)
  y(13) = ASIN(-1.0)
  y(14) = ASIN(-0.5)
  y(15) = ASIN(0.0)
  y(16) = ASIN(0.5)
  y(17) = ATAN(-1.0)
  y(18) = ATAN(-0.5)
  y(19) = ATAN(0.0)
  y(20) = ATAN(0.5)
  y(21) = ATAN2(-1.0,1.5)
  y(22) = ATAN2(-0.5,1.0)
  y(23) = ATAN2(1.0,-0.5)
  y(24) = ATAN2(1.5,-1.0)
  y(25) = COS(-1.0)
  y(26) = COS(-0.5)
  y(27) = COS(0.0)
  y(28) = COS(0.5)
  y(29) = SIN(-1.0)
  y(30) = SIN(-0.5)
  y(31) = SIN(0.0)
  y(32) = SIN(0.5)
  y(33) = TAN(-1.0)
  y(34) = TAN(-0.5)
  y(35) = TAN(0.0)
  y(36) = TAN(0.5)
  !
  !     Exercise routines in Category C4B.
  !
  y(37) = ALOG(10.0)
  y(38) = ALOG(20.0)
  y(39) = ALOG(30.0)
  y(40) = ALOG(40.0)
  y(41) = ALOG10(10.0)
  y(42) = ALOG10(20.0)
  y(43) = ALOG10(30.0)
  y(44) = ALOG10(40.0)
  y(45) = EXP(1.0053E-06)
  y(46) = EXP(-1.57000E-04)
  y(47) = EXP(3.87650E-05)
  y(48) = EXP(-8.03000E-03)
  !
  !     Exercise routines in Category C4C.
  !
  y(49) = COSH(-1.0)
  y(50) = COSH(-0.5)
  y(51) = COSH(0.0)
  y(52) = COSH(0.5)
  y(53) = SINH(-1.00000)
  y(54) = SINH(-0.50000)
  y(55) = SINH(0.000000)
  y(56) = SINH(0.500000)
  y(57) = TANH(-1.00000)
  y(58) = TANH(-0.50000)
  y(59) = TANH(0.000000)
  y(60) = TANH(0.500000)
  !
  !     Check for possible errors.
  !
  Ipass = 1
  errtol = SQRT(R1MACH(4))
  DO i = 1, 60
    IF ( ABS(v(i)-y(i))>=errtol*ABS(v(i))+errtol ) THEN
      Ipass = 0
      IF ( Kprint>=2 ) WRITE (UNIT=Lun,FMT=99003) i, y(i), v(i)
    ENDIF
  ENDDO
  IF ( Ipass/=0.AND.Kprint>=2 ) WRITE (UNIT=Lun,FMT=99002)
  RETURN
  99001 FORMAT (//' Test of single precision Fullerton intrinsic ','routines')
  99002 FORMAT (' Single precision Fullerton intrinsic function ','routines o.k.')
  99003 FORMAT (' For I  = ',I3,'  test fails with ',/' computed result = ',1P,&
    E22.14,/' and true result = ',E22.14)
END SUBROUTINE QCINTS
