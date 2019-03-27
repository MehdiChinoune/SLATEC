MODULE TEST04_MOD
  IMPLICIT NONE

CONTAINS
  !** CFNCK
  SUBROUTINE CFNCK(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for the complex Fullerton special functions.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Boland, W. Robert, (LANL)
    !           Chow, Jeff, (LANL)
    !           Rivera, Shawn, (LANL)
    !***
    ! **Description:**
    !
    !     This subroutine does a quick check for the complex
    !     routines in the Fullerton special function library.
    !
    !     Parameter list-
    !
    !     LUN      input integer value to designate the external
    !              device unit for message output
    !     KPRINT   input integer value to specify amount of
    !              printing to be done by quick check
    !     IPASS    output value indicating whether tests passed or
    !              failed
    !
    !***
    ! **Routines called:**  C0LGMC, CACOS, CACOSH, CASIN, CASINH, CATAN,
    !                    CATAN2, CATANH, CBETA, CCBRT, CCOSH, CCOT, CEXPRL,
    !                    CGAMMA, CGAMR, CLBETA, CLNGAM, CLNREL, CLOG10,
    !                    CPSI, CSINH, CTAN, CTANH, R1MACH

    !* REVISION HISTORY  (YYMMDD)
    !   800901  DATE WRITTEN
    !   891115  REVISION DATE from Version 3.2
    !   891120  Checks of remainder of FNLIB routines added and code
    !           reorganized.  (WRB)
    !   900330  Prologue converted to Version 4.0 format.  (BAB)
    !   900727  Added EXTERNAL statement.  (WRB)

    INTEGER i, Lun, Kprint, Ipass
    REAL sqrt2, sqrt3, pi, errmax, errtol, abserr, relerr
    REAL, EXTERNAL :: R1MACH
    COMPLEX w(48)
    COMPLEX, EXTERNAL :: C0LGMC, CACOS, CACOSH, CASIN, CASINH, CATAN, CATAN2, &
      CATANH, CBETA, CCBRT, CCOSH, CCOT, CEXPRL, CGAMMA, CGAMR, CLBETA, CLNGAM, &
      CLNREL, CLOG10, CPSI, CSINH, CTAN, CTANH
    !
    !     Constants to be used
    !
    COMPLEX, PARAMETER :: c1 = (1.E0,0.E0), ci = (0.E0,1.E0)
    DATA sqrt2 /.14142135623730950488E1/
    DATA sqrt3 /.17320508075688772935E1/
    DATA pi /3.14159265358979323846E0/
    !
    !     Complex values through different calculations are stored in C(*)
    !
    COMPLEX, PARAMETER :: c(48) = [ (.121699028117870E1,.326091563038355E0), &
      (.866025403784438E0,.500000000000000E0), &
      (.520802437952465E0,-.196048071390002E1), &
      (.599865470357589E0,.113287925945897E1), &
      (.970930856437313E0,-.113287925945897E1), &
      (.104999388884240E1,.196048071389998E1), &
      (.313314753080534E-1,.541264220944095E-1), &
      (-.785398163397449E0,.658478948462413E0), &
      (-.785398163397449E0,-.658478948462413E0), &
      (.785398163397449E0,-.658478948462413E0), &
      (.313314753080534E-1,.541264220944095E-1), &
      (-.313314753080534E-1,.541264220944095E-1), &
      (.183048772171245E1,.000000000000000E0), &
      (-.757236713834364E-1,-.961745759068982E0), &
      (-.813630257280238E-1,.103336966511721E1), &
      (.546302489843789E0,.000000000000000E0), &
      (.150514997831990E0,-.341094088460459E0), &
      (.301029995663980E0,.227396058973639E0), &
      (.000000000000000E0,.636619772367581E0), &
      (.137802461354738E1,.909330673631480E0), &
      (.303123109082158E-1,-.244978663126864E0), &
      (.693147180559947E0,.523598775598298E0), &
      (-.152857091948100E1,.114371774040242E1), &
      (.144363547517882E1,.157079632679490E1), &
      (-.100000000000000E1,.000000000000000E0), &
      (.181878614736412E1,.586225017697977E0), &
      (.402359478108525E0,.101722196789785E1), &
      (.549306144334055E0,-.157079632679490E1), &
      (.000000000000000E0,-.117520119364380E1), &
      (-.642148124715515E0,-.106860742138277E1), &
      (.397515306849130E0,.104467701612914E1), &
      (-.117520119364380E1,.000000000000000E0), &
      (-.116673625724091E1,-.243458201185722E0), &
      (.761594155955766E0,.000000000000000E0), &
      (.365427607174532E-1,-.612881308922810E-1), &
      (.896860330225849E-2,.244804656578857E-1), &
      (.177245385090552E1,.000000000000000E0), &
      (.300694617260656E0,-.424967879433124E0), &
      (.110951302025214E1,-.156806064476794E1), &
      (.183074439659052E1,.569607641036682E0), &
      (-.340863758923258E1,.142127515954291E1), &
      (-.156059525546301E1,.152533527872833E1), &
      (-.211272372936533E0,-.765528316537801E0), &
      (.380273164249058E-1,-.286343074460341E0), &
      (-.268079774264798E1,.130151697855085E1), &
      (-.164841998888369E1,.785398163397448E0), &
      (-.196351002602143E1,.000000000000000E0), &
      (.161278484461574E1,.147079632679497E1) ]
    !* FIRST EXECUTABLE STATEMENT  CFNCK
    !
    !     Compute functional values
    !
    !     Exercise routines in Category C2.
    !
    w(1) = CCBRT(sqrt2*(1.E0+ci))
    w(2) = CCBRT(ci)
    !
    !     Exercise routines in Category C4A.
    !
    w(3) = CACOS(pi+sqrt3*ci)
    w(4) = CACOS(sqrt2-.25E0*pi*ci)
    w(5) = CASIN(sqrt2-.25E0*pi*ci)
    w(6) = CASIN(pi+sqrt3*ci)
    w(7) = CATAN(.3125E-1+.541265877365273E-1*ci)
    w(8) = CATAN(-.5E0+.866025403784438E0*ci)
    w(9) = CATAN2(-.5E0-.866025403784438E0*ci,c1)
    w(10) = CATAN2(.5E0-.866025403784438E0*ci,c1)
    w(11) = CATAN2(.3125E-1+.541265877365273E-1*ci,c1)
    w(12) = CATAN2(-.3125E-1+.541265877365273E-1*ci,c1)
    w(13) = CCOT(.5E0+0.E0*ci)
    w(14) = CCOT(-1.E0+.5E0*pi*ci)
    w(15) = CTAN(-1.E0+.5E0*pi*ci)
    w(16) = CTAN(.5E0+0.E0*ci)
    !
    !     Exercise routines in Category C4B.
    !
    w(17) = CLOG10(1.E0-ci)
    w(18) = CLOG10(sqrt3+ci)
    w(19) = CEXPRL(pi*ci)
    w(20) = CEXPRL(1.E0+ci)
    w(21) = CLNREL(-.25E0*ci)
    w(22) = CLNREL(sqrt3-1.E0+ci)
    !
    !     Exercise routines in Category C4C.
    !
    w(23) = CACOSH(1.E0-2.E0*ci)
    w(24) = CACOSH(2.E0*ci)
    w(25) = CASINH(-.117520119364380E1+0.E0*ci)
    w(26) = CASINH(2.5E0+1.75E0*ci)
    w(27) = CATANH(1.E0+1.E0*ci)
    w(28) = CATANH(2.E0+0.E0*ci)
    w(29) = CCOSH(1.E0-.5E0*pi*ci)
    w(30) = CCOSH(-1.E0+2.E0*ci)
    w(31) = CSINH(1.E0-1.E0/pi+ci)
    w(32) = CSINH(1.E0+pi*ci)
    w(33) = CTANH(-1.E0+2.E0*ci)
    w(34) = CTANH(1.E0+pi*ci)
    !
    !     Exercise routines in Category C7A.
    !
    w(35) = C0LGMC(.5E0+.5E0*ci)
    w(36) = C0LGMC(1.E0-1.E0*ci)
    w(37) = CGAMMA(.5E0+0.E0*ci)
    w(38) = CGAMMA(.5E0+ci)
    w(39) = CGAMR(.5E0-ci)
    w(40) = CGAMR(1.E0+ci)
    w(41) = CLNGAM(1.1E0+3.2E0*ci)
    w(42) = CLNGAM(1.9E0+2.4E0*ci)
    !
    !     Exercise routines in Category C7B.
    !
    w(43) = CBETA(1.E0+ci,1.E0+ci)
    w(44) = CBETA(2.E0-ci,.5E0+ci)
    w(45) = CLBETA(2.E0+ci,1.E0-2.E0*ci)
    w(46) = CLBETA(1.E0-ci,2.E0+ci)
    !
    !     Exercise routines in Category C7C.
    !
    w(47) = CPSI(.5E0+0.E0*ci)
    w(48) = CPSI(1.E0+5.E0*ci)
    !
    !     Check for possible errors
    !
    errmax = R1MACH(4)
    errtol = SQRT(errmax)
    DO i = 1, 48
      abserr = ABS(c(i)-w(i))
      relerr = abserr/ABS(c(i))
      errmax = MAX(relerr,errmax)
      IF ( relerr>errtol.AND.Kprint>=2 ) WRITE (Lun,99001) i, relerr, abserr
      99001 FORMAT (' For I  = ',I3,'  test fails with RELERR  = ',E38.30,&
        '  and ABSERR  = ',E38.30)
    ENDDO
    Ipass = 0
    IF ( errmax<=errtol ) Ipass = 1
    IF ( Ipass/=0.AND.Kprint>=2 ) WRITE (Lun,99002)
    99002 FORMAT (' Complex Fullerton special function routines o.k.')
    RETURN
  END SUBROUTINE CFNCK
END MODULE TEST04_MOD
!** TEST04
PROGRAM TEST04
  USE TEST04_MOD
  IMPLICIT NONE
  !>
  !***
  !  Driver for testing SLATEC subprograms
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  C
  !***
  ! **Type:**      COMPLEX (TEST02-S, TEST03-D, TEST04-C)
  !***
  ! **Keywords:**  QUICK CHECK DRIVER
  !***
  ! **Author:**  SLATEC Common Mathematical Library Committee
  !***
  ! **Description:**
  !
  !- Usage:
  !     One input data record is required
  !         READ (LIN, '(I1)') KPRINT
  !
  !- Arguments:
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
  !- Description:
  !     Driver for testing SLATEC subprograms
  !        complex Fullerton routines
  !
  !***
  ! **References:**  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***
  ! **Routines called:**  CFNCK, I1MACH, XERMAX, XSETF, XSETUN

  !* REVISION HISTORY  (YYMMDD)
  !   890618  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900524  Cosmetic changes to code.  (WRB)

  INTEGER I1MACH
  INTEGER ipass, kprint, lin, lun, nfail
  !* FIRST EXECUTABLE STATEMENT  TEST04
  lun = I1MACH(2)
  lin = I1MACH(1)
  nfail = 0
  !
  !     Read KPRINT parameter
  !
  CALL GET_ARGUMENT(kprint)
  CALL XERMAX(1000)
  CALL XSETUN(lun)
  IF ( kprint<=1 ) THEN
    CALL XSETF(0)
  ELSE
    CALL XSETF(1)
  ENDIF
  !
  !     Test complex Fullerton routines
  !
  CALL CFNCK(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF ( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST04 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST04  *************')
  ENDIF
  STOP
END PROGRAM TEST04
