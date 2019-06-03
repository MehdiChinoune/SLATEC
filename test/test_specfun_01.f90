MODULE TEST02_MOD
  USE service, ONLY : SP, DP
  IMPLICIT NONE

CONTAINS
  !** SFNCK
  SUBROUTINE SFNCK(Lun,Kprint,Ipass)
    !>
    !  Quick check for the single precision Fullerton
    !            special functions.
    !***
    ! **Library:**   SLATEC
    !***
    ! **Keywords:**  QUICK CHECK
    !***
    ! **Author:**  Boland, W. Robert, (LANL)
    !           Chow, Jeff, (LANL)
    !***
    ! **Description:**
    !
    !     This subroutine does a quick check for the single precision
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
    ! **Routines called:**  ACOSH, AI, AIE, ALI, ALNREL, ASINH, ATANH, BESI0,
    !                    BESI0E, BESI1, BESI1E, BESJ0, BESJ1, BESK0, BESK0E,
    !                    BESK1, BESK1E, BESKES, BESKS, BESY0, BESY1, BETA,
    !                    BETAI, BI, BIE, BINOM, CBRT, CHU, COSDG, COT, DAWS,
    !                    E1, EI, ERF, EXPREL, FAC, GAMI, GAMIC, GAMIT,
    !                    GAMR, POCH, POCH1, PSI, R1MACH, R9ATN1,
    !                    R9LN2R, SINDG, SPENC

    !* REVISION HISTORY  (YYMMDD)
    !   800901  DATE WRITTEN
    !   891115  REVISION DATE from Version 3.2
    !   891120  Checks of remainder of FNLIB routines added and code
    !           reorganized.  (WRB)
    !   900330  Prologue converted to Version 4.0 format.  (BAB)
    !   900727  Added EXTERNAL statement.  (WRB)
    USE slatec, ONLY : AI, AIE, ALI, ALNREL, BESI0, BESI0E, BESI1, BESI1E, BESK0, &
      BESK0E, BESK1, BESK1E, BESKES, BESKS, BETA, BETAI, BI, BIE, BINOM, CBRT, CHU, &
      COSDG, COT, DAWS, E1, EI, EXPREL, FAC, GAMI, GAMIC, GAMIT, GAMR, POCH, POCH1, &
      PSI, R1MACH, R9ATN1, R9LN2R, SINDG, SPENC
    INTEGER i, Lun, Kprint, Ipass
    REAL(SP) y(105), errmax, errtol, abserr, relerr
    !
    !     Correct values through different calculations are stored in V(*)
    !
    REAL, PARAMETER :: v(87) = [ .834451800000000000000000000000E+09, &
      .225082957512000000000000000000E+13,  .130767436800000000000000000000E+13, &
      .822283865417792281772556288000E+34, -.200000000000000000000000000000E+01, &
      .998340790000000000000000000000E+02,  .866025403784438646763723170753E+00, &
      -.707106781186547524400844362105E+00, .642092615934330703006419986594E+00, &
      -.183048772171245191926801943897E+01,-.290819127993551070285950148310E+00, &
      -.111606410275738687122866817478E+00, .500000000000000000000000000000E+00, &
      .707106781186547524400844362105E+00,  .137149838147233638243285631505E+00, &
      -.100000050000033333358333416027E-05, .100125104231803398984880296644E+01, &
      .995016625083194642609402280122E+00,  .243720864865315055824104923715E+00, &
      .193147180559945309417232121458E+00, -.378671043061087976727207184637E+00, &
      .104516378011749278484458888919E+01,  .559773594776160811746795939295E+00, &
      .100019582406632651901909339800E+00,  .454219904863173579920523812663E+00, &
      .189511781635593675546652093433E+01,  .582240526465012505902656320160E+00, &
      .164493406684822643647241516665E+01, .318309886183790671537767526733E+00, &
      .882395720020380090550940262394E-06, -.282094791773878143474039725759E+00, &
      .187500000000000000000000000000E+01,  .513516668382050295584635612122E-01, &
      .598750000000000000000000000000E+02,  .157079632679489661923132169164E+01, &
      .755006169037464042751871235437E-03,  .422784335098467139393487909918E+00, &
      .230300103429768637527259355045E+01,  .999856618263723706885830759463E+00, &
      .888290707183956735878281870759E+00,  .135335283236612691893999494971E+00, &
      .346930306295801456170933128256E-03,  .786938680574733152792400930048E+00, &
      .631673391775258123291222663623E-01,  .381281566461770916149261183171E+00, &
      .265625000000000000000000000000E+00,  .424436383502022295934042352455E+00, &
      .337000659742093423383019719632E+00,  .227958530233606726743720444020E+01, &
      .272398718236044468945442320700E+02,  .159063685463732906338225442450E+01, &
      .243356421424505271991430504400E+02,  .113893872749533435652719574910E+00, &
      .369109833404259427473526100740E-02,  .139865881816522427284598806997E+00, &
      .404461344545216420836502183700E-02,  .308508322553671039533384319255E+00, &
      .183540812609328353073650751820E+00,  .163972266944542356926122903850E+00, &
      .215269289248937659158505143243E+00,  .841568215070771417919124867127E+00, &
      .547807564313518986868201568700E+00,  .600273858788312582936045656600E+00, &
      .103347684706868857317535710603E+01,  .886226925452758013649083741000E+00, &
      .132934038817913702047362561200E+01,  .288023750772146354435952215970E+01, &
      .560499121639792869931128243359E+00,  .672598945967751443917353892000E+00, &
      .964058489220443736281540578570E+00,  .461068504447894558439575873876E+00, &
      .922137008895789116879151747751E+00,  .231693606480833489769125254500E+00, &
      .157259233804704899952660465400E-01,  .293277159129947362450897433147E+00, &
      .219322205128712060862850888400E+00,  .854277043103155493300048798776E+00, &
      .187894150374789500090933504950E+01,  .674892411115630212865414309867E+00, &
      .464750480196092515019775411670E+00,  .249999999999999999999999999880E+00, &
      .735008609300377745369706799000E+00,  .406961787650672979742685260000E+00, &
      .448256669291582953916931735480E+00,  .596347362323194074341078499290E+00, &
      .757342086122175953454414369190E+00,  .757872156141312106043351240000E+00 ]
    !* FIRST EXECUTABLE STATEMENT  SFNCK
    !
    !     Exercise routines in Category C1.
    !
    y(1) = BINOM(35,12)
    y(2) = BINOM(50,15)
    y(3) = FAC(15)
    y(4) = FAC(31)
    !
    !     Exercise routines in Category C2
    !
    y(5) = CBRT(-8.E0)
    y(6) = CBRT(.995030624365703964475039000000E6)
    !
    !     Exercise routines in Category C4A.
    !
    y(7) = COSDG(30.E0)
    y(8) = COSDG(135.E0)
    y(9) = COT(1.E0)
    y(10) = COT(-.5E0)
    y(11) = R9ATN1(.5E0)
    y(12) = R9ATN1(2.E0)
    y(13) = SINDG(30.E0)
    y(14) = SINDG(135.E0)
    !
    !     Exercise routines in Category C4B.
    !
    y(15) = ALNREL(.147E0)
    y(16) = ALNREL(-.1E-5)
    y(17) = EXPREL(.25E-2)
    y(18) = EXPREL(-.1E-1)
    y(19) = R9LN2R(.5E0)
    y(20) = R9LN2R(1.E0)
    !
    !     Exercise routines in Category C5.
    !
    y(21) = ALI(.5E0)
    y(22) = ALI(2.E0)
    y(23) = E1(.5E0)
    y(24) = E1(1.5E0)
    y(25) = EI(.5E0)
    y(26) = EI(1.E0)
    y(27) = SPENC(.5E0)
    y(28) = SPENC(1.E0)
    y(29) = GAMR(-1.5E0)*GAMR(2.5E0)
    y(30) = GAMR(10.5E0)
    !
    !     Exercise routines in Category C7A.
    !
    y(31) = POCH(-.5E0,1.5E0)
    y(32) = POCH(.5E0,3.E0)
    y(33) = POCH1(.5E0,2.5E0)
    y(34) = POCH1(10.5E0,2.E0)
    !
    !     Exercise routines in Category C7B.
    !
    y(35) = BETA(.5E0,1.5E0)
    y(36) = BETA(5.5E0,5.5E0)
    !
    !     Exercise routines in Category C7C.
    !
    y(37) = PSI(2.E0)
    y(38) = PSI(10.5E0)
    !
    !     Exercise routines in Category C7E.
    !
    y(39) = GAMI(1.E0,8.85E0)
    y(40) = GAMI(2.E0,3.75E0)
    y(41) = GAMIC(1.E0,2.E0)
    y(42) = GAMIC(2.E0,10.4E0)
    y(43) = GAMIT(1.E0,.5E0)
    y(44) = GAMIT(2.E0,3.75E0)
    !
    !     Exercise routines in Category C7F.
    !
    y(45) = BETAI(.5E0,2.E0,1.5E0)
    y(46) = BETAI(.25E0,1.5E0,2.E0)
    !
    !     Exercise routines in Category C8C.
    !
    y(47) = DAWS(.5E0)
    y(48) = DAWS(1.84E0)
    !
    !     Exercise routines in Category C10B1.
    !
    y(49) = BESI0(2.E0)
    y(50) = BESI0(5.E0)
    y(51) = BESI1(2.E0)
    y(52) = BESI1(5.E0)
    y(53) = BESK0(2.E0)
    y(54) = BESK0(5.E0)
    y(55) = BESK1(2.E0)
    y(56) = BESK1(5.E0)
    y(57) = BESI0E(2.E0)
    y(58) = BESI0E(5.E0)
    y(59) = BESI1E(5.E0)
    y(60) = BESI1E(2.E0)
    y(61) = BESK0E(2.E0)
    y(62) = BESK0E(5.E0)
    y(63) = BESK1E(5.E0)
    y(64) = BESK1E(2.E0)
    !
    !     Exercise routines in Category C10B3.
    !
    CALL BESKES(.5E0,2.E0,3,y(65:67))
    CALL BESKES(.5E0,5.E0,3,y(68:70))
    CALL BESKS(.5E0,1.E0,2,y(71:72))
    !
    !     Exercise routines in Category C10D.
    !
    y(73) = AI(.5E0)
    y(74) = AI(2.5E0)
    y(75) = AIE(.5E0)
    y(76) = AIE(2.5E0)
    y(77) = BI(.5E0)
    y(78) = BI(1.5E0)
    y(79) = BIE(.5E0)
    y(80) = BIE(2.5E0)
    !
    !     Exercise routines in Category C11.
    !
    y(81) = CHU(1.E0,2.E0,4.E0)
    y(82) = CHU(5.E0/6.E0,5.E0/3.E0,4.E0/3.E0)
    y(83) = CHU(.75E0,.75E0,2.5E0)
    y(84) = CHU(1.E0,1.E0,1.5E0)
    y(85) = CHU(1.E0,1.E0,1.E0)
    y(86) = CHU(1.E0,1.E0,-LOG(.5E0))
    y(87) = CHU(.5E0,.5E0,1.E0)
    !
    !     Check for possible errors
    !
    errmax = R1MACH(4)
    errtol = SQRT(errmax)
    DO i = 1, 87
      abserr = ABS(v(i)-y(i))
      relerr = abserr/ABS(v(i))
      errmax = MAX(relerr,errmax)
      IF ( relerr>errtol.AND.Kprint>=2 ) WRITE (Lun,99001) i, relerr, abserr
      99001 FORMAT (' For I  = ',I3,'  test fails with RELERR  = ',E38.30,&
        '  and ABSERR  = ',E38.30)
    END DO
    Ipass = 0
    IF ( errmax<=errtol ) Ipass = 1
    IF ( Ipass/=0.AND.Kprint>=2 ) WRITE (Lun,99002)
    !
    99002 FORMAT (' Single precision Fullerton special function ',' routines o.k.')
    RETURN
  END SUBROUTINE SFNCK
END MODULE TEST02_MOD
!** TEST02
PROGRAM TEST02
  USE TEST02_MOD, ONLY : SFNCK
  USE slatec, ONLY : I1MACH, XSETF, XSETUN, XERMAX
  USE common_mod, ONLY : GET_ARGUMENT
  IMPLICIT NONE
  !>
  !  Driver for testing SLATEC subprograms
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  C
  !***
  ! **Type:**      SINGLE PRECISION (TEST02-S, TEST03-D, TEST04-C)
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
  !        single precision Fullerton routines
  !
  !***
  ! **References:**  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***
  ! **Routines called:**  I1MACH, SFNCK, XERMAX, XSETF, XSETUN

  !* REVISION HISTORY  (YYMMDD)
  !   890618  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900524  Cosmetic changes to code.  (WRB)
  INTEGER ipass, kprint, lin, lun, nfail
  !* FIRST EXECUTABLE STATEMENT  TEST02
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
  END IF
  !
  !     Test single precision Fullerton routines
  !
  CALL SFNCK(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF ( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST02 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST02  *************')
  END IF
  STOP
END PROGRAM TEST02
