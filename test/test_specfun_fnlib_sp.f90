MODULE TEST02_MOD
  USE service, ONLY : SP, DP
  IMPLICIT NONE

CONTAINS
  !** SFNCK
  SUBROUTINE SFNCK(Lun,Kprint,Ipass)
    !> Quick check for the single precision Fullerton
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
    !   891120  Checks of remainder of FNLIB routines added and code reorganized.  (WRB)
    !   900330  Prologue converted to Version 4.0 format.  (BAB)
    !   900727  Added EXTERNAL statement.  (WRB)
    USE service, ONLY : eps_sp
    USE special_functions, ONLY : AI, AIE, ALI, ALNREL, BESI0, BESI0E, BESI1, BESI1E, BESK0, &
      BESK0E, BESK1, BESK1E, BESKES, BESKS, BETA, BETAI, BI, BIE, BINOM, CBRT, CHU, &
      COSDG, COT, DAWS, E1, EI, EXPREL, FAC, GAMI, GAMIC, GAMIT, GAMR, POCH, POCH1, &
      PSI, R9ATN1, R9LN2R, SINDG, SPENC
    !
    INTEGER :: i, Lun, Kprint, Ipass
    REAL(SP) :: y(105), errmax, errtol, abserr, relerr
    !
    !     Correct values through different calculations are stored in V(*)
    !
    REAL(SP), PARAMETER :: v(87) = [ .834451800000000000000000000000E+09_SP, &
      .225082957512000000000000000000E+13_SP,  .130767436800000000000000000000E+13_SP, &
      .822283865417792281772556288000E+34_SP, -.200000000000000000000000000000E+01_SP, &
      .998340790000000000000000000000E+02_SP,  .866025403784438646763723170753E+00_SP, &
      -.707106781186547524400844362105E+00_SP, .642092615934330703006419986594E+00_SP, &
      -.183048772171245191926801943897E+01_SP,-.290819127993551070285950148310E+00_SP, &
      -.111606410275738687122866817478E+00_SP, .500000000000000000000000000000E+00_SP, &
      .707106781186547524400844362105E+00_SP,  .137149838147233638243285631505E+00_SP, &
      -.100000050000033333358333416027E-05_SP, .100125104231803398984880296644E+01_SP, &
      .995016625083194642609402280122E+00_SP,  .243720864865315055824104923715E+00_SP, &
      .193147180559945309417232121458E+00_SP, -.378671043061087976727207184637E+00_SP, &
      .104516378011749278484458888919E+01_SP,  .559773594776160811746795939295E+00_SP, &
      .100019582406632651901909339800E+00_SP,  .454219904863173579920523812663E+00_SP, &
      .189511781635593675546652093433E+01_SP,  .582240526465012505902656320160E+00_SP, &
      .164493406684822643647241516665E+01_SP, .318309886183790671537767526733E+00_SP, &
      .882395720020380090550940262394E-06_SP, -.282094791773878143474039725759E+00_SP, &
      .187500000000000000000000000000E+01_SP,  .513516668382050295584635612122E-01_SP, &
      .598750000000000000000000000000E+02_SP,  .157079632679489661923132169164E+01_SP, &
      .755006169037464042751871235437E-03_SP,  .422784335098467139393487909918E+00_SP, &
      .230300103429768637527259355045E+01_SP,  .999856618263723706885830759463E+00_SP, &
      .888290707183956735878281870759E+00_SP,  .135335283236612691893999494971E+00_SP, &
      .346930306295801456170933128256E-03_SP,  .786938680574733152792400930048E+00_SP, &
      .631673391775258123291222663623E-01_SP,  .381281566461770916149261183171E+00_SP, &
      .265625000000000000000000000000E+00_SP,  .424436383502022295934042352455E+00_SP, &
      .337000659742093423383019719632E+00_SP,  .227958530233606726743720444020E+01_SP, &
      .272398718236044468945442320700E+02_SP,  .159063685463732906338225442450E+01_SP, &
      .243356421424505271991430504400E+02_SP,  .113893872749533435652719574910E+00_SP, &
      .369109833404259427473526100740E-02_SP,  .139865881816522427284598806997E+00_SP, &
      .404461344545216420836502183700E-02_SP,  .308508322553671039533384319255E+00_SP, &
      .183540812609328353073650751820E+00_SP,  .163972266944542356926122903850E+00_SP, &
      .215269289248937659158505143243E+00_SP,  .841568215070771417919124867127E+00_SP, &
      .547807564313518986868201568700E+00_SP,  .600273858788312582936045656600E+00_SP, &
      .103347684706868857317535710603E+01_SP,  .886226925452758013649083741000E+00_SP, &
      .132934038817913702047362561200E+01_SP,  .288023750772146354435952215970E+01_SP, &
      .560499121639792869931128243359E+00_SP,  .672598945967751443917353892000E+00_SP, &
      .964058489220443736281540578570E+00_SP,  .461068504447894558439575873876E+00_SP, &
      .922137008895789116879151747751E+00_SP,  .231693606480833489769125254500E+00_SP, &
      .157259233804704899952660465400E-01_SP,  .293277159129947362450897433147E+00_SP, &
      .219322205128712060862850888400E+00_SP,  .854277043103155493300048798776E+00_SP, &
      .187894150374789500090933504950E+01_SP,  .674892411115630212865414309867E+00_SP, &
      .464750480196092515019775411670E+00_SP,  .249999999999999999999999999880E+00_SP, &
      .735008609300377745369706799000E+00_SP,  .406961787650672979742685260000E+00_SP, &
      .448256669291582953916931735480E+00_SP,  .596347362323194074341078499290E+00_SP, &
      .757342086122175953454414369190E+00_SP,  .757872156141312106043351240000E+00_SP ]
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
    y(5) = CBRT(-8._SP)
    y(6) = CBRT(.995030624365703964475039000000E6_SP)
    !
    !     Exercise routines in Category C4A.
    !
    y(7) = COSDG(30._SP)
    y(8) = COSDG(135._SP)
    y(9) = COT(1._SP)
    y(10) = COT(-.5_SP)
    y(11) = R9ATN1(.5_SP)
    y(12) = R9ATN1(2._SP)
    y(13) = SINDG(30._SP)
    y(14) = SINDG(135._SP)
    !
    !     Exercise routines in Category C4B.
    !
    y(15) = ALNREL(.147_SP)
    y(16) = ALNREL(-.1E-5_SP)
    y(17) = EXPREL(.25E-2_SP)
    y(18) = EXPREL(-.1E-1_SP)
    y(19) = R9LN2R(.5_SP)
    y(20) = R9LN2R(1._SP)
    !
    !     Exercise routines in Category C5.
    !
    y(21) = ALI(.5_SP)
    y(22) = ALI(2._SP)
    y(23) = E1(.5_SP)
    y(24) = E1(1.5_SP)
    y(25) = EI(.5_SP)
    y(26) = EI(1._SP)
    y(27) = SPENC(.5_SP)
    y(28) = SPENC(1._SP)
    y(29) = GAMR(-1.5_SP)*GAMR(2.5_SP)
    y(30) = GAMR(10.5_SP)
    !
    !     Exercise routines in Category C7A.
    !
    y(31) = POCH(-.5_SP,1.5_SP)
    y(32) = POCH(.5_SP,3._SP)
    y(33) = POCH1(.5_SP,2.5_SP)
    y(34) = POCH1(10.5_SP,2._SP)
    !
    !     Exercise routines in Category C7B.
    !
    y(35) = BETA(.5_SP,1.5_SP)
    y(36) = BETA(5.5_SP,5.5_SP)
    !
    !     Exercise routines in Category C7C.
    !
    y(37) = PSI(2._SP)
    y(38) = PSI(10.5_SP)
    !
    !     Exercise routines in Category C7E.
    !
    y(39) = GAMI(1._SP,8.85_SP)
    y(40) = GAMI(2._SP,3.75_SP)
    y(41) = GAMIC(1._SP,2._SP)
    y(42) = GAMIC(2._SP,10.4_SP)
    y(43) = GAMIT(1._SP,.5_SP)
    y(44) = GAMIT(2._SP,3.75_SP)
    !
    !     Exercise routines in Category C7F.
    !
    y(45) = BETAI(.5_SP,2._SP,1.5_SP)
    y(46) = BETAI(.25_SP,1.5_SP,2._SP)
    !
    !     Exercise routines in Category C8C.
    !
    y(47) = DAWS(.5_SP)
    y(48) = DAWS(1.84_SP)
    !
    !     Exercise routines in Category C10B1.
    !
    y(49) = BESI0(2._SP)
    y(50) = BESI0(5._SP)
    y(51) = BESI1(2._SP)
    y(52) = BESI1(5._SP)
    y(53) = BESK0(2._SP)
    y(54) = BESK0(5._SP)
    y(55) = BESK1(2._SP)
    y(56) = BESK1(5._SP)
    y(57) = BESI0E(2._SP)
    y(58) = BESI0E(5._SP)
    y(59) = BESI1E(5._SP)
    y(60) = BESI1E(2._SP)
    y(61) = BESK0E(2._SP)
    y(62) = BESK0E(5._SP)
    y(63) = BESK1E(5._SP)
    y(64) = BESK1E(2._SP)
    !
    !     Exercise routines in Category C10B3.
    !
    CALL BESKES(.5_SP,2._SP,3,y(65:67))
    CALL BESKES(.5_SP,5._SP,3,y(68:70))
    CALL BESKS(.5_SP,1._SP,2,y(71:72))
    !
    !     Exercise routines in Category C10D.
    !
    y(73) = AI(.5_SP)
    y(74) = AI(2.5_SP)
    y(75) = AIE(.5_SP)
    y(76) = AIE(2.5_SP)
    y(77) = BI(.5_SP)
    y(78) = BI(1.5_SP)
    y(79) = BIE(.5_SP)
    y(80) = BIE(2.5_SP)
    !
    !     Exercise routines in Category C11.
    !
    y(81) = CHU(1._SP,2._SP,4._SP)
    y(82) = CHU(5._SP/6._SP,5._SP/3._SP,4._SP/3._SP)
    y(83) = CHU(.75_SP,.75_SP,2.5_SP)
    y(84) = CHU(1._SP,1._SP,1.5_SP)
    y(85) = CHU(1._SP,1._SP,1._SP)
    y(86) = CHU(1._SP,1._SP,-LOG(.5_SP))
    y(87) = CHU(.5_SP,.5_SP,1._SP)
    !
    !     Check for possible errors
    !
    errmax = eps_sp
    errtol = SQRT(errmax)
    DO i = 1, 87
      abserr = ABS(v(i)-y(i))
      relerr = abserr/ABS(v(i))
      errmax = MAX(relerr,errmax)
      IF( relerr>errtol .AND. Kprint>=2 ) WRITE (Lun,99001) i, relerr, abserr
      99001 FORMAT (' For I  = ',I3,'  test fails with RELERR  = ',E38.30,&
        '  and ABSERR  = ',E38.30)
    END DO
    Ipass = 0
    IF( errmax<=errtol ) Ipass = 1
    IF( Ipass/=0 .AND. Kprint>=2 ) WRITE (Lun,99002)
    !
    99002 FORMAT (' Single precision Fullerton special function ',' routines o.k.')
    RETURN
  END SUBROUTINE SFNCK
END MODULE TEST02_MOD
!** TEST02
PROGRAM TEST02
  USE TEST02_MOD, ONLY : SFNCK
  USE ISO_FORTRAN_ENV, ONLY : INPUT_UNIT, OUTPUT_UNIT
  USE common_mod, ONLY : GET_ARGUMENT
  IMPLICIT NONE
  !> Driver for testing SLATEC subprograms
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
  ! **Routines called:**  I1MACH, SFNCK, XERMAX, XSETF

  !* REVISION HISTORY  (YYMMDD)
  !   890618  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900524  Cosmetic changes to code.  (WRB)
  INTEGER :: ipass, kprint, lin, lun, nfail
  !* FIRST EXECUTABLE STATEMENT  TEST02
  lun = OUTPUT_UNIT
  lin = INPUT_UNIT
  nfail = 0
  !
  !     Read KPRINT parameter
  !
  CALL GET_ARGUMENT(kprint)
  !
  !     Test single precision Fullerton routines
  !
  CALL SFNCK(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST02 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST02  *************')
  END IF
  STOP
END PROGRAM TEST02
