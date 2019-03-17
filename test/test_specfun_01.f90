MODULE TEST02_MOD
  IMPLICIT NONE

CONTAINS
  !DECK SFNCK
  SUBROUTINE SFNCK(Lun,Kprint,Ipass)
    IMPLICIT NONE
    REAL ACOSH, AI, AIE, ALI, ALNREL, ASINH, ATANH, BESI0, BESI0E
    !***BEGIN PROLOGUE  SFNCK
    !***PURPOSE  Quick check for the single precision Fullerton
    !            special functions.
    !***LIBRARY   SLATEC
    !***KEYWORDS  QUICK CHECK
    !***AUTHOR  Boland, W. Robert, (LANL)
    !           Chow, Jeff, (LANL)
    !***DESCRIPTION
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
    !***ROUTINES CALLED  ACOSH, AI, AIE, ALI, ALNREL, ASINH, ATANH, BESI0,
    !                    BESI0E, BESI1, BESI1E, BESJ0, BESJ1, BESK0, BESK0E,
    !                    BESK1, BESK1E, BESKES, BESKS, BESY0, BESY1, BETA,
    !                    BETAI, BI, BIE, BINOM, CBRT, CHU, COSDG, COT, DAWS,
    !                    E1, EI, ERF, EXPREL, FAC, GAMI, GAMIC, GAMIT,
    !                    GAMMA, GAMR, POCH, POCH1, PSI, R1MACH, R9ATN1,
    !                    R9LN2R, SINDG, SPENC
    !***REVISION HISTORY  (YYMMDD)
    !   800901  DATE WRITTEN
    !   891115  REVISION DATE from Version 3.2
    !   891120  Checks of remainder of FNLIB routines added and code
    !           reorganized.  (WRB)
    !   900330  Prologue converted to Version 4.0 format.  (BAB)
    !   900727  Added EXTERNAL statement.  (WRB)
    !***END PROLOGUE  SFNCK
    INTEGER i, Lun, Kprint, Ipass
    REAL R1MACH, y(105), v(105), errmax, errtol, abserr, relerr, &
      BESI1, BESI1E, BESJ0, BESJ1, BESK0, BESK0E, BESK1, BESK1E, &
      BESY0, BESY1, BETA, BETAI, BI, BIE, BINOM, CBRT, CHU, &
      COSDG, COT, DAWS, E1, EI, ERF, EXPREL, FAC, GAMI, GAMIC, &
      GAMIT, GAMMA, GAMR, POCH, POCH1, PSI, R9ATN1, R9LN2R, SINDG, &
      SPENC
    EXTERNAL COT, ERF, GAMMA
    !
    !     Correct values through different calculations are stored in V(*)
    !
    DATA v(1)/.834451800000000000000000000000E+09/
    DATA v(2)/.225082957512000000000000000000E+13/
    DATA v(3)/.130767436800000000000000000000E+13/
    DATA v(4)/.822283865417792281772556288000E+34/
    DATA v(5)/ - .200000000000000000000000000000E+01/
    DATA v(6)/.998340790000000000000000000000E+02/
    DATA v(7)/.866025403784438646763723170753E+00/
    DATA v(8)/ - .707106781186547524400844362105E+00/
    DATA v(9)/.642092615934330703006419986594E+00/
    DATA v(10)/ - .183048772171245191926801943897E+01/
    DATA v(11)/ - .290819127993551070285950148310E+00/
    DATA v(12)/ - .111606410275738687122866817478E+00/
    DATA v(13)/.500000000000000000000000000000E+00/
    DATA v(14)/.707106781186547524400844362105E+00/
    DATA v(15)/.137149838147233638243285631505E+00/
    DATA v(16)/ - .100000050000033333358333416027E-05/
    DATA v(17)/.100125104231803398984880296644E+01/
    DATA v(18)/.995016625083194642609402280122E+00/
    DATA v(19)/.243720864865315055824104923715E+00/
    DATA v(20)/.193147180559945309417232121458E+00/
    DATA v(21)/.111112222233333444440000000000E+00/
    DATA v(22)/.314159265359000000000000000000E+01/
    DATA v(23)/.998340790000000000000000000000E-01/
    DATA v(24)/ - .119476321700000000000000000000E+01/
    DATA v(25)/ - .111112222233333444440000000000E+00/
    DATA v(26)/.264665241200000000000000000000E+01/
    DATA v(27)/ - .378671043061087976727207184637E+00/
    DATA v(28)/.104516378011749278484458888919E+01/
    DATA v(29)/.559773594776160811746795939295E+00/
    DATA v(30)/.100019582406632651901909339800E+00/
    DATA v(31)/.454219904863173579920523812663E+00/
    DATA v(32)/.189511781635593675546652093433E+01/
    DATA v(33)/.582240526465012505902656320160E+00/
    DATA v(34)/.164493406684822643647241516665E+01/
    DATA v(35)/.886226925452758013649083741687E+00/
    DATA v(36)/ - .314159265358979323846264338328E+01/
    DATA v(37)/.318309886183790671537767526733E+00/
    DATA v(38)/.882395720020380090550940262394E-06/
    DATA v(39)/ - .282094791773878143474039725759E+00/
    DATA v(40)/.187500000000000000000000000000E+01/
    DATA v(41)/.513516668382050295584635612122E-01/
    DATA v(42)/.598750000000000000000000000000E+02/
    DATA v(43)/.157079632679489661923132169164E+01/
    DATA v(44)/.755006169037464042751871235437E-03/
    DATA v(45)/.422784335098467139393487909918E+00/
    DATA v(46)/.230300103429768637527259355045E+01/
    DATA v(47)/.999856618263723706885830759463E+00/
    DATA v(48)/.888290707183956735878281870759E+00/
    DATA v(49)/.135335283236612691893999494971E+00/
    DATA v(50)/.346930306295801456170933128256E-03/
    DATA v(51)/.786938680574733152792400930048E+00/
    DATA v(52)/.631673391775258123291222663623E-01/
    DATA v(53)/.381281566461770916149261183171E+00/
    DATA v(54)/.265625000000000000000000000000E+00/
    DATA v(55)/.520499877813046537682746653770E+00/
    DATA v(56)/.888388231701707764069578446749E+00/
    DATA v(57)/.424436383502022295934042352455E+00/
    DATA v(58)/.337000659742093423383019719632E+00/
    DATA v(59)/ - .177596771314338304347397013056E+00/
    DATA v(60)/.223890779141235668051827454628E+00/
    DATA v(61)/ - .327579137591465222037734321812E+00/
    DATA v(62)/.576724807756873387202448242187E+00/
    DATA v(63)/.510375672649745119596606592612E+00/
    DATA v(64)/ - .308517625249033780073648984210E+00/
    DATA v(65)/.147863143391226844801050675510E+00/
    DATA v(66)/ - .107032431540937546888370772230E+00/
    DATA v(67)/.227958530233606726743720444020E+01/
    DATA v(68)/.272398718236044468945442320700E+02/
    DATA v(69)/.159063685463732906338225442450E+01/
    DATA v(70)/.243356421424505271991430504400E+02/
    DATA v(71)/.113893872749533435652719574910E+00/
    DATA v(72)/.369109833404259427473526100740E-02/
    DATA v(73)/.139865881816522427284598806997E+00/
    DATA v(74)/.404461344545216420836502183700E-02/
    DATA v(75)/.308508322553671039533384319255E+00/
    DATA v(76)/.183540812609328353073650751820E+00/
    DATA v(77)/.163972266944542356926122903850E+00/
    DATA v(78)/.215269289248937659158505143243E+00/
    DATA v(79)/.841568215070771417919124867127E+00/
    DATA v(80)/.547807564313518986868201568700E+00/
    DATA v(81)/.600273858788312582936045656600E+00/
    DATA v(82)/.103347684706868857317535710603E+01/
    DATA v(83)/.886226925452758013649083741000E+00/
    DATA v(84)/.132934038817913702047362561200E+01/
    DATA v(85)/.288023750772146354435952215970E+01/
    DATA v(86)/.560499121639792869931128243359E+00/
    DATA v(87)/.672598945967751443917353892000E+00/
    DATA v(88)/.964058489220443736281540578570E+00/
    DATA v(89)/.461068504447894558439575873876E+00/
    DATA v(90)/.922137008895789116879151747751E+00/
    DATA v(91)/.231693606480833489769125254500E+00/
    DATA v(92)/.157259233804704899952660465400E-01/
    DATA v(93)/.293277159129947362450897433147E+00/
    DATA v(94)/.219322205128712060862850888400E+00/
    DATA v(95)/.854277043103155493300048798776E+00/
    DATA v(96)/.187894150374789500090933504950E+01/
    DATA v(97)/.674892411115630212865414309867E+00/
    DATA v(98)/.464750480196092515019775411670E+00/
    DATA v(99)/.249999999999999999999999999880E+00/
    DATA v(100)/.735008609300377745369706799000E+00/
    DATA v(101)/.406961787650672979742685260000E+00/
    DATA v(102)/.448256669291582953916931735480E+00/
    DATA v(103)/.596347362323194074341078499290E+00/
    DATA v(104)/.757342086122175953454414369190E+00/
    DATA v(105)/.757872156141312106043351240000E+00/
    !***FIRST EXECUTABLE STATEMENT  SFNCK
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
    !     Exercise routines in Category C4C.
    !
    y(21) = ACOSH(.100617931649094823747218929626E1)
    y(22) = ACOSH(.115919532755239084628557897777E2)
    y(23) = ASINH(.100000000101295145211538706587E0)
    y(24) = ASINH(-.149999999948240634124264852207E1)
    y(25) = ATANH(-.110657208041383998066515207788E0)
    y(26) = ATANH(.989999999992791300663084082410E0)
    !
    !     Exercise routines in Category C5.
    !
    y(27) = ALI(.5E0)
    y(28) = ALI(2.E0)
    y(29) = E1(.5E0)
    y(30) = E1(1.5E0)
    y(31) = EI(.5E0)
    y(32) = EI(1.E0)
    y(33) = SPENC(.5E0)
    y(34) = SPENC(1.E0)
    y(35) = GAMMA(1.5E0)
    y(36) = GAMMA(-.5E0)*GAMMA(1.5E0)
    y(37) = GAMR(-1.5E0)*GAMR(2.5E0)
    y(38) = GAMR(10.5E0)
    !
    !     Exercise routines in Category C7A.
    !
    y(39) = POCH(-.5E0,1.5E0)
    y(40) = POCH(.5E0,3.E0)
    y(41) = POCH1(.5E0,2.5E0)
    y(42) = POCH1(10.5E0,2.E0)
    !
    !     Exercise routines in Category C7B.
    !
    y(43) = BETA(.5E0,1.5E0)
    y(44) = BETA(5.5E0,5.5E0)
    !
    !     Exercise routines in Category C7C.
    !
    y(45) = PSI(2.E0)
    y(46) = PSI(10.5E0)
    !
    !     Exercise routines in Category C7E.
    !
    y(47) = GAMI(1.E0,8.85E0)
    y(48) = GAMI(2.E0,3.75E0)
    y(49) = GAMIC(1.E0,2.E0)
    y(50) = GAMIC(2.E0,10.4E0)
    y(51) = GAMIT(1.E0,.5E0)
    y(52) = GAMIT(2.E0,3.75E0)
    !
    !     Exercise routines in Category C7F.
    !
    y(53) = BETAI(.5E0,2.E0,1.5E0)
    y(54) = BETAI(.25E0,1.5E0,2.E0)
    !
    !     Exercise routines in Category C8A.
    !
    y(55) = ERF(.5E0)
    y(56) = ERF(1.125E0)
    !
    !     Exercise routines in Category C8C.
    !
    y(57) = DAWS(.5E0)
    y(58) = DAWS(1.84E0)
    !
    !     Exercise routines in Category C10A1.
    !
    y(59) = BESJ0(5.E0)
    y(60) = BESJ0(2.E0)
    y(61) = BESJ1(5.E0)
    y(62) = BESJ1(2.E0)
    y(63) = BESY0(2.E0)
    y(64) = BESY0(5.E0)
    y(65) = BESY1(5.E0)
    y(66) = BESY1(2.E0)
    !
    !     Exercise routines in Category C10B1.
    !
    y(67) = BESI0(2.E0)
    y(68) = BESI0(5.E0)
    y(69) = BESI1(2.E0)
    y(70) = BESI1(5.E0)
    y(71) = BESK0(2.E0)
    y(72) = BESK0(5.E0)
    y(73) = BESK1(2.E0)
    y(74) = BESK1(5.E0)
    y(75) = BESI0E(2.E0)
    y(76) = BESI0E(5.E0)
    y(77) = BESI1E(5.E0)
    y(78) = BESI1E(2.E0)
    y(79) = BESK0E(2.E0)
    y(80) = BESK0E(5.E0)
    y(81) = BESK1E(5.E0)
    y(82) = BESK1E(2.E0)
    !
    !     Exercise routines in Category C10B3.
    !
    CALL BESKES(.5E0,2.E0,3,y(83))
    CALL BESKES(.5E0,5.E0,3,y(86))
    CALL BESKS(.5E0,1.E0,2,y(89))
    !
    !     Exercise routines in Category C10D.
    !
    y(91) = AI(.5E0)
    y(92) = AI(2.5E0)
    y(93) = AIE(.5E0)
    y(94) = AIE(2.5E0)
    y(95) = BI(.5E0)
    y(96) = BI(1.5E0)
    y(97) = BIE(.5E0)
    y(98) = BIE(2.5E0)
    !
    !     Exercise routines in Category C11.
    !
    y(99) = CHU(1.E0,2.E0,4.E0)
    y(100) = CHU(5.E0/6.E0,5.E0/3.E0,4.E0/3.E0)
    y(101) = CHU(.75E0,.75E0,2.5E0)
    y(102) = CHU(1.E0,1.E0,1.5E0)
    y(103) = CHU(1.E0,1.E0,1.E0)
    y(104) = CHU(1.E0,1.E0,-LOG(.5E0))
    y(105) = CHU(.5E0,.5E0,1.E0)
    !
    !     Check for possible errors
    !
    errmax = R1MACH(4)
    errtol = SQRT(errmax)
    DO i = 1, 105
      abserr = ABS(v(i)-y(i))
      relerr = abserr/ABS(v(i))
      errmax = MAX(relerr,errmax)
      IF ( relerr>errtol.AND.Kprint>=2 ) WRITE (Lun,99001) i, relerr, abserr
      99001 FORMAT (' For I  = ',I3,'  test fails with RELERR  = ',E38.30,&
        '  and ABSERR  = ',E38.30)
    ENDDO
    Ipass = 0
    IF ( errmax<=errtol ) Ipass = 1
    IF ( Ipass/=0.AND.Kprint>=2 ) WRITE (Lun,99002)
    !
    99002 FORMAT (' Single precision Fullerton special function ',' routines o.k.')
    RETURN
  END SUBROUTINE SFNCK
END MODULE TEST02_MOD
!DECK TEST02
PROGRAM TEST02
  USE TEST02_MOD
  IMPLICIT NONE
  !***BEGIN PROLOGUE  TEST02
  !***PURPOSE  Driver for testing SLATEC subprograms
  !***LIBRARY   SLATEC
  !***CATEGORY  C
  !***TYPE      SINGLE PRECISION (TEST02-S, TEST03-D, TEST04-C)
  !***KEYWORDS  QUICK CHECK DRIVER
  !***AUTHOR  SLATEC Common Mathematical Library Committee
  !***DESCRIPTION
  !
  ! *Usage:
  !     One input data record is required
  !         READ (LIN, '(I1)') KPRINT
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
  !     Driver for testing SLATEC subprograms
  !        single precision Fullerton routines
  !
  !***REFERENCES  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***ROUTINES CALLED  I1MACH, SFNCK, XERMAX, XSETF, XSETUN
  !***REVISION HISTORY  (YYMMDD)
  !   890618  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900524  Cosmetic changes to code.  (WRB)
  !***END PROLOGUE  TEST02
  INTEGER I1MACH
  INTEGER ipass, kprint, lin, lun, nfail
  !***FIRST EXECUTABLE STATEMENT  TEST02
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
  ENDIF
  STOP
END PROGRAM TEST02