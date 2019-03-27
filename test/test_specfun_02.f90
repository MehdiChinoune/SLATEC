MODULE TEST03_MOD
  IMPLICIT NONE

CONTAINS
  !** DFNCK
  SUBROUTINE DFNCK(Lun,Kprint,Ipass)
    IMPLICIT NONE
    !>
    !***
    !  Quick check for the double precision Fullerton
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
    !     This subroutine does a quick check for the double precision
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
    ! **Routines called:**  D1MACH, D9ATN1, D9LN2R, DACOSH, DAI, DAIE, DASINH,
    !                    DATANH, DBESI0, DBESI1, DBESJ0, DBESJ1, DBESK0,
    !                    DBESK1, DBESKS, DBESY0, DBESY1, DBETA, DBETAI, DBI,
    !                    DBIE, DBINOM, DBSI0E, DBSI1E, DBSK0E, DBSK1E,
    !                    DBSKES, DCBRT, DCHU, DCOSDG, DCOT, DDAWS, DE1, DEI,
    !                    DERF, DEXPRL, DFAC, DGAMI, DGAMIC, DGAMIT, DGAMMA,
    !                    DGAMR, DLI, DLNREL, DPOCH, DPOCH1, DPSI, DSINDG,
    !                    DSPENC

    !* REVISION HISTORY  (YYMMDD)
    !   800801  DATE WRITTEN
    !   891115  REVISION DATE from Version 3.2
    !   891120  Checks of remainder of FNLIB routines added and code
    !           reorganized.  (WRB)
    !   900330  Prologue converted to Version 4.0 format.  (BAB)
    !   900727  Added EXTERNAL statement.  (WRB)

    INTEGER i, Lun, Kprint, Ipass
    REAL(8) :: y(105), errmax, errtol, abserr, relerr
    REAL(8), EXTERNAL :: D1MACH, D9ATN1, D9LN2R, DACOSH, DAI, DAIE, DASINH, &
      DATANH, DBESI0, DBESI1, DBESJ0, DBESJ1, DBESK0, DBESK1, DBESY0, DBESY1, &
      DBETA, DBETAI, DBI, DBIE, DBINOM, DBSI0E, DBSI1E, DBSK0E, DBSK1E, DCBRT, &
      DCHU, DCOSDG, DCOT, DDAWS, DE1, DEI, DERF, DEXPRL, DFAC, DGAMI, DGAMIC, &
      DGAMIT, DGAMMA, DGAMR, DLI, DLNREL, DPOCH, DPOCH1, DPSI, DSINDG, DSPENC
    EXTERNAL :: DBESKS, DBSKES
    !
    !     Correct values through different calculations are stored in V(*)
    !
    REAL(8), PARAMETER :: v(105) = [ .834451800000000000000000000000D+09, &
      .225082957512000000000000000000D+13,  .130767436800000000000000000000D+13, &
      .822283865417792281772556288000D+34, -.200000000000000000000000000000D+01, &
      .998340790000000000000000000000D+02,  .866025403784438646763723170753D+00, &
      -.707106781186547524400844362105D+00, .642092615934330703006419986594D+00, &
      -.183048772171245191926801943897D+01,-.290819127993551070285950148310D+00, &
      -.111606410275738687122866817478D+00, .500000000000000000000000000000D+00, &
      .707106781186547524400844362105D+00,  .137149838147233638243285631505D+00, &
      -.100000050000033333358333416027D-05, .100125104231803398984880296644D+01, &
      .995016625083194642609402280122D+00,  .243720864865315055824104923715D+00, &
      .193147180559945309417232121458D+00,  .111112222233333444440000000000D+00, &
      .314159265359000000000000000000D+01,  .998340790000000000000000000000D-01, &
      -.119476321700000000000000000000D+01,-.111112222233333444440000000000D+00, &
      .264665241200000000000000000000D+01, -.378671043061087976727207184637D+00, &
      .104516378011749278484458888919D+01,  .559773594776160811746795939295D+00, &
      .100019582406632651901909339800D+00,  .454219904863173579920523812663D+00, &
      .189511781635593675546652093433D+01,  .582240526465012505902656320160D+00, &
      .164493406684822643647241516665D+01,  .886226925452758013649083741687D+00, &
      -.314159265358979323846264338328D+01, .318309886183790671537767526733D+00, &
      .882395720020380090550940262394D-06, -.282094791773878143474039725759D+00, &
      .187500000000000000000000000000D+01,  .513516668382050295584635612122D-01, &
      .598750000000000000000000000000D+02,  .157079632679489661923132169164D+01, &
      .755006169037464042751871235437D-03,  .422784335098467139393487909918D+00, &
      .230300103429768637527259355045D+01,  .999856618263723706885830759463D+00, &
      .888290707183956735878281870759D+00,  .135335283236612691893999494971D+00, &
      .346930306295801456170933128256D-03,  .786938680574733152792400930048D+00, &
      .631673391775258123291222663623D-01,  .381281566461770916149261183171D+00, &
      .265625000000000000000000000000D+00,  .520499877813046537682746653770D+00, &
      .888388231701707764069578446749D+00,  .424436383502022295934042352455D+00, &
      .337000659742093423383019719632D+00, -.177596771314338304347397013056D+00, &
      .223890779141235668051827454628D+00, -.327579137591465222037734321812D+00, &
      .576724807756873387202448242187D+00,  .510375672649745119596606592612D+00, &
      -.308517625249033780073648984210D+00, .147863143391226844801050675510D+00, &
      -.107032431540937546888370772230D+00, .227958530233606726743720444020D+01, &
      .272398718236044468945442320700D+02,  .159063685463732906338225442450D+01, &
      .243356421424505271991430504400D+02,  .113893872749533435652719574910D+00, &
      .369109833404259427473526100740D-02,  .139865881816522427284598806997D+00, &
      .404461344545216420836502183700D-02,  .308508322553671039533384319255D+00, &
      .183540812609328353073650751820D+00,  .163972266944542356926122903850D+00, &
      .215269289248937659158505143243D+00,  .841568215070771417919124867127D+00, &
      .547807564313518986868201568700D+00,  .600273858788312582936045656600D+00, &
      .103347684706868857317535710603D+01,  .886226925452758013649083741000D+00, &
      .132934038817913702047362561200D+01,  .288023750772146354435952215970D+01, &
      .560499121639792869931128243359D+00,  .672598945967751443917353892000D+00, &
      .964058489220443736281540578570D+00,  .461068504447894558439575873876D+00, &
      .922137008895789116879151747751D+00,  .231693606480833489769125254500D+00, &
      .157259233804704899952660465400D-01,  .293277159129947362450897433147D+00, &
      .219322205128712060862850888400D+00,  .854277043103155493300048798776D+00, &
      .187894150374789500090933504950D+01,  .674892411115630212865414309867D+00, &
      .464750480196092515019775411670D+00,  .249999999999999999999999999880D+00, &
      .735008609300377745369706799000D+00,  .406961787650672979742685260000D+00, &
      .448256669291582953916931735480D+00,  .596347362323194074341078499290D+00, &
      .757342086122175953454414369190D+00,  .757872156141312106043351240000D+00 ]
    !* FIRST EXECUTABLE STATEMENT  DFNCK
    !
    !     Compute functional values
    !
    !     Exercise routines in Category C1.
    !
    y(1) = DBINOM(35,12)
    y(2) = DBINOM(50,15)
    y(3) = DFAC(15)
    y(4) = DFAC(31)
    !
    !     Exercise routines in Category C2
    !
    y(5) = DCBRT(-8.D0)
    y(6) = DCBRT(.995030624365703964475039000000D6)
    !
    !     Exercise routines in Category C4A.
    !
    y(7) = DCOSDG(30.D0)
    y(8) = DCOSDG(135.D0)
    y(9) = DCOT(1.D0)
    y(10) = DCOT(-.5D0)
    y(11) = D9ATN1(.5D0)
    y(12) = D9ATN1(2.D0)
    y(13) = DSINDG(30.D0)
    y(14) = DSINDG(135.D0)
    !
    !     Exercise routines in Category C4B.
    !
    y(15) = DLNREL(.147D0)
    y(16) = DLNREL(-.1D-5)
    y(17) = DEXPRL(.25D-2)
    y(18) = DEXPRL(-.1D-1)
    y(19) = D9LN2R(.5D0)
    y(20) = D9LN2R(1.D0)
    !
    !     Exercise routines in Category C4C.
    !
    y(21) = DACOSH(.100617931649094823747218929626D1)
    y(22) = DACOSH(.115919532755239084628557897777D2)
    y(23) = DASINH(.100000000101295145211538706587D0)
    y(24) = DASINH(-.149999999948240634124264852207D1)
    y(25) = DATANH(-.110657208041383998066515207788D0)
    y(26) = DATANH(.989999999992791300663084082410D0)
    !
    !     Exercise routines in Category C5.
    !
    y(27) = DLI(.5D0)
    y(28) = DLI(2.D0)
    y(29) = DE1(.5D0)
    y(30) = DE1(1.5D0)
    y(31) = DEI(.5D0)
    y(32) = DEI(1.D0)
    y(33) = DSPENC(.5D0)
    y(34) = DSPENC(1.D0)
    y(35) = DGAMMA(1.5D0)
    y(36) = DGAMMA(-.5D0)*DGAMMA(1.5D0)
    y(37) = DGAMR(-1.5D0)*DGAMR(2.5D0)
    y(38) = DGAMR(10.5D0)
    !
    !     Exercise routines in Category C7A.
    !
    y(39) = DPOCH(-.5D0,1.5D0)
    y(40) = DPOCH(.5D0,3.D0)
    y(41) = DPOCH1(.5D0,2.5D0)
    y(42) = DPOCH1(10.5D0,2.D0)
    !
    !     Exercise routines in Category C7B.
    !
    y(43) = DBETA(.5D0,1.5D0)
    y(44) = DBETA(5.5D0,5.5D0)
    !
    !     Exercise routines in Category C7C.
    !
    y(45) = DPSI(2.D0)
    y(46) = DPSI(10.5D0)
    !
    !     Exercise routines in Category C7E.
    !
    y(47) = DGAMI(1.D0,8.85D0)
    y(48) = DGAMI(2.D0,3.75D0)
    y(49) = DGAMIC(1.D0,2.D0)
    y(50) = DGAMIC(2.D0,10.4D0)
    y(51) = DGAMIT(1.D0,.5D0)
    y(52) = DGAMIT(2.D0,3.75D0)
    !
    !     Exercise routines in Category C7F.
    !
    y(53) = DBETAI(.5D0,2.D0,1.5D0)
    y(54) = DBETAI(.25D0,1.5D0,2.D0)
    !
    !     Exercise routines in Category C8A.
    !
    y(55) = DERF(.5D0)
    y(56) = DERF(1.125D0)
    !
    !     Exercise routines in Category C8C.
    !
    y(57) = DDAWS(.5D0)
    y(58) = DDAWS(1.84D0)
    !
    !     Exercise routines in Category C10A1.
    !
    y(59) = DBESJ0(5.D0)
    y(60) = DBESJ0(2.D0)
    y(61) = DBESJ1(5.D0)
    y(62) = DBESJ1(2.D0)
    y(63) = DBESY0(2.D0)
    y(64) = DBESY0(5.D0)
    y(65) = DBESY1(5.D0)
    y(66) = DBESY1(2.D0)
    !
    !     Exercise routines in Category C10B1.
    !
    y(67) = DBESI0(2.D0)
    y(68) = DBESI0(5.D0)
    y(69) = DBESI1(2.D0)
    y(70) = DBESI1(5.D0)
    y(71) = DBESK0(2.D0)
    y(72) = DBESK0(5.D0)
    y(73) = DBESK1(2.D0)
    y(74) = DBESK1(5.D0)
    y(75) = DBSI0E(2.D0)
    y(76) = DBSI0E(5.D0)
    y(77) = DBSI1E(5.D0)
    y(78) = DBSI1E(2.D0)
    y(79) = DBSK0E(2.D0)
    y(80) = DBSK0E(5.D0)
    y(81) = DBSK1E(5.D0)
    y(82) = DBSK1E(2.D0)
    !
    !     Exercise routines in Category C10B3.
    !
    CALL DBSKES(.5D0,2.D0,3,y(83))
    CALL DBSKES(.5D0,5.D0,3,y(86))
    CALL DBESKS(.5D0,1.D0,2,y(89))
    !
    !     Exercise routines in Category C10D.
    !
    y(91) = DAI(.5D0)
    y(92) = DAI(2.5D0)
    y(93) = DAIE(.5D0)
    y(94) = DAIE(2.5D0)
    y(95) = DBI(.5D0)
    y(96) = DBI(1.5D0)
    y(97) = DBIE(.5D0)
    y(98) = DBIE(2.5D0)
    !
    !     Exercise routines in Category C11.
    !
    y(99) = DCHU(1.D0,2.D0,4.D0)
    y(100) = DCHU(5.D0/6.D0,5.D0/3.D0,4.D0/3.D0)
    y(101) = DCHU(.75D0,.75D0,2.5D0)
    y(102) = DCHU(1.D0,1.D0,1.5D0)
    y(103) = DCHU(1.D0,1.D0,1.D0)
    y(104) = DCHU(1.D0,1.D0,-LOG(.5D0))
    y(105) = DCHU(.5D0,.5D0,1.D0)
    !
    !   Check for possible errors
    !
    errmax = D1MACH(4)
    errtol = SQRT(errmax)
    DO i = 1, 105
      abserr = ABS(v(i)-y(i))
      relerr = abserr/ABS(v(i))
      errmax = MAX(relerr,errmax)
      IF ( relerr>errtol.AND.Kprint>=2 ) WRITE (Lun,99001) i, relerr, abserr
      99001 FORMAT (' For I  = ',I3,'  test fails with RELERR  = ',D38.30,&
        '  and ABSERR  = ',D38.30)
    ENDDO
    Ipass = 0
    IF ( errmax<=errtol ) Ipass = 1
    IF ( Ipass/=0.AND.Kprint>=2 ) WRITE (Lun,99002)
    99002 FORMAT (' Double precision Fullerton special function ',' routines o.k.')
    RETURN
  END SUBROUTINE DFNCK
END MODULE TEST03_MOD
!** TEST03
PROGRAM TEST03
  USE TEST03_MOD
  IMPLICIT NONE
  !>
  !***
  !  Driver for testing SLATEC subprograms
  !***
  ! **Library:**   SLATEC
  !***
  ! **Category:**  C
  !***
  ! **Type:**      DOUBLE PRECISION (TEST02-S, TEST03-D, TEST04-C)
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
  !        double precision Fullerton routines
  !
  !***
  ! **References:**  Kirby W. Fong, Thomas H. Jefferson, Tokihiko Suyehiro
  !                 and Lee Walton, Guide to the SLATEC Common Mathema-
  !                 tical Library, April 10, 1990.
  !***
  ! **Routines called:**  DFNCK, I1MACH, XERMAX, XSETF, XSETUN

  !* REVISION HISTORY  (YYMMDD)
  !   890618  DATE WRITTEN
  !   890618  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900524  Cosmetic changes to code.  (WRB)

  INTEGER I1MACH
  INTEGER ipass, kprint, lin, lun, nfail
  !* FIRST EXECUTABLE STATEMENT  TEST03
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
  !     Test double precision Fullerton routines
  !
  CALL DFNCK(lun,kprint,ipass)
  IF ( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF ( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST03 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST03  *************')
  ENDIF
  STOP
END PROGRAM TEST03
