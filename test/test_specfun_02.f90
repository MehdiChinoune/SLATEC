MODULE TEST03_MOD
  USE service, ONLY : SP, DP
  IMPLICIT NONE

CONTAINS
  !** DFNCK
  SUBROUTINE DFNCK(Lun,Kprint,Ipass)
    !> Quick check for the double precision Fullerton
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
    !                    DERF, DEXPRL, DFAC, DGAMI, DGAMIC, DGAMIT,
    !                    DGAMR, DLI, DLNREL, DPOCH, DPOCH1, DPSI, DSINDG,
    !                    DSPENC

    !* REVISION HISTORY  (YYMMDD)
    !   800801  DATE WRITTEN
    !   891115  REVISION DATE from Version 3.2
    !   891120  Checks of remainder of FNLIB routines added and code
    !           reorganized.  (WRB)
    !   900330  Prologue converted to Version 4.0 format.  (BAB)
    !   900727  Added EXTERNAL statement.  (WRB)
    USE slatec, ONLY : D1MACH, D9ATN1, D9LN2R, DAI, DAIE, DBESI0, DBESI1, DBESK0, &
      DBESK1, DBESKS, DBETA, DBETAI, DBI, DBIE, DBINOM, DBSI0E, DBSI1E, DBSK0E, &
      DBSK1E, DBSKES, DCBRT, DCHU, DCOSDG, DCOT, DDAWS, DE1, DEI, DEXPRL, DFAC, &
      DGAMI, DGAMIC, DGAMIT, DGAMR, DLI, DLNREL, DPOCH, DPOCH1, DPSI, DSINDG, DSPENC
    INTEGER :: i, Lun, Kprint, Ipass
    REAL(DP) :: y(105), errmax, errtol, abserr, relerr
    !
    !     Correct values through different calculations are stored in V(*)
    !
    REAL(DP), PARAMETER :: v(87) = [ .834451800000000000000000000000E+09_DP, &
      .225082957512000000000000000000E+13_DP,  .130767436800000000000000000000E+13_DP, &
      .822283865417792281772556288000E+34_DP, -.200000000000000000000000000000E+01_DP, &
      .998340790000000000000000000000E+02_DP,  .866025403784438646763723170753E+00_DP, &
      -.707106781186547524400844362105E+00_DP, .642092615934330703006419986594E+00_DP, &
      -.183048772171245191926801943897E+01_DP,-.290819127993551070285950148310E+00_DP, &
      -.111606410275738687122866817478E+00_DP, .500000000000000000000000000000E+00_DP, &
      .707106781186547524400844362105E+00_DP,  .137149838147233638243285631505E+00_DP, &
      -.100000050000033333358333416027E-05_DP, .100125104231803398984880296644E+01_DP, &
      .995016625083194642609402280122E+00_DP,  .243720864865315055824104923715E+00_DP, &
      .193147180559945309417232121458E+00_DP, -.378671043061087976727207184637E+00_DP, &
      .104516378011749278484458888919E+01_DP,  .559773594776160811746795939295E+00_DP, &
      .100019582406632651901909339800E+00_DP,  .454219904863173579920523812663E+00_DP, &
      .189511781635593675546652093433E+01_DP,  .582240526465012505902656320160E+00_DP, &
      .164493406684822643647241516665E+01_DP, .318309886183790671537767526733E+00_DP, &
      .882395720020380090550940262394E-06_DP, -.282094791773878143474039725759E+00_DP, &
      .187500000000000000000000000000E+01_DP,  .513516668382050295584635612122E-01_DP, &
      .598750000000000000000000000000E+02_DP,  .157079632679489661923132169164E+01_DP, &
      .755006169037464042751871235437E-03_DP,  .422784335098467139393487909918E+00_DP, &
      .230300103429768637527259355045E+01_DP,  .999856618263723706885830759463E+00_DP, &
      .888290707183956735878281870759E+00_DP,  .135335283236612691893999494971E+00_DP, &
      .346930306295801456170933128256E-03_DP,  .786938680574733152792400930048E+00_DP, &
      .631673391775258123291222663623E-01_DP,  .381281566461770916149261183171E+00_DP, &
      .265625000000000000000000000000E+00_DP,  .424436383502022295934042352455E+00_DP, &
      .337000659742093423383019719632E+00_DP, .227958530233606726743720444020E+01_DP, &
      .272398718236044468945442320700E+02_DP,  .159063685463732906338225442450E+01_DP, &
      .243356421424505271991430504400E+02_DP,  .113893872749533435652719574910E+00_DP, &
      .369109833404259427473526100740E-02_DP,  .139865881816522427284598806997E+00_DP, &
      .404461344545216420836502183700E-02_DP,  .308508322553671039533384319255E+00_DP, &
      .183540812609328353073650751820E+00_DP,  .163972266944542356926122903850E+00_DP, &
      .215269289248937659158505143243E+00_DP,  .841568215070771417919124867127E+00_DP, &
      .547807564313518986868201568700E+00_DP,  .600273858788312582936045656600E+00_DP, &
      .103347684706868857317535710603E+01_DP,  .886226925452758013649083741000E+00_DP, &
      .132934038817913702047362561200E+01_DP,  .288023750772146354435952215970E+01_DP, &
      .560499121639792869931128243359E+00_DP,  .672598945967751443917353892000E+00_DP, &
      .964058489220443736281540578570E+00_DP,  .461068504447894558439575873876E+00_DP, &
      .922137008895789116879151747751E+00_DP,  .231693606480833489769125254500E+00_DP, &
      .157259233804704899952660465400E-01_DP,  .293277159129947362450897433147E+00_DP, &
      .219322205128712060862850888400E+00_DP,  .854277043103155493300048798776E+00_DP, &
      .187894150374789500090933504950E+01_DP,  .674892411115630212865414309867E+00_DP, &
      .464750480196092515019775411670E+00_DP,  .249999999999999999999999999880E+00_DP, &
      .735008609300377745369706799000E+00_DP,  .406961787650672979742685260000E+00_DP, &
      .448256669291582953916931735480E+00_DP,  .596347362323194074341078499290E+00_DP, &
      .757342086122175953454414369190E+00_DP,  .757872156141312106043351240000E+00_DP ]
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
    y(5) = DCBRT(-8._DP)
    y(6) = DCBRT(.995030624365703964475039000000E6_DP)
    !
    !     Exercise routines in Category C4A.
    !
    y(7) = DCOSDG(30._DP)
    y(8) = DCOSDG(135._DP)
    y(9) = DCOT(1._DP)
    y(10) = DCOT(-.5_DP)
    y(11) = D9ATN1(.5_DP)
    y(12) = D9ATN1(2._DP)
    y(13) = DSINDG(30._DP)
    y(14) = DSINDG(135._DP)
    !
    !     Exercise routines in Category C4B.
    !
    y(15) = DLNREL(.147_DP)
    y(16) = DLNREL(-.1E-5_DP)
    y(17) = DEXPRL(.25E-2_DP)
    y(18) = DEXPRL(-.1E-1_DP)
    y(19) = D9LN2R(.5_DP)
    y(20) = D9LN2R(1._DP)
    !
    !     Exercise routines in Category C5.
    !
    y(21) = DLI(.5_DP)
    y(22) = DLI(2._DP)
    y(23) = DE1(.5_DP)
    y(24) = DE1(1.5_DP)
    y(25) = DEI(.5_DP)
    y(26) = DEI(1._DP)
    y(27) = DSPENC(.5_DP)
    y(28) = DSPENC(1._DP)
    y(29) = DGAMR(-1.5_DP)*DGAMR(2.5_DP)
    y(30) = DGAMR(10.5_DP)
    !
    !     Exercise routines in Category C7A.
    !
    y(31) = DPOCH(-.5_DP,1.5_DP)
    y(32) = DPOCH(.5_DP,3._DP)
    y(33) = DPOCH1(.5_DP,2.5_DP)
    y(34) = DPOCH1(10.5_DP,2._DP)
    !
    !     Exercise routines in Category C7B.
    !
    y(35) = DBETA(.5_DP,1.5_DP)
    y(36) = DBETA(5.5_DP,5.5_DP)
    !
    !     Exercise routines in Category C7C.
    !
    y(37) = DPSI(2._DP)
    y(38) = DPSI(10.5_DP)
    !
    !     Exercise routines in Category C7E.
    !
    y(39) = DGAMI(1._DP,8.85_DP)
    y(40) = DGAMI(2._DP,3.75_DP)
    y(41) = DGAMIC(1._DP,2._DP)
    y(42) = DGAMIC(2._DP,10.4_DP)
    y(43) = DGAMIT(1._DP,.5_DP)
    y(44) = DGAMIT(2._DP,3.75_DP)
    !
    !     Exercise routines in Category C7F.
    !
    y(45) = DBETAI(.5_DP,2._DP,1.5_DP)
    y(46) = DBETAI(.25_DP,1.5_DP,2._DP)
    !
    !     Exercise routines in Category C8C.
    !
    y(47) = DDAWS(.5_DP)
    y(48) = DDAWS(1.84_DP)
    !
    !     Exercise routines in Category C10B1.
    !
    y(49) = DBESI0(2._DP)
    y(50) = DBESI0(5._DP)
    y(51) = DBESI1(2._DP)
    y(52) = DBESI1(5._DP)
    y(53) = DBESK0(2._DP)
    y(54) = DBESK0(5._DP)
    y(55) = DBESK1(2._DP)
    y(56) = DBESK1(5._DP)
    y(57) = DBSI0E(2._DP)
    y(58) = DBSI0E(5._DP)
    y(59) = DBSI1E(5._DP)
    y(60) = DBSI1E(2._DP)
    y(61) = DBSK0E(2._DP)
    y(62) = DBSK0E(5._DP)
    y(63) = DBSK1E(5._DP)
    y(64) = DBSK1E(2._DP)
    !
    !     Exercise routines in Category C10B3.
    !
    CALL DBSKES(.5_DP,2._DP,3,y(65:67))
    CALL DBSKES(.5_DP,5._DP,3,y(68:70))
    CALL DBESKS(.5_DP,1._DP,2,y(71:72))
    !
    !     Exercise routines in Category C10D.
    !
    y(73) = DAI(.5_DP)
    y(74) = DAI(2.5_DP)
    y(75) = DAIE(.5_DP)
    y(76) = DAIE(2.5_DP)
    y(77) = DBI(.5_DP)
    y(78) = DBI(1.5_DP)
    y(79) = DBIE(.5_DP)
    y(80) = DBIE(2.5_DP)
    !
    !     Exercise routines in Category C11.
    !
    y(81) = DCHU(1._DP,2._DP,4._DP)
    y(82) = DCHU(5._DP/6._DP,5._DP/3._DP,4._DP/3._DP)
    y(83) = DCHU(.75_DP,.75_DP,2.5_DP)
    y(84) = DCHU(1._DP,1._DP,1.5_DP)
    y(85) = DCHU(1._DP,1._DP,1._DP)
    y(86) = DCHU(1._DP,1._DP,-LOG(.5_DP))
    y(87) = DCHU(.5_DP,.5_DP,1._DP)
    !
    !   Check for possible errors
    !
    errmax = D1MACH(4)
    errtol = SQRT(errmax)
    DO i = 1, 87
      abserr = ABS(v(i)-y(i))
      relerr = abserr/ABS(v(i))
      errmax = MAX(relerr,errmax)
      IF( relerr>errtol .AND. Kprint>=2 ) WRITE (Lun,99001) i, relerr, abserr
      99001 FORMAT (' For I  = ',I3,'  test fails with RELERR  = ',D38.30,&
        '  and ABSERR  = ',D38.30)
    END DO
    Ipass = 0
    IF( errmax<=errtol ) Ipass = 1
    IF( Ipass/=0 .AND. Kprint>=2 ) WRITE (Lun,99002)
    99002 FORMAT (' Double precision Fullerton special function ',' routines o.k.')
    RETURN
  END SUBROUTINE DFNCK
END MODULE TEST03_MOD
!** TEST03
PROGRAM TEST03
  USE TEST03_MOD, ONLY : DFNCK
  USE slatec, ONLY : I1MACH, XSETF, XSETUN, XERMAX
  USE common_mod, ONLY : GET_ARGUMENT
  IMPLICIT NONE
  !> Driver for testing SLATEC subprograms
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
  INTEGER :: ipass, kprint, lin, lun, nfail
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
  IF( kprint<=1 ) THEN
    CALL XSETF(0)
  ELSE
    CALL XSETF(1)
  END IF
  !
  !     Test double precision Fullerton routines
  !
  CALL DFNCK(lun,kprint,ipass)
  IF( ipass==0 ) nfail = nfail + 1
  !
  !     Write PASS or FAIL message
  !
  IF( nfail==0 ) THEN
    WRITE (lun,99001)
    99001 FORMAT (/' --------------TEST03 PASSED ALL TESTS----------------')
  ELSE
    WRITE (lun,99002) nfail
    99002 FORMAT (/' ************* WARNING -- ',I5,&
      ' TEST(S) FAILED IN PROGRAM TEST03  *************')
  END IF
  STOP
END PROGRAM TEST03
