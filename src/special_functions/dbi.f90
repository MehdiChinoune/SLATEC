!** DBI
REAL(DP) FUNCTION DBI(X)
  !> Evaluate the Bairy function (the Airy function of the
  !            second kind).
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C10D
  !***
  ! **Type:**      DOUBLE PRECISION (BI-S, DBI-D)
  !***
  ! **Keywords:**  BAIRY FUNCTION, FNLIB, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! DBI(X) calculates the double precision Airy function of the
  ! second kind for double precision argument X.
  !
  ! Series for BIF        on the interval -1.00000E+00 to  1.00000E+00
  !                                        with weighted error   1.45E-32
  !                                         log weighted error  31.84
  !                               significant figures required  30.85
  !                                    decimal places required  32.40
  !
  ! Series for BIG        on the interval -1.00000E+00 to  1.00000E+00
  !                                        with weighted error   1.29E-33
  !                                         log weighted error  32.89
  !                               significant figures required  31.48
  !                                    decimal places required  33.45
  !
  ! Series for BIF2       on the interval  1.00000E+00 to  8.00000E+00
  !                                        with weighted error   6.08E-32
  !                                         log weighted error  31.22
  !                        approx significant figures required  30.8
  !                                    decimal places required  31.80
  !
  ! Series for BIG2       on the interval  1.00000E+00 to  8.00000E+00
  !                                        with weighted error   4.91E-33
  !                                         log weighted error  32.31
  !                        approx significant figures required  31.6
  !                                    decimal places required  32.90
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  D1MACH, D9AIMP, DBIE, DCSEVL, INITDS, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  USE service, ONLY : XERMSG, D1MACH
  REAL(DP) :: X
  REAL(DP) :: theta, xm, z
  INTEGER, SAVE :: nbif, nbig, nbif2, nbig2
  REAL(DP), PARAMETER :: eta = 0.1_DP*D1MACH(3), x3sml = eta**0.3333_DP, &
    xmax = (1.5_SP*LOG(D1MACH(2)))**0.6666_DP
  REAL(DP), PARAMETER :: bifcs(13) = [ -.16730216471986649483537423928176E-1_DP, &
    +.10252335834249445611426362777757E+0_DP, +.17083092507381516539429650242013E-2_DP, &
    +.11862545467744681179216459210040E-4_DP, +.44932907017792133694531887927242E-7_DP, &
    +.10698207143387889067567767663628E-9_DP, +.17480643399771824706010517628573E-12_DP, &
    +.20810231071761711025881891834399E-15_DP, +.18849814695665416509927971733333E-18_DP, &
    +.13425779173097804625882666666666E-21_DP, +.77159593429658887893333333333333E-25_DP, &
    +.36533879617478566399999999999999E-28_DP, +.14497565927953066666666666666666E-31_DP ]
  REAL(DP), PARAMETER :: bigcs(13) = [ +.22466223248574522283468220139024E-1_DP, &
    +.37364775453019545441727561666752E-1_DP, +.44476218957212285696215294326639E-3_DP, &
    +.24708075636329384245494591948882E-5_DP, +.79191353395149635134862426285596E-8_DP, &
    +.16498079851827779880887872402706E-10_DP, +.24119906664835455909247501122841E-13_DP, &
    +.26103736236091436985184781269333E-16_DP, +.21753082977160323853123792000000E-19_DP, &
    +.14386946400390433219483733333333E-22_DP, +.77349125612083468629333333333333E-26_DP, &
    +.34469292033849002666666666666666E-29_DP, +.12938919273216000000000000000000E-32_DP ]
  REAL(DP), PARAMETER :: bif2cs(15) = [ +.0998457269381604104468284257993E+0_DP, &
    +.47862497786300553772211467318231E+0_DP, +.25155211960433011771324415436675E-1_DP, &
    +.58206938852326456396515697872216E-3_DP, +.74997659644377865943861457378217E-5_DP, &
    +.61346028703493836681403010356474E-7_DP, +.34627538851480632900434268733359E-9_DP, &
    +.14288910080270254287770846748931E-11_DP, +.44962704298334641895056472179200E-14_DP, &
    +.11142323065833011708428300106666E-16_DP, +.22304791066175002081517866666666E-19_DP, &
    +.36815778736393142842922666666666E-22_DP, +.50960868449338261333333333333333E-25_DP, &
    +.60003386926288554666666666666666E-28_DP, +.60827497446570666666666666666666E-31_DP ]
  REAL(DP), PARAMETER :: big2cs(15) = [ +.033305662145514340465176188111647E+0_DP, &
    +.161309215123197067613287532084943E+0_DP, +.631900730961342869121615634921173E-2_DP, &
    +.118790456816251736389780192304567E-3_DP, +.130453458862002656147116485012843E-5_DP, &
    +.937412599553521729546809615508936E-8_DP, +.474580188674725153788510169834595E-10_DP, &
    +.178310726509481399800065667560946E-12_DP, +.516759192784958180374276356640000E-15_DP, &
    +.119004508386827125129496251733333E-17_DP, +.222982880666403517277063466666666E-20_DP, &
    +.346551923027689419722666666666666E-23_DP, +.453926336320504514133333333333333E-26_DP, &
    +.507884996513522346666666666666666E-29_DP, +.491020674696533333333333333333333E-32_DP ]
  LOGICAL, SAVE :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  DBI
  IF( first ) THEN
    nbif = INITDS(bifcs,13,eta)
    nbig = INITDS(bigcs,13,eta)
    nbif2 = INITDS(bif2cs,15,eta)
    nbig2 = INITDS(big2cs,15,eta)
    first = .FALSE.
  END IF
  !
  IF( X<(-1._DP) ) THEN
    CALL D9AIMP(X,xm,theta)
    DBI = xm*SIN(theta)
    RETURN
    !
  ELSEIF( X<=1._DP ) THEN
    z = 0._DP
    IF( ABS(X)>x3sml ) z = X**3
    DBI = 0.625_DP + DCSEVL(z,bifcs,nbif) + X*(0.4375_DP+DCSEVL(z,bigcs,nbig))
    RETURN
    !
  ELSEIF( X>2._DP ) THEN
    !
    IF( X>xmax ) CALL XERMSG('DBI','X SO BIG THAT BI OVERFLOWS',1,2)
    !
    DBI = DBIE(X)*EXP(2._DP*X*SQRT(X)/3._DP)
    RETURN
  END IF
  z = (2._DP*X**3-9._DP)/7._DP
  DBI = 1.125_DP + DCSEVL(z,bif2cs,nbif2)&
    + X*(0.625_DP+DCSEVL(z,big2cs,nbig2))
  RETURN
END FUNCTION DBI
