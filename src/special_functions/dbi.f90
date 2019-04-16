!** DBI
REAL(8) FUNCTION DBI(X)
  !>
  !***
  !  Evaluate the Bairy function (the Airy function of the
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

  REAL eta
  REAL(8) :: X, theta, xm, z
  INTEGER, SAVE :: nbif, nbig, nbif2, nbig2
  REAL(8), SAVE :: x3sml, xmax
  REAL(8), PARAMETER :: bifcs(13) = [ -.16730216471986649483537423928176D-1, &
    +.10252335834249445611426362777757D+0, +.17083092507381516539429650242013D-2, &
    +.11862545467744681179216459210040D-4, +.44932907017792133694531887927242D-7, &
    +.10698207143387889067567767663628D-9, +.17480643399771824706010517628573D-12, &
    +.20810231071761711025881891834399D-15, +.18849814695665416509927971733333D-18, &
    +.13425779173097804625882666666666D-21, +.77159593429658887893333333333333D-25, &
    +.36533879617478566399999999999999D-28, +.14497565927953066666666666666666D-31 ]
  REAL(8), PARAMETER :: bigcs(13) = [ +.22466223248574522283468220139024D-1, &
    +.37364775453019545441727561666752D-1, +.44476218957212285696215294326639D-3, &
    +.24708075636329384245494591948882D-5, +.79191353395149635134862426285596D-8, &
    +.16498079851827779880887872402706D-10, +.24119906664835455909247501122841D-13, &
    +.26103736236091436985184781269333D-16, +.21753082977160323853123792000000D-19, &
    +.14386946400390433219483733333333D-22, +.77349125612083468629333333333333D-26, &
    +.34469292033849002666666666666666D-29, +.12938919273216000000000000000000D-32 ]
  REAL(8), PARAMETER :: bif2cs(15) = [ +.0998457269381604104468284257993D+0, &
    +.47862497786300553772211467318231D+0, +.25155211960433011771324415436675D-1, &
    +.58206938852326456396515697872216D-3, +.74997659644377865943861457378217D-5, &
    +.61346028703493836681403010356474D-7, +.34627538851480632900434268733359D-9, &
    +.14288910080270254287770846748931D-11, +.44962704298334641895056472179200D-14, &
    +.11142323065833011708428300106666D-16, +.22304791066175002081517866666666D-19, &
    +.36815778736393142842922666666666D-22, +.50960868449338261333333333333333D-25, &
    +.60003386926288554666666666666666D-28, +.60827497446570666666666666666666D-31 ]
  REAL(8), PARAMETER :: big2cs(15) = [ +.033305662145514340465176188111647D+0, &
    +.161309215123197067613287532084943D+0, +.631900730961342869121615634921173D-2, &
    +.118790456816251736389780192304567D-3, +.130453458862002656147116485012843D-5, &
    +.937412599553521729546809615508936D-8, +.474580188674725153788510169834595D-10, &
    +.178310726509481399800065667560946D-12, +.516759192784958180374276356640000D-15, &
    +.119004508386827125129496251733333D-17, +.222982880666403517277063466666666D-20, &
    +.346551923027689419722666666666666D-23, +.453926336320504514133333333333333D-26, &
    +.507884996513522346666666666666666D-29, +.491020674696533333333333333333333D-32 ]
  LOGICAL :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  DBI
  IF ( first ) THEN
    eta = 0.1*REAL(D1MACH(3))
    nbif = INITDS(bifcs,13,eta)
    nbig = INITDS(bigcs,13,eta)
    nbif2 = INITDS(bif2cs,15,eta)
    nbig2 = INITDS(big2cs,15,eta)
    !
    x3sml = eta**0.3333
    xmax = (1.5*LOG(D1MACH(2)))**0.6666D0
    first = .FALSE.
  END IF
  !
  IF ( X<(-1.0D0) ) THEN
    CALL D9AIMP(X,xm,theta)
    DBI = xm*SIN(theta)
    RETURN
    !
  ELSEIF ( X<=1.0D0 ) THEN
    z = 0.D0
    IF ( ABS(X)>x3sml ) z = X**3
    DBI = 0.625 + DCSEVL(z,bifcs,nbif) + X*(0.4375D0+DCSEVL(z,bigcs,nbig))
    RETURN
    !
  ELSEIF ( X>2.0D0 ) THEN
    !
    IF ( X>xmax ) CALL XERMSG('SLATEC','DBI','X SO BIG THAT BI OVERFLOWS',1,2)
    !
    DBI = DBIE(X)*EXP(2.0D0*X*SQRT(X)/3.0D0)
    RETURN
  END IF
  z = (2.0D0*X**3-9.0D0)/7.D0
  DBI = 1.125D0 + DCSEVL(z,bif2cs,nbif2)&
    + X*(0.625D0+DCSEVL(z,big2cs,nbig2))
  RETURN
END FUNCTION DBI
