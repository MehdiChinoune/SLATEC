!*==DBI.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK DBI
REAL(8) FUNCTION DBI(X)
  IMPLICIT NONE
  !*--DBI5
  !*** Start of declarations inserted by SPAG
  REAL eta
  INTEGER INITDS , nbif , nbif2 , nbig , nbig2
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  DBI
  !***PURPOSE  Evaluate the Bairy function (the Airy function of the
  !            second kind).
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C10D
  !***TYPE      DOUBLE PRECISION (BI-S, DBI-D)
  !***KEYWORDS  BAIRY FUNCTION, FNLIB, SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
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
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  D1MACH, D9AIMP, DBIE, DCSEVL, INITDS, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !***END PROLOGUE  DBI
  REAL(8) :: X , bifcs(13) , bigcs(13) , bif2cs(15) , big2cs(15) , &
    theta , xm , xmax , x3sml , z , D1MACH , DCSEVL , DBIE
  LOGICAL first
  SAVE bifcs , bigcs , bif2cs , big2cs , nbif , nbig , nbif2 , nbig2 , &
    x3sml , xmax , first
  DATA bifcs(1)/ - .16730216471986649483537423928176D-1/
  DATA bifcs(2)/ + .10252335834249445611426362777757D+0/
  DATA bifcs(3)/ + .17083092507381516539429650242013D-2/
  DATA bifcs(4)/ + .11862545467744681179216459210040D-4/
  DATA bifcs(5)/ + .44932907017792133694531887927242D-7/
  DATA bifcs(6)/ + .10698207143387889067567767663628D-9/
  DATA bifcs(7)/ + .17480643399771824706010517628573D-12/
  DATA bifcs(8)/ + .20810231071761711025881891834399D-15/
  DATA bifcs(9)/ + .18849814695665416509927971733333D-18/
  DATA bifcs(10)/ + .13425779173097804625882666666666D-21/
  DATA bifcs(11)/ + .77159593429658887893333333333333D-25/
  DATA bifcs(12)/ + .36533879617478566399999999999999D-28/
  DATA bifcs(13)/ + .14497565927953066666666666666666D-31/
  DATA bigcs(1)/ + .22466223248574522283468220139024D-1/
  DATA bigcs(2)/ + .37364775453019545441727561666752D-1/
  DATA bigcs(3)/ + .44476218957212285696215294326639D-3/
  DATA bigcs(4)/ + .24708075636329384245494591948882D-5/
  DATA bigcs(5)/ + .79191353395149635134862426285596D-8/
  DATA bigcs(6)/ + .16498079851827779880887872402706D-10/
  DATA bigcs(7)/ + .24119906664835455909247501122841D-13/
  DATA bigcs(8)/ + .26103736236091436985184781269333D-16/
  DATA bigcs(9)/ + .21753082977160323853123792000000D-19/
  DATA bigcs(10)/ + .14386946400390433219483733333333D-22/
  DATA bigcs(11)/ + .77349125612083468629333333333333D-26/
  DATA bigcs(12)/ + .34469292033849002666666666666666D-29/
  DATA bigcs(13)/ + .12938919273216000000000000000000D-32/
  DATA bif2cs(1)/ + .0998457269381604104468284257993D+0/
  DATA bif2cs(2)/ + .47862497786300553772211467318231D+0/
  DATA bif2cs(3)/ + .25155211960433011771324415436675D-1/
  DATA bif2cs(4)/ + .58206938852326456396515697872216D-3/
  DATA bif2cs(5)/ + .74997659644377865943861457378217D-5/
  DATA bif2cs(6)/ + .61346028703493836681403010356474D-7/
  DATA bif2cs(7)/ + .34627538851480632900434268733359D-9/
  DATA bif2cs(8)/ + .14288910080270254287770846748931D-11/
  DATA bif2cs(9)/ + .44962704298334641895056472179200D-14/
  DATA bif2cs(10)/ + .11142323065833011708428300106666D-16/
  DATA bif2cs(11)/ + .22304791066175002081517866666666D-19/
  DATA bif2cs(12)/ + .36815778736393142842922666666666D-22/
  DATA bif2cs(13)/ + .50960868449338261333333333333333D-25/
  DATA bif2cs(14)/ + .60003386926288554666666666666666D-28/
  DATA bif2cs(15)/ + .60827497446570666666666666666666D-31/
  DATA big2cs(1)/ + .033305662145514340465176188111647D+0/
  DATA big2cs(2)/ + .161309215123197067613287532084943D+0/
  DATA big2cs(3)/ + .631900730961342869121615634921173D-2/
  DATA big2cs(4)/ + .118790456816251736389780192304567D-3/
  DATA big2cs(5)/ + .130453458862002656147116485012843D-5/
  DATA big2cs(6)/ + .937412599553521729546809615508936D-8/
  DATA big2cs(7)/ + .474580188674725153788510169834595D-10/
  DATA big2cs(8)/ + .178310726509481399800065667560946D-12/
  DATA big2cs(9)/ + .516759192784958180374276356640000D-15/
  DATA big2cs(10)/ + .119004508386827125129496251733333D-17/
  DATA big2cs(11)/ + .222982880666403517277063466666666D-20/
  DATA big2cs(12)/ + .346551923027689419722666666666666D-23/
  DATA big2cs(13)/ + .453926336320504514133333333333333D-26/
  DATA big2cs(14)/ + .507884996513522346666666666666666D-29/
  DATA big2cs(15)/ + .491020674696533333333333333333333D-32/
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  DBI
  IF ( first ) THEN
    eta = 0.1*REAL(D1MACH(3))
    nbif = INITDS(bifcs,13,eta)
    nbig = INITDS(bigcs,13,eta)
    nbif2 = INITDS(bif2cs,15,eta)
    nbig2 = INITDS(big2cs,15,eta)
    !
    x3sml = eta**0.3333
    xmax = (1.5*LOG(D1MACH(2)))**0.6666D0
  ENDIF
  first = .FALSE.
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
    IF ( X>xmax ) CALL XERMSG('SLATEC','DBI','X SO BIG THAT BI OVERFLOWS',1,&
      2)
    !
    DBI = DBIE(X)*EXP(2.0D0*X*SQRT(X)/3.0D0)
    GOTO 99999
  ENDIF
  z = (2.0D0*X**3-9.0D0)/7.D0
  DBI = 1.125D0 + DCSEVL(z,bif2cs,nbif2)&
    + X*(0.625D0+DCSEVL(z,big2cs,nbig2))
  RETURN
  !
  99999 END FUNCTION DBI
