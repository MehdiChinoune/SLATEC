!*==DAWS.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK DAWS
FUNCTION DAWS(X)
  IMPLICIT NONE
  !*--DAWS5
  !*** Start of declarations inserted by SPAG
  REAL CSEVL , daw2cs , dawacs , dawcs , DAWS , eps , R1MACH , X , xbig , &
    xmax , xsml , y
  INTEGER INITS , ntdaw , ntdaw2 , ntdawa
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  DAWS
  !***PURPOSE  Compute Dawson's function.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C8C
  !***TYPE      SINGLE PRECISION (DAWS-S, DDAWS-D)
  !***KEYWORDS  DAWSON'S FUNCTION, FNLIB, SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! DAWS(X) calculates Dawson's integral for real argument X.
  !
  ! Series for DAW        on the interval  0.          to  1.00000D+00
  !                                        with weighted error   3.83E-17
  !                                         log weighted error  16.42
  !                               significant figures required  15.78
  !                                    decimal places required  16.97
  !
  ! Series for DAW2       on the interval  0.          to  1.60000D+01
  !                                        with weighted error   5.17E-17
  !                                         log weighted error  16.29
  !                               significant figures required  15.90
  !                                    decimal places required  17.02
  !
  ! Series for DAWA       on the interval  0.          to  6.25000D-02
  !                                        with weighted error   2.24E-17
  !                                         log weighted error  16.65
  !                               significant figures required  14.73
  !                                    decimal places required  17.36
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  CSEVL, INITS, R1MACH, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   780401  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   920618  Removed space from variable names.  (RWC, WRB)
  !***END PROLOGUE  DAWS
  DIMENSION dawcs(13) , daw2cs(29) , dawacs(26)
  LOGICAL first
  SAVE dawcs , daw2cs , dawacs , ntdaw , ntdaw2 , ntdawa , xsml , xbig , &
    xmax , first
  DATA dawcs(1)/ - .006351734375145949E0/
  DATA dawcs(2)/ - .22940714796773869E0/
  DATA dawcs(3)/.022130500939084764E0/
  DATA dawcs(4)/ - .001549265453892985E0/
  DATA dawcs(5)/.000084973277156849E0/
  DATA dawcs(6)/ - .000003828266270972E0/
  DATA dawcs(7)/.000000146285480625E0/
  DATA dawcs(8)/ - .000000004851982381E0/
  DATA dawcs(9)/.000000000142146357E0/
  DATA dawcs(10)/ - .000000000003728836E0/
  DATA dawcs(11)/.000000000000088549E0/
  DATA dawcs(12)/ - .000000000000001920E0/
  DATA dawcs(13)/.000000000000000038E0/
  DATA daw2cs(1)/ - .056886544105215527E0/
  DATA daw2cs(2)/ - .31811346996168131E0/
  DATA daw2cs(3)/.20873845413642237E0/
  DATA daw2cs(4)/ - .12475409913779131E0/
  DATA daw2cs(5)/.067869305186676777E0/
  DATA daw2cs(6)/ - .033659144895270940E0/
  DATA daw2cs(7)/.015260781271987972E0/
  DATA daw2cs(8)/ - .006348370962596214E0/
  DATA daw2cs(9)/.002432674092074852E0/
  DATA daw2cs(10)/ - .000862195414910650E0/
  DATA daw2cs(11)/.000283765733363216E0/
  DATA daw2cs(12)/ - .000087057549874170E0/
  DATA daw2cs(13)/.000024986849985481E0/
  DATA daw2cs(14)/ - .000006731928676416E0/
  DATA daw2cs(15)/.000001707857878557E0/
  DATA daw2cs(16)/ - .000000409175512264E0/
  DATA daw2cs(17)/.000000092828292216E0/
  DATA daw2cs(18)/ - .000000019991403610E0/
  DATA daw2cs(19)/.000000004096349064E0/
  DATA daw2cs(20)/ - .000000000800324095E0/
  DATA daw2cs(21)/.000000000149385031E0/
  DATA daw2cs(22)/ - .000000000026687999E0/
  DATA daw2cs(23)/.000000000004571221E0/
  DATA daw2cs(24)/ - .000000000000751873E0/
  DATA daw2cs(25)/.000000000000118931E0/
  DATA daw2cs(26)/ - .000000000000018116E0/
  DATA daw2cs(27)/.000000000000002661E0/
  DATA daw2cs(28)/ - .000000000000000377E0/
  DATA daw2cs(29)/.000000000000000051E0/
  DATA dawacs(1)/.01690485637765704E0/
  DATA dawacs(2)/.00868325227840695E0/
  DATA dawacs(3)/.00024248640424177E0/
  DATA dawacs(4)/.00001261182399572E0/
  DATA dawacs(5)/.00000106645331463E0/
  DATA dawacs(6)/.00000013581597947E0/
  DATA dawacs(7)/.00000002171042356E0/
  DATA dawacs(8)/.00000000286701050E0/
  DATA dawacs(9)/ - .00000000019013363E0/
  DATA dawacs(10)/ - .00000000030977804E0/
  DATA dawacs(11)/ - .00000000010294148E0/
  DATA dawacs(12)/ - .00000000000626035E0/
  DATA dawacs(13)/.00000000000856313E0/
  DATA dawacs(14)/.00000000000303304E0/
  DATA dawacs(15)/ - .00000000000025236E0/
  DATA dawacs(16)/ - .00000000000042106E0/
  DATA dawacs(17)/ - .00000000000004431E0/
  DATA dawacs(18)/.00000000000004911E0/
  DATA dawacs(19)/.00000000000001235E0/
  DATA dawacs(20)/ - .00000000000000578E0/
  DATA dawacs(21)/ - .00000000000000228E0/
  DATA dawacs(22)/.00000000000000076E0/
  DATA dawacs(23)/.00000000000000038E0/
  DATA dawacs(24)/ - .00000000000000011E0/
  DATA dawacs(25)/ - .00000000000000006E0/
  DATA dawacs(26)/.00000000000000002E0/
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  DAWS
  IF ( first ) THEN
    eps = R1MACH(3)
    ntdaw = INITS(dawcs,13,0.1*eps)
    ntdaw2 = INITS(daw2cs,29,0.1*eps)
    ntdawa = INITS(dawacs,26,0.1*eps)
    !
    xsml = SQRT(1.5*eps)
    xbig = SQRT(0.5/eps)
    xmax = EXP(MIN(-LOG(2.*R1MACH(1)),LOG(R1MACH(2)))-1.0)
  ENDIF
  first = .FALSE.
  !
  y = ABS(X)
  IF ( y<=1.0 ) THEN
    !
    DAWS = X
    IF ( y<=xsml ) RETURN
    !
    DAWS = X*(0.75+CSEVL(2.0*y*y-1.0,dawcs,ntdaw))
    RETURN
    !
  ELSEIF ( y<=4.0 ) THEN
    DAWS = X*(0.25+CSEVL(0.125*y*y-1.0,daw2cs,ntdaw2))
    RETURN
    !
  ELSEIF ( y>xmax ) THEN
    !
    CALL XERMSG('SLATEC','DAWS','ABS(X) SO LARGE DAWS UNDERFLOWS',1,1)
    DAWS = 0.0
    GOTO 99999
  ENDIF
  DAWS = 0.5/X
  IF ( y>xbig ) RETURN
  !
  DAWS = (0.5+CSEVL(32.0/y**2-1.0,dawacs,ntdawa))/X
  RETURN
  !
  99999 END FUNCTION DAWS
