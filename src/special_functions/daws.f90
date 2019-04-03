!** DAWS
REAL FUNCTION DAWS(X)
  IMPLICIT NONE
  !>
  !***
  !  Compute Dawson's function.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C8C
  !***
  ! **Type:**      SINGLE PRECISION (DAWS-S, DDAWS-D)
  !***
  ! **Keywords:**  DAWSON'S FUNCTION, FNLIB, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
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
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  CSEVL, INITS, R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   780401  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   920618  Removed space from variable names.  (RWC, WRB)

  REAL CSEVL, eps, R1MACH, X, y
  INTEGER INITS
  INTEGER, SAVE :: ntdaw, ntdaw2, ntdawa
  REAL, SAVE :: xsml, xbig, xmax
  REAL, PARAMETER :: dawcs(13) = [ -.006351734375145949E0,-.22940714796773869E0, &
    .022130500939084764E0, -.001549265453892985E0, .000084973277156849E0, &
    -.000003828266270972E0, .000000146285480625E0,-.000000004851982381E0, &
    .000000000142146357E0, -.000000000003728836E0, .000000000000088549E0, &
    -.000000000000001920E0, .000000000000000038E0 ]
  REAL, PARAMETER :: daw2cs(29) = [ -.056886544105215527E0, -.31811346996168131E0, &
    .20873845413642237E0,  -.12475409913779131E0,  .067869305186676777E0, &
    -.033659144895270940E0, .015260781271987972E0,-.006348370962596214E0, &
    .002432674092074852E0, -.000862195414910650E0, .000283765733363216E0, &
    -.000087057549874170E0, .000024986849985481E0,-.000006731928676416E0, &
    .000001707857878557E0, -.000000409175512264E0, .000000092828292216E0, &
    -.000000019991403610E0, .000000004096349064E0,-.000000000800324095E0, &
    .000000000149385031E0, -.000000000026687999E0, .000000000004571221E0, &
    -.000000000000751873E0, .000000000000118931E0,-.000000000000018116E0, &
    .000000000000002661E0, -.000000000000000377E0, .000000000000000051E0 ]
  REAL, PARAMETER :: dawacs(26) = [ .01690485637765704E0, .00868325227840695E0, &
    .00024248640424177E0,  .00001261182399572E0, .00000106645331463E0, &
    .00000013581597947E0,  .00000002171042356E0, .00000000286701050E0, &
    -.00000000019013363E0,-.00000000030977804E0,-.00000000010294148E0, &
    -.00000000000626035E0, .00000000000856313E0, .00000000000303304E0, &
    -.00000000000025236E0,-.00000000000042106E0,-.00000000000004431E0, &
    .00000000000004911E0,  .00000000000001235E0,-.00000000000000578E0, &
    -.00000000000000228E0, .00000000000000076E0, .00000000000000038E0, &
    -.00000000000000011E0,-.00000000000000006E0, .00000000000000002E0 ]
  LOGICAL :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  DAWS
  IF ( first ) THEN
    eps = R1MACH(3)
    ntdaw = INITS(dawcs,13,0.1*eps)
    ntdaw2 = INITS(daw2cs,29,0.1*eps)
    ntdawa = INITS(dawacs,26,0.1*eps)
    !
    xsml = SQRT(1.5*eps)
    xbig = SQRT(0.5/eps)
    xmax = EXP(MIN(-LOG(2.*R1MACH(1)),LOG(R1MACH(2)))-1.0)
    first = .FALSE.
  ENDIF
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
    RETURN
  ENDIF
  DAWS = 0.5/X
  IF ( y>xbig ) RETURN
  !
  DAWS = (0.5+CSEVL(32.0/y**2-1.0,dawacs,ntdawa))/X
  RETURN
END FUNCTION DAWS
