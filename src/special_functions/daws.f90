!** DAWS
REAL(SP) ELEMENTAL FUNCTION DAWS(X)
  !> Compute Dawson's function.
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
  USE service, ONLY : eps_2_sp, tiny_sp, huge_sp
  !
  REAL(SP), INTENT(IN) :: X
  !
  REAL(SP) :: y
  INTEGER, PARAMETER :: ntdaw = 7, ntdaw2 = 18, ntdawa = 7
  REAL(SP), PARAMETER :: eps = eps_2_sp, xsml = SQRT(1.5_SP*eps), xbig = SQRT(0.5_SP/eps), &
    xmax = EXP(MIN(-LOG(2._SP*tiny_sp),LOG(huge_sp))-1._SP)
  REAL(SP), PARAMETER :: dawcs(13) = [ -.006351734375145949_SP,-.22940714796773869_SP, &
    .022130500939084764_SP, -.001549265453892985_SP, .000084973277156849_SP, &
    -.000003828266270972_SP, .000000146285480625_SP,-.000000004851982381_SP, &
    .000000000142146357_SP, -.000000000003728836_SP, .000000000000088549_SP, &
    -.000000000000001920_SP, .000000000000000038_SP ]
  REAL(SP), PARAMETER :: daw2cs(29) = [ -.056886544105215527_SP, -.31811346996168131_SP, &
    .20873845413642237_SP,  -.12475409913779131_SP,  .067869305186676777_SP, &
    -.033659144895270940_SP, .015260781271987972_SP,-.006348370962596214_SP, &
    .002432674092074852_SP, -.000862195414910650_SP, .000283765733363216_SP, &
    -.000087057549874170_SP, .000024986849985481_SP,-.000006731928676416_SP, &
    .000001707857878557_SP, -.000000409175512264_SP, .000000092828292216_SP, &
    -.000000019991403610_SP, .000000004096349064_SP,-.000000000800324095_SP, &
    .000000000149385031_SP, -.000000000026687999_SP, .000000000004571221_SP, &
    -.000000000000751873_SP, .000000000000118931_SP,-.000000000000018116_SP, &
    .000000000000002661_SP, -.000000000000000377_SP, .000000000000000051_SP ]
  REAL(SP), PARAMETER :: dawacs(26) = [ .01690485637765704_SP, .00868325227840695_SP, &
    .00024248640424177_SP,  .00001261182399572_SP, .00000106645331463_SP, &
    .00000013581597947_SP,  .00000002171042356_SP, .00000000286701050_SP, &
    -.00000000019013363_SP,-.00000000030977804_SP,-.00000000010294148_SP, &
    -.00000000000626035_SP, .00000000000856313_SP, .00000000000303304_SP, &
    -.00000000000025236_SP,-.00000000000042106_SP,-.00000000000004431_SP, &
    .00000000000004911_SP,  .00000000000001235_SP,-.00000000000000578_SP, &
    -.00000000000000228_SP, .00000000000000076_SP, .00000000000000038_SP, &
    -.00000000000000011_SP,-.00000000000000006_SP, .00000000000000002_SP ]
  !* FIRST EXECUTABLE STATEMENT  DAWS
  ! ntdaw = INITS(dawcs,0.1_SP*eps)
  ! ntdaw2 = INITS(daw2cs,0.1_SP*eps)
  ! ntdawa = INITS(dawacs,0.1_SP*eps)
  !
  y = ABS(X)
  IF( y<=xsml ) THEN
    DAWS = X
  ELSEIF( y<=1._SP ) THEN
    DAWS = X*(0.75_SP+CSEVL(2._SP*y*y-1._SP,dawcs(1:ntdaw)))
  ELSEIF( y<=4._SP ) THEN
    DAWS = X*(0.25_SP+CSEVL(0.125_SP*y*y-1._SP,daw2cs(1:ntdaw2)))
  ELSEIF( y<=xbig ) THEN
    DAWS = (0.5_SP+CSEVL(32._SP/y**2-1._SP,dawacs(1:ntdawa)))/X
  ELSEIF( y<=xmax ) THEN
    DAWS = 0.5_SP/X
  ELSE
    ! CALL XERMSG('DAWS','ABS(X) SO LARGE DAWS UNDERFLOWS',1,1)
    DAWS = 0._SP
  END IF

  RETURN
END FUNCTION DAWS