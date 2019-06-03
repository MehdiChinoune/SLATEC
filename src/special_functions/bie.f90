!** BIE
REAL FUNCTION BIE(X)
  !>
  !  Calculate the Bairy function for a negative argument and an
  !            exponentially scaled Bairy function for a non-negative
  !            argument.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C10D
  !***
  ! **Type:**      SINGLE PRECISION (BIE-S, DBIE-D)
  !***
  ! **Keywords:**  BAIRY FUNCTION, EXPONENTIALLY SCALED, FNLIB,
  !             SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! Evaluate BI(X) for X .LE. 0  and  BI(X)*EXP(ZETA)  where
  ! ZETA = 2/3 * X**(3/2)  for X .GE. 0.0
  !
  ! Series for BIF        on the interval -1.00000D+00 to  1.00000D+00
  !                                        with weighted error   1.88E-19
  !                                         log weighted error  18.72
  !                               significant figures required  17.74
  !                                    decimal places required  19.20
  !
  ! Series for BIG        on the interval -1.00000D+00 to  1.00000D+00
  !                                        with weighted error   2.61E-17
  !                                         log weighted error  16.58
  !                               significant figures required  15.17
  !                                    decimal places required  17.03
  !
  ! Series for BIF2       on the interval  1.00000D+00 to  8.00000D+00
  !                                        with weighted error   1.11E-17
  !                                         log weighted error  16.95
  !                        approx significant figures required  16.5
  !                                    decimal places required  17.45
  !
  ! Series for BIG2       on the interval  1.00000D+00 to  8.00000D+00
  !                                        with weighted error   1.19E-18
  !                                         log weighted error  17.92
  !                        approx significant figures required  17.2
  !                                    decimal places required  18.42
  !
  ! Series for BIP        on the interval  1.25000D-01 to  3.53553D-01
  !                                        with weighted error   1.91E-17
  !                                         log weighted error  16.72
  !                               significant figures required  15.35
  !                                    decimal places required  17.41
  !
  ! Series for BIP2       on the interval  0.          to  1.25000D-01
  !                                        with weighted error   1.05E-18
  !                                         log weighted error  17.98
  !                               significant figures required  16.74
  !                                    decimal places required  18.71
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  CSEVL, INITS, R1MACH, R9AIMP

  !* REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890206  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  USE service, ONLY : R1MACH
  REAL :: X
  REAL :: sqrtx, theta, xm, z
  INTEGER, SAVE :: nbif, nbig, nbif2, nbig2, nbip, nbip2
  REAL, PARAMETER :: eta = 0.1*R1MACH(3), x3sml = eta**0.3333, &
    x32sml = 1.3104*x3sml**2, xbig = R1MACH(2)**0.6666
  REAL, PARAMETER :: bifcs(9) = [ -.01673021647198664948E0, .1025233583424944561E0, &
    .00170830925073815165E0, .00001186254546774468E0, .00000004493290701779E0, &
    .00000000010698207143E0, .00000000000017480643E0, .00000000000000020810E0, &
    .00000000000000000018E0 ]
  REAL, PARAMETER :: bigcs(8) = [ .02246622324857452E0, .03736477545301955E0, &
    .00044476218957212E0, .00000247080756363E0, .00000000791913533E0, &
    .00000000001649807E0, .00000000000002411E0, .00000000000000002E0 ]
  REAL, PARAMETER :: bif2cs(10) = [ 0.09984572693816041E0, .478624977863005538E0, &
    .0251552119604330118E0, .0005820693885232645E0, .0000074997659644377E0, &
    .0000000613460287034E0, .0000000003462753885E0, .0000000000014288910E0, &
    .0000000000000044962E0, .0000000000000000111E0 ]
  REAL, PARAMETER :: big2cs(10) = [ .033305662145514340E0, .161309215123197068E0, &
    .0063190073096134286E0, .0001187904568162517E0, .0000013045345886200E0, &
    .0000000093741259955E0, .0000000000474580188E0, .0000000000001783107E0, &
    .0000000000000005167E0, .0000000000000000011E0 ]
  REAL, PARAMETER :: bipcs(24) = [ -.08322047477943447E0, .01146118927371174E0, &
    .00042896440718911E0, -.00014906639379950E0,-.00001307659726787E0, &
    .00000632759839610E0, -.00000042226696982E0,-.00000019147186298E0, &
    .00000006453106284E0, -.00000000784485467E0,-.00000000096077216E0, &
    .00000000070004713E0, -.00000000017731789E0, .00000000002272089E0, &
    .00000000000165404E0, -.00000000000185171E0, .00000000000059576E0, &
    -.00000000000012194E0, .00000000000001334E0, .00000000000000172E0, &
    -.00000000000000145E0, .00000000000000049E0,-.00000000000000011E0, &
    .00000000000000001E0 ]
  REAL, PARAMETER :: bip2cs(29) = [ -.113596737585988679E0, .0041381473947881595E0, &
    .0001353470622119332E0,  .0000104273166530153E0, .0000013474954767849E0, &
    .0000001696537405438E0, -.0000000100965008656E0,-.0000000167291194937E0, &
    -.0000000045815364485E0, .0000000003736681366E0, .0000000005766930320E0, &
    .0000000000621812650E0, -.0000000000632941202E0,-.0000000000149150479E0, &
    .0000000000078896213E0,  .0000000000024960513E0,-.0000000000012130075E0, &
    -.0000000000003740493E0, .0000000000002237727E0, .0000000000000474902E0, &
    -.0000000000000452616E0,-.0000000000000030172E0, .0000000000000091058E0, &
    -.0000000000000009814E0,-.0000000000000016429E0, .0000000000000005533E0, &
    .0000000000000002175E0, -.0000000000000001737E0,-.0000000000000000010E0 ]
  REAL, PARAMETER :: atr = 8.7506905708484345E0
  REAL, PARAMETER :: btr = -2.093836321356054E0
  LOGICAL, SAVE :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  BIE
  IF ( first ) THEN
    nbif = INITS(bifcs,9,eta)
    nbig = INITS(bigcs,8,eta)
    nbif2 = INITS(bif2cs,10,eta)
    nbig2 = INITS(big2cs,10,eta)
    nbip = INITS(bipcs,24,eta)
    nbip2 = INITS(bip2cs,29,eta)
    first = .FALSE.
  END IF
  !
  IF ( X<(-1.0) ) THEN
    CALL R9AIMP(X,xm,theta)
    BIE = xm*SIN(theta)
    RETURN
    !
  ELSEIF ( X<=1.0 ) THEN
    z = 0.0
    IF ( ABS(X)>x3sml ) z = X**3
    BIE = 0.625 + CSEVL(z,bifcs,nbif) + X*(0.4375+CSEVL(z,bigcs,nbig))
    IF ( X>x32sml ) BIE = BIE*EXP(-2.0*X*SQRT(X)/3.0)
    RETURN
    !
  ELSEIF ( X<=2.0 ) THEN
    z = (2.0*X**3-9.0)/7.0
    BIE = EXP(-2.0*X*SQRT(X)/3.0)&
      *(1.125+CSEVL(z,bif2cs,nbif2)+X*(0.625+CSEVL(z,big2cs,nbig2)))
    RETURN
    !
  ELSEIF ( X>4.0 ) THEN
    !
    sqrtx = SQRT(X)
    z = -1.0
    IF ( X<xbig ) z = 16.0/(X*sqrtx) - 1.0
    BIE = (0.625+CSEVL(z,bip2cs,nbip2))/SQRT(sqrtx)
    RETURN
  END IF
  sqrtx = SQRT(X)
  z = atr/(X*sqrtx) + btr
  BIE = (0.625+CSEVL(z,bipcs,nbip))/SQRT(sqrtx)
  RETURN
END FUNCTION BIE
