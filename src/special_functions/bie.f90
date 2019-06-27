!** BIE
REAL(SP) ELEMENTAL FUNCTION BIE(X)
  !> Calculate the Bairy function for a negative argument and an
  !  exponentially scaled Bairy function for a non-negative argument.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C10D
  !***
  ! **Type:**      SINGLE PRECISION (BIE-S, DBIE-D)
  !***
  ! **Keywords:**  BAIRY FUNCTION, EXPONENTIALLY SCALED, FNLIB, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! Evaluate BI(X) for X <= 0  and  BI(X)*EXP(ZETA)  where
  ! ZETA = 2/3 * X**(3/2)  for X >= 0.0
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
  REAL(SP), INTENT(IN) :: X
  REAL(SP) :: sqrtx, theta, xm, z
  INTEGER, PARAMETER :: nbif = 5, nbig = 5, nbif2 = 6, nbig2 = 6, nbip = 10, nbip2 = 8
  REAL(SP), PARAMETER :: eta = 0.1_SP*R1MACH(3), x3sml = eta**0.3333_SP, &
    x32sml = 1.3104_SP*x3sml**2, xbig = R1MACH(2)**0.6666_SP
  REAL(SP), PARAMETER :: bifcs(9) = [ -.01673021647198664948_SP, .1025233583424944561_SP, &
    .00170830925073815165_SP, .00001186254546774468_SP, .00000004493290701779_SP, &
    .00000000010698207143_SP, .00000000000017480643_SP, .00000000000000020810_SP, &
    .00000000000000000018_SP ]
  REAL(SP), PARAMETER :: bigcs(8) = [ .02246622324857452_SP, .03736477545301955_SP, &
    .00044476218957212_SP, .00000247080756363_SP, .00000000791913533_SP, &
    .00000000001649807_SP, .00000000000002411_SP, .00000000000000002_SP ]
  REAL(SP), PARAMETER :: bif2cs(10) = [ 0.09984572693816041_SP, .478624977863005538_SP, &
    .0251552119604330118_SP, .0005820693885232645_SP, .0000074997659644377_SP, &
    .0000000613460287034_SP, .0000000003462753885_SP, .0000000000014288910_SP, &
    .0000000000000044962_SP, .0000000000000000111_SP ]
  REAL(SP), PARAMETER :: big2cs(10) = [ .033305662145514340_SP, .161309215123197068_SP, &
    .0063190073096134286_SP, .0001187904568162517_SP, .0000013045345886200_SP, &
    .0000000093741259955_SP, .0000000000474580188_SP, .0000000000001783107_SP, &
    .0000000000000005167_SP, .0000000000000000011_SP ]
  REAL(SP), PARAMETER :: bipcs(24) = [ -.08322047477943447_SP, .01146118927371174_SP, &
    .00042896440718911_SP, -.00014906639379950_SP,-.00001307659726787_SP, &
    .00000632759839610_SP, -.00000042226696982_SP,-.00000019147186298_SP, &
    .00000006453106284_SP, -.00000000784485467_SP,-.00000000096077216_SP, &
    .00000000070004713_SP, -.00000000017731789_SP, .00000000002272089_SP, &
    .00000000000165404_SP, -.00000000000185171_SP, .00000000000059576_SP, &
    -.00000000000012194_SP, .00000000000001334_SP, .00000000000000172_SP, &
    -.00000000000000145_SP, .00000000000000049_SP,-.00000000000000011_SP, &
    .00000000000000001_SP ]
  REAL(SP), PARAMETER :: bip2cs(29) = [ -.113596737585988679_SP, .0041381473947881595_SP, &
    .0001353470622119332_SP,  .0000104273166530153_SP, .0000013474954767849_SP, &
    .0000001696537405438_SP, -.0000000100965008656_SP,-.0000000167291194937_SP, &
    -.0000000045815364485_SP, .0000000003736681366_SP, .0000000005766930320_SP, &
    .0000000000621812650_SP, -.0000000000632941202_SP,-.0000000000149150479_SP, &
    .0000000000078896213_SP,  .0000000000024960513_SP,-.0000000000012130075_SP, &
    -.0000000000003740493_SP, .0000000000002237727_SP, .0000000000000474902_SP, &
    -.0000000000000452616_SP,-.0000000000000030172_SP, .0000000000000091058_SP, &
    -.0000000000000009814_SP,-.0000000000000016429_SP, .0000000000000005533_SP, &
    .0000000000000002175_SP, -.0000000000000001737_SP,-.0000000000000000010_SP ]
  REAL(SP), PARAMETER :: atr = 8.7506905708484345_SP
  REAL(SP), PARAMETER :: btr = -2.093836321356054_SP
  !* FIRST EXECUTABLE STATEMENT  BIE
  ! nbif = INITS(bifcs,eta)
  ! nbig = INITS(bigcs,eta)
  ! nbif2 = INITS(bif2cs,eta)
  ! nbig2 = INITS(big2cs,eta)
  ! nbip = INITS(bipcs,eta)
  ! nbip2 = INITS(bip2cs,eta)
  !
  IF( X<(-1._SP) ) THEN
    CALL R9AIMP(X,xm,theta)
    BIE = xm*SIN(theta)
  ELSEIF( X<=1._SP ) THEN
    z = 0._SP
    IF( ABS(X)>x3sml ) z = X**3
    BIE = 0.625_SP + CSEVL(z,bifcs(1:nbif)) + X*(0.4375_SP+CSEVL(z,bigcs(1:nbig)))
    IF( X>x32sml ) BIE = BIE*EXP(-2._SP*X*SQRT(X)/3._SP)
  ELSEIF( X<=2._SP ) THEN
    z = (2._SP*X**3-9._SP)/7._SP
    BIE = EXP(-2._SP*X*SQRT(X)/3._SP)&
      *(1.125_SP+CSEVL(z,bif2cs(1:nbif2))+X*(0.625_SP+CSEVL(z,big2cs(1:nbig2))))
  ELSEIF( X<=4._SP ) THEN
    sqrtx = SQRT(X)
    z = atr/(X*sqrtx) + btr
    BIE = (0.625_SP+CSEVL(z,bipcs(1:nbip)))/SQRT(sqrtx)
  ELSE
    sqrtx = SQRT(X)
    z = -1._SP
    IF( X<xbig ) z = 16._SP/(X*sqrtx) - 1._SP
    BIE = (0.625_SP+CSEVL(z,bip2cs(1:nbip2)))/SQRT(sqrtx)
  END IF

  RETURN
END FUNCTION BIE