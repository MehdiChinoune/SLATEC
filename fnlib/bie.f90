!*==BIE.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK BIE
FUNCTION BIE(X)
  IMPLICIT NONE
  !*--BIE5
  !*** Start of declarations inserted by SPAG
  REAL atr, BIE, bif2cs, bifcs, big2cs, bigcs, bip2cs, bipcs, btr, &
    CSEVL, eta, R1MACH, sqrtx, theta, X, x32sml, x3sml, xbig, &
    xm, z
  INTEGER INITS, nbif, nbif2, nbig, nbig2, nbip, nbip2
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  BIE
  !***PURPOSE  Calculate the Bairy function for a negative argument and an
  !            exponentially scaled Bairy function for a non-negative
  !            argument.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C10D
  !***TYPE      SINGLE PRECISION (BIE-S, DBIE-D)
  !***KEYWORDS  BAIRY FUNCTION, EXPONENTIALLY SCALED, FNLIB,
  !             SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
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
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  CSEVL, INITS, R1MACH, R9AIMP
  !***REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890206  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !***END PROLOGUE  BIE
  LOGICAL first
  DIMENSION bifcs(9), bigcs(8), bif2cs(10), big2cs(10), bipcs(24), &
    bip2cs(29)
  SAVE bifcs, bigcs, bif2cs, big2cs, bipcs, bip2cs, atr, btr, nbif, &
    nbig, nbif2, nbig2, nbip, nbip2, x3sml, x32sml, xbig, first
  DATA bifcs(1)/ - .01673021647198664948E0/
  DATA bifcs(2)/.1025233583424944561E0/
  DATA bifcs(3)/.00170830925073815165E0/
  DATA bifcs(4)/.00001186254546774468E0/
  DATA bifcs(5)/.00000004493290701779E0/
  DATA bifcs(6)/.00000000010698207143E0/
  DATA bifcs(7)/.00000000000017480643E0/
  DATA bifcs(8)/.00000000000000020810E0/
  DATA bifcs(9)/.00000000000000000018E0/
  DATA bigcs(1)/.02246622324857452E0/
  DATA bigcs(2)/.03736477545301955E0/
  DATA bigcs(3)/.00044476218957212E0/
  DATA bigcs(4)/.00000247080756363E0/
  DATA bigcs(5)/.00000000791913533E0/
  DATA bigcs(6)/.00000000001649807E0/
  DATA bigcs(7)/.00000000000002411E0/
  DATA bigcs(8)/.00000000000000002E0/
  DATA bif2cs(1)/0.09984572693816041E0/
  DATA bif2cs(2)/.478624977863005538E0/
  DATA bif2cs(3)/.0251552119604330118E0/
  DATA bif2cs(4)/.0005820693885232645E0/
  DATA bif2cs(5)/.0000074997659644377E0/
  DATA bif2cs(6)/.0000000613460287034E0/
  DATA bif2cs(7)/.0000000003462753885E0/
  DATA bif2cs(8)/.0000000000014288910E0/
  DATA bif2cs(9)/.0000000000000044962E0/
  DATA bif2cs(10)/.0000000000000000111E0/
  DATA big2cs(1)/.033305662145514340E0/
  DATA big2cs(2)/.161309215123197068E0/
  DATA big2cs(3)/.0063190073096134286E0/
  DATA big2cs(4)/.0001187904568162517E0/
  DATA big2cs(5)/.0000013045345886200E0/
  DATA big2cs(6)/.0000000093741259955E0/
  DATA big2cs(7)/.0000000000474580188E0/
  DATA big2cs(8)/.0000000000001783107E0/
  DATA big2cs(9)/.0000000000000005167E0/
  DATA big2cs(10)/.0000000000000000011E0/
  DATA bipcs(1)/ - .08322047477943447E0/
  DATA bipcs(2)/.01146118927371174E0/
  DATA bipcs(3)/.00042896440718911E0/
  DATA bipcs(4)/ - .00014906639379950E0/
  DATA bipcs(5)/ - .00001307659726787E0/
  DATA bipcs(6)/.00000632759839610E0/
  DATA bipcs(7)/ - .00000042226696982E0/
  DATA bipcs(8)/ - .00000019147186298E0/
  DATA bipcs(9)/.00000006453106284E0/
  DATA bipcs(10)/ - .00000000784485467E0/
  DATA bipcs(11)/ - .00000000096077216E0/
  DATA bipcs(12)/.00000000070004713E0/
  DATA bipcs(13)/ - .00000000017731789E0/
  DATA bipcs(14)/.00000000002272089E0/
  DATA bipcs(15)/.00000000000165404E0/
  DATA bipcs(16)/ - .00000000000185171E0/
  DATA bipcs(17)/.00000000000059576E0/
  DATA bipcs(18)/ - .00000000000012194E0/
  DATA bipcs(19)/.00000000000001334E0/
  DATA bipcs(20)/.00000000000000172E0/
  DATA bipcs(21)/ - .00000000000000145E0/
  DATA bipcs(22)/.00000000000000049E0/
  DATA bipcs(23)/ - .00000000000000011E0/
  DATA bipcs(24)/.00000000000000001E0/
  DATA bip2cs(1)/ - .113596737585988679E0/
  DATA bip2cs(2)/.0041381473947881595E0/
  DATA bip2cs(3)/.0001353470622119332E0/
  DATA bip2cs(4)/.0000104273166530153E0/
  DATA bip2cs(5)/.0000013474954767849E0/
  DATA bip2cs(6)/.0000001696537405438E0/
  DATA bip2cs(7)/ - .0000000100965008656E0/
  DATA bip2cs(8)/ - .0000000167291194937E0/
  DATA bip2cs(9)/ - .0000000045815364485E0/
  DATA bip2cs(10)/.0000000003736681366E0/
  DATA bip2cs(11)/.0000000005766930320E0/
  DATA bip2cs(12)/.0000000000621812650E0/
  DATA bip2cs(13)/ - .0000000000632941202E0/
  DATA bip2cs(14)/ - .0000000000149150479E0/
  DATA bip2cs(15)/.0000000000078896213E0/
  DATA bip2cs(16)/.0000000000024960513E0/
  DATA bip2cs(17)/ - .0000000000012130075E0/
  DATA bip2cs(18)/ - .0000000000003740493E0/
  DATA bip2cs(19)/.0000000000002237727E0/
  DATA bip2cs(20)/.0000000000000474902E0/
  DATA bip2cs(21)/ - .0000000000000452616E0/
  DATA bip2cs(22)/ - .0000000000000030172E0/
  DATA bip2cs(23)/.0000000000000091058E0/
  DATA bip2cs(24)/ - .0000000000000009814E0/
  DATA bip2cs(25)/ - .0000000000000016429E0/
  DATA bip2cs(26)/.0000000000000005533E0/
  DATA bip2cs(27)/.0000000000000002175E0/
  DATA bip2cs(28)/ - .0000000000000001737E0/
  DATA bip2cs(29)/ - .0000000000000000010E0/
  DATA atr/8.7506905708484345E0/
  DATA btr/ - 2.093836321356054E0/
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  BIE
  IF ( first ) THEN
    eta = 0.1*R1MACH(3)
    nbif = INITS(bifcs,9,eta)
    nbig = INITS(bigcs,8,eta)
    nbif2 = INITS(bif2cs,10,eta)
    nbig2 = INITS(big2cs,10,eta)
    nbip = INITS(bipcs,24,eta)
    nbip2 = INITS(bip2cs,29,eta)
    !
    x3sml = eta**0.3333
    x32sml = 1.3104*x3sml**2
    xbig = R1MACH(2)**0.6666
  ENDIF
  first = .FALSE.
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
    GOTO 99999
  ENDIF
  sqrtx = SQRT(X)
  z = atr/(X*sqrtx) + btr
  BIE = (0.625+CSEVL(z,bipcs,nbip))/SQRT(sqrtx)
  RETURN
  !
  99999 CONTINUE
  END FUNCTION BIE
