!DECK AIE
FUNCTION AIE(X)
  IMPLICIT NONE
  REAL AIE, aifcs, aigcs, aipcs, CSEVL, eta, R1MACH, sqrtx, theta, &
    X, x32sml, x3sml, xbig, xm, z
  INTEGER INITS, naif, naig, naip
  !***BEGIN PROLOGUE  AIE
  !***PURPOSE  Calculate the Airy function for a negative argument and an
  !            exponentially scaled Airy function for a non-negative
  !            argument.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C10D
  !***TYPE      SINGLE PRECISION (AIE-S, DAIE-D)
  !***KEYWORDS  EXPONENTIALLY SCALED AIRY FUNCTION, FNLIB,
  !             SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! AIE(X) computes the exponentially scaled Airy function for
  ! non-negative X.  It evaluates AI(X) for X .LE. 0.0 and
  ! EXP(ZETA)*AI(X) for X .GE. 0.0 where ZETA = (2.0/3.0)*(X**1.5).
  !
  ! Series for AIF        on the interval -1.00000D+00 to  1.00000D+00
  !                                        with weighted error   1.09E-19
  !                                         log weighted error  18.96
  !                               significant figures required  17.76
  !                                    decimal places required  19.44
  !
  ! Series for AIG        on the interval -1.00000D+00 to  1.00000D+00
  !                                        with weighted error   1.51E-17
  !                                         log weighted error  16.82
  !                               significant figures required  15.19
  !                                    decimal places required  17.27
  !
  ! Series for AIP        on the interval  0.          to  1.00000D+00
  !                                        with weighted error   5.10E-17
  !                                         log weighted error  16.29
  !                               significant figures required  14.41
  !                                    decimal places required  17.06
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  CSEVL, INITS, R1MACH, R9AIMP
  !***REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890206  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920618  Removed space from variable names.  (RWC, WRB)
  !***END PROLOGUE  AIE
  DIMENSION aifcs(9), aigcs(8), aipcs(34)
  LOGICAL first
  SAVE aifcs, aigcs, aipcs, naif, naig, naip, x3sml, x32sml, xbig, &
    first
  DATA aifcs(1)/ - .03797135849666999750E0/
  DATA aifcs(2)/.05919188853726363857E0/
  DATA aifcs(3)/.00098629280577279975E0/
  DATA aifcs(4)/.00000684884381907656E0/
  DATA aifcs(5)/.00000002594202596219E0/
  DATA aifcs(6)/.00000000006176612774E0/
  DATA aifcs(7)/.00000000000010092454E0/
  DATA aifcs(8)/.00000000000000012014E0/
  DATA aifcs(9)/.00000000000000000010E0/
  DATA aigcs(1)/.01815236558116127E0/
  DATA aigcs(2)/.02157256316601076E0/
  DATA aigcs(3)/.00025678356987483E0/
  DATA aigcs(4)/.00000142652141197E0/
  DATA aigcs(5)/.00000000457211492E0/
  DATA aigcs(6)/.00000000000952517E0/
  DATA aigcs(7)/.00000000000001392E0/
  DATA aigcs(8)/.00000000000000001E0/
  DATA aipcs(1)/ - .0187519297793868E0/
  DATA aipcs(2)/ - .0091443848250055E0/
  DATA aipcs(3)/.0009010457337825E0/
  DATA aipcs(4)/ - .0001394184127221E0/
  DATA aipcs(5)/.0000273815815785E0/
  DATA aipcs(6)/ - .0000062750421119E0/
  DATA aipcs(7)/.0000016064844184E0/
  DATA aipcs(8)/ - .0000004476392158E0/
  DATA aipcs(9)/.0000001334635874E0/
  DATA aipcs(10)/ - .0000000420735334E0/
  DATA aipcs(11)/.0000000139021990E0/
  DATA aipcs(12)/ - .0000000047831848E0/
  DATA aipcs(13)/.0000000017047897E0/
  DATA aipcs(14)/ - .0000000006268389E0/
  DATA aipcs(15)/.0000000002369824E0/
  DATA aipcs(16)/ - .0000000000918641E0/
  DATA aipcs(17)/.0000000000364278E0/
  DATA aipcs(18)/ - .0000000000147475E0/
  DATA aipcs(19)/.0000000000060851E0/
  DATA aipcs(20)/ - .0000000000025552E0/
  DATA aipcs(21)/.0000000000010906E0/
  DATA aipcs(22)/ - .0000000000004725E0/
  DATA aipcs(23)/.0000000000002076E0/
  DATA aipcs(24)/ - .0000000000000924E0/
  DATA aipcs(25)/.0000000000000417E0/
  DATA aipcs(26)/ - .0000000000000190E0/
  DATA aipcs(27)/.0000000000000087E0/
  DATA aipcs(28)/ - .0000000000000040E0/
  DATA aipcs(29)/.0000000000000019E0/
  DATA aipcs(30)/ - .0000000000000009E0/
  DATA aipcs(31)/.0000000000000004E0/
  DATA aipcs(32)/ - .0000000000000002E0/
  DATA aipcs(33)/.0000000000000001E0/
  DATA aipcs(34)/ - .0000000000000000E0/
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  AIE
  IF ( first ) THEN
    eta = 0.1*R1MACH(3)
    naif = INITS(aifcs,9,eta)
    naig = INITS(aigcs,8,eta)
    naip = INITS(aipcs,34,eta)
    !
    x3sml = eta**0.3333
    x32sml = 1.3104*x3sml**2
    xbig = R1MACH(2)**0.6666
  ENDIF
  first = .FALSE.
  !
  IF ( X<(-1.0) ) THEN
    CALL R9AIMP(X,xm,theta)
    AIE = xm*COS(theta)
    RETURN
    !
  ELSEIF ( X>1.0 ) THEN
    !
    sqrtx = SQRT(X)
    z = -1.0
    IF ( X<xbig ) z = 2.0/(X*sqrtx) - 1.0
    AIE = (.28125+CSEVL(z,aipcs,naip))/SQRT(sqrtx)
    GOTO 99999
  ENDIF
  z = 0.0
  IF ( ABS(X)>x3sml ) z = X**3
  AIE = 0.375 + (CSEVL(z,aifcs,naif)-X*(0.25+CSEVL(z,aigcs,naig)))
  IF ( X>x32sml ) AIE = AIE*EXP(2.0*X*SQRT(X)/3.0)
  RETURN
  !
  99999 CONTINUE
  END FUNCTION AIE
