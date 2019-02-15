!DECK AI
FUNCTION AI(X)
  IMPLICIT NONE
  REAL AI, AIE, aifcs, aigcs, CSEVL, R1MACH, theta, X, x3sml, xm, &
    xmax, xmaxt, z
  INTEGER INITS, naif, naig
  !***BEGIN PROLOGUE  AI
  !***PURPOSE  Evaluate the Airy function.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C10D
  !***TYPE      SINGLE PRECISION (AI-S, DAI-D)
  !***KEYWORDS  AIRY FUNCTION, FNLIB, SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! AI(X) computes the Airy function Ai(X)
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
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  AIE, CSEVL, INITS, R1MACH, R9AIMP, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920618  Removed space from variable names.  (RWC, WRB)
  !***END PROLOGUE  AI
  DIMENSION aifcs(9), aigcs(8)
  LOGICAL first
  SAVE aifcs, aigcs, naif, naig, x3sml, xmax, first
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
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  AI
  IF ( first ) THEN
    naif = INITS(aifcs,9,0.1*R1MACH(3))
    naig = INITS(aigcs,8,0.1*R1MACH(3))
    !
    x3sml = R1MACH(3)**0.3334
    xmaxt = (-1.5*LOG(R1MACH(1)))**0.6667
    xmax = xmaxt - xmaxt*LOG(xmaxt)/(4.0*SQRT(xmaxt)+1.0) - 0.01
  ENDIF
  first = .FALSE.
  !
  IF ( X<(-1.0) ) THEN
    CALL R9AIMP(X,xm,theta)
    AI = xm*COS(theta)
    RETURN
    !
  ELSEIF ( X<=1.0 ) THEN
    z = 0.0
    IF ( ABS(X)>x3sml ) z = X**3
    AI = 0.375 + (CSEVL(z,aifcs,naif)-X*(0.25+CSEVL(z,aigcs,naig)))
    RETURN
    !
  ELSEIF ( X>xmax ) THEN
    !
    AI = 0.0
    CALL XERMSG('SLATEC','AI','X SO BIG AI UNDERFLOWS',1,1)
    GOTO 99999
  ENDIF
  AI = AIE(X)*EXP(-2.0*X*SQRT(X)/3.0)
  RETURN
  !
  99999 CONTINUE
  END FUNCTION AI
