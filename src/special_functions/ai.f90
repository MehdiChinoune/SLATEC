!** AI
REAL FUNCTION AI(X)
  !>
  !  Evaluate the Airy function.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C10D
  !***
  ! **Type:**      SINGLE PRECISION (AI-S, DAI-D)
  !***
  ! **Keywords:**  AIRY FUNCTION, FNLIB, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
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
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  AIE, CSEVL, INITS, R1MACH, R9AIMP, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   920618  Removed space from variable names.  (RWC, WRB)
  USE service, ONLY : XERMSG, R1MACH
  REAL theta, X, xm, xmaxt, z
  INTEGER, SAVE :: naif, naig
  REAL, SAVE :: x3sml, xmax
  REAL, PARAMETER :: aifcs(9) = [ -.03797135849666999750E0, .05919188853726363857E0, &
    .00098629280577279975E0, .00000684884381907656E0, .00000002594202596219E0, &
    .00000000006176612774E0, .00000000000010092454E0, .00000000000000012014E0, &
    .00000000000000000010E0 ]
  REAL, PARAMETER :: aigcs(8) = [ .01815236558116127E0, .02157256316601076E0, &
    .00025678356987483E0, .00000142652141197E0, .00000000457211492E0, &
    .00000000000952517E0, .00000000000001392E0, .00000000000000001E0 ]
  LOGICAL :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  AI
  IF ( first ) THEN
    naif = INITS(aifcs,9,0.1*R1MACH(3))
    naig = INITS(aigcs,8,0.1*R1MACH(3))
    !
    x3sml = R1MACH(3)**0.3334
    xmaxt = (-1.5*LOG(R1MACH(1)))**0.6667
    xmax = xmaxt - xmaxt*LOG(xmaxt)/(4.0*SQRT(xmaxt)+1.0) - 0.01
    first = .FALSE.
  END IF
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
    RETURN
  END IF
  AI = AIE(X)*EXP(-2.0*X*SQRT(X)/3.0)
  RETURN
END FUNCTION AI
