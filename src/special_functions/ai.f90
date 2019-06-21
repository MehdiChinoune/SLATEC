!** AI
REAL(SP) FUNCTION AI(X)
  !> Evaluate the Airy function.
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
  REAL(SP) :: X
  REAL(SP) :: theta, xm, z
  INTEGER, SAVE :: naif, naig
  REAL(SP), PARAMETER :: x3sml = R1MACH(3)**0.3334_SP, &
    xmaxt = (-1.5_SP*LOG(R1MACH(1)))**0.6667_SP, &
    xmax = xmaxt - xmaxt*LOG(xmaxt)/(4._SP*SQRT(xmaxt)+1._SP) - 0.01_SP
  REAL(SP), PARAMETER :: aifcs(9) = [ -.03797135849666999750_SP, .05919188853726363857_SP, &
    .00098629280577279975_SP, .00000684884381907656_SP, .00000002594202596219_SP, &
    .00000000006176612774_SP, .00000000000010092454_SP, .00000000000000012014_SP, &
    .00000000000000000010_SP ]
  REAL(SP), PARAMETER :: aigcs(8) = [ .01815236558116127_SP, .02157256316601076_SP, &
    .00025678356987483_SP, .00000142652141197_SP, .00000000457211492_SP, &
    .00000000000952517_SP, .00000000000001392_SP, .00000000000000001_SP ]
  LOGICAL, SAVE :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  AI
  IF( first ) THEN
    naif = INITS(aifcs,9,0.1_SP*R1MACH(3))
    naig = INITS(aigcs,8,0.1_SP*R1MACH(3))
    first = .FALSE.
  END IF
  !
  IF( X<(-1._SP) ) THEN
    CALL R9AIMP(X,xm,theta)
    AI = xm*COS(theta)
    RETURN
    !
  ELSEIF( X<=1._SP ) THEN
    z = 0._SP
    IF( ABS(X)>x3sml ) z = X**3
    AI = 0.375_SP + (CSEVL(z,aifcs,naif)-X*(0.25_SP+CSEVL(z,aigcs,naig)))
    RETURN
    !
  ELSEIF( X>xmax ) THEN
    !
    AI = 0._SP
    CALL XERMSG('AI','X SO BIG AI UNDERFLOWS',1,1)
    RETURN
  END IF
  AI = AIE(X)*EXP(-2._SP*X*SQRT(X)/3._SP)
  RETURN
END FUNCTION AI
