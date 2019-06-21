!** AIE
REAL(SP) FUNCTION AIE(X)
  !> Calculate the Airy function for a negative argument and an
  !            exponentially scaled Airy function for a non-negative
  !            argument.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C10D
  !***
  ! **Type:**      SINGLE PRECISION (AIE-S, DAIE-D)
  !***
  ! **Keywords:**  EXPONENTIALLY SCALED AIRY FUNCTION, FNLIB,
  !             SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! AIE(X) computes the exponentially scaled Airy function for
  ! non-negative X.  It evaluates AI(X) for X <= 0.0 and
  ! EXP(ZETA)*AI(X) for X >= 0.0 where ZETA = (2.0/3.0)*(X**1.5).
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
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  CSEVL, INITS, R1MACH, R9AIMP

  !* REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890206  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   920618  Removed space from variable names.  (RWC, WRB)
  USE service, ONLY : R1MACH
  REAL(SP) :: X
  REAL(SP) :: sqrtx, theta, xm, z
  INTEGER, SAVE :: naif, naig, naip
  REAL(SP), PARAMETER :: eta = 0.1_SP*R1MACH(3), x3sml = eta**0.3333_SP, &
    x32sml = 1.3104_SP*x3sml**2, xbig = R1MACH(2)**0.6666_SP
  REAL(SP), PARAMETER :: aifcs(9) = [ -.03797135849666999750_SP, .05919188853726363857_SP, &
    .00098629280577279975_SP, .00000684884381907656_SP, .00000002594202596219_SP, &
    .00000000006176612774_SP, .00000000000010092454_SP, .00000000000000012014_SP, &
    .00000000000000000010_SP ]
  REAL(SP), PARAMETER :: aigcs(8) = [ .01815236558116127_SP, .02157256316601076_SP, &
    .00025678356987483_SP, .00000142652141197_SP, .00000000457211492_SP, &
    .00000000000952517_SP, .00000000000001392_SP, .00000000000000001_SP ]
  REAL(SP), PARAMETER :: aipcs(34) = [ -.0187519297793868_SP, -.0091443848250055_SP, &
    .0009010457337825_SP, -.0001394184127221_SP, .0000273815815785_SP, &
    -.0000062750421119_SP, .0000016064844184_SP,-.0000004476392158_SP, &
    .0000001334635874_SP, -.0000000420735334_SP, .0000000139021990_SP, &
    -.0000000047831848_SP, .0000000017047897_SP,-.0000000006268389_SP, &
    .0000000002369824_SP, -.0000000000918641_SP, .0000000000364278_SP, &
    -.0000000000147475_SP, .0000000000060851_SP,-.0000000000025552_SP, &
    .0000000000010906_SP, -.0000000000004725_SP, .0000000000002076_SP, &
    -.0000000000000924_SP, .0000000000000417_SP,-.0000000000000190_SP, &
    .0000000000000087_SP, -.0000000000000040_SP, .0000000000000019_SP, &
    -.0000000000000009_SP, .0000000000000004_SP,-.0000000000000002_SP, &
    .0000000000000001_SP, -.0000000000000000_SP ]
  LOGICAL, SAVE :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  AIE
  IF( first ) THEN
    naif = INITS(aifcs,9,eta)
    naig = INITS(aigcs,8,eta)
    naip = INITS(aipcs,34,eta)
    first = .FALSE.
  END IF
  !
  IF( X<(-1._SP) ) THEN
    CALL R9AIMP(X,xm,theta)
    AIE = xm*COS(theta)
    RETURN
    !
  ELSEIF( X>1._SP ) THEN
    !
    sqrtx = SQRT(X)
    z = -1._SP
    IF( X<xbig ) z = 2._SP/(X*sqrtx) - 1._SP
    AIE = (.28125_SP+CSEVL(z,aipcs,naip))/SQRT(sqrtx)
    RETURN
  END IF
  z = 0._SP
  IF( ABS(X)>x3sml ) z = X**3
  AIE = 0.375_SP + (CSEVL(z,aifcs,naif)-X*(0.25_SP+CSEVL(z,aigcs,naig)))
  IF( X>x32sml ) AIE = AIE*EXP(2._SP*X*SQRT(X)/3._SP)
  RETURN
END FUNCTION AIE
