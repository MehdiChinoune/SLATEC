!** BESJ0
REAL FUNCTION BESJ0(X)
  IMPLICIT NONE
  !>
  !***
  !  Compute the Bessel function of the first kind of order
  !            zero.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C10A1
  !***
  ! **Type:**      SINGLE PRECISION (BESJ0-S, DBESJ0-D)
  !***
  ! **Keywords:**  BESSEL FUNCTION, FIRST KIND, FNLIB, ORDER ZERO,
  !             SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! BESJ0(X) calculates the Bessel function of the first kind of
  ! order zero for real argument X.
  !
  ! Series for BJ0        on the interval  0.          to  1.60000D+01
  !                                        with weighted error   7.47E-18
  !                                         log weighted error  17.13
  !                               significant figures required  16.98
  !                                    decimal places required  17.68
  !
  ! Series for BM0        on the interval  0.          to  6.25000D-02
  !                                        with weighted error   4.98E-17
  !                                         log weighted error  16.30
  !                               significant figures required  14.97
  !                                    decimal places required  16.96
  !
  ! Series for BTH0       on the interval  0.          to  6.25000D-02
  !                                        with weighted error   3.67E-17
  !                                         log weighted error  16.44
  !                               significant figures required  15.53
  !                                    decimal places required  17.13
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  CSEVL, INITS, R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   890210  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)

  REAL ampl, CSEVL, R1MACH, theta, X, y, z
  INTEGER INITS
  INTEGER, SAVE :: ntj0, ntm0, ntth0
  REAL, SAVE :: xsml, xmax
  REAL, PARAMETER :: bj0cs(13) = [ .100254161968939137E0, -.665223007764405132E0, &
    .248983703498281314E0,  -.0332527231700357697E0, .0023114179304694015E0, &
    -.0000991127741995080E0, .0000028916708643998E0,-.0000000612108586630E0, &
    .0000000009838650793E0, -.0000000000124235515E0, .0000000000001265433E0, &
    -.0000000000000010619E0, .0000000000000000074E0 ]
  REAL, PARAMETER :: bm0cs(21) = [ .09284961637381644E0,-.00142987707403484E0, &
    .00002830579271257E0, -.00000143300611424E0, .00000012028628046E0, &
    -.00000001397113013E0, .00000000204076188E0,-.00000000035399669E0, &
    .00000000007024759E0, -.00000000001554107E0, .00000000000376226E0, &
    -.00000000000098282E0, .00000000000027408E0,-.00000000000008091E0, &
    .00000000000002511E0, -.00000000000000814E0, .00000000000000275E0, &
    -.00000000000000096E0, .00000000000000034E0,-.00000000000000012E0, &
    .00000000000000004E0 ]
  REAL, PARAMETER :: bth0cs(24) = [ -.24639163774300119E0, .001737098307508963E0, &
    -.000062183633402968E0, .000004368050165742E0,-.000000456093019869E0, &
    .000000062197400101E0, -.000000010300442889E0, .000000001979526776E0, &
    -.000000000428198396E0, .000000000102035840E0,-.000000000026363898E0, &
    .000000000007297935E0, -.000000000002144188E0, .000000000000663693E0, &
    -.000000000000215126E0, .000000000000072659E0,-.000000000000025465E0, &
    .000000000000009229E0, -.000000000000003448E0, .000000000000001325E0, &
    -.000000000000000522E0, .000000000000000210E0,-.000000000000000087E0, &
    .000000000000000036E0 ]
  REAL, PARAMETER ::  pi4 = 0.78539816339744831E0
  LOGICAL :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  BESJ0
  IF ( first ) THEN
    ntj0 = INITS(bj0cs,13,0.1*R1MACH(3))
    ntm0 = INITS(bm0cs,21,0.1*R1MACH(3))
    ntth0 = INITS(bth0cs,24,0.1*R1MACH(3))
    !
    xsml = SQRT(8.0*R1MACH(3))
    xmax = 1.0/R1MACH(4)
    first = .FALSE.
  ENDIF
  !
  y = ABS(X)
  IF ( y>4.0 ) THEN
    !
    IF ( y>xmax ) CALL XERMSG('SLATEC','BESJ0',&
      'NO PRECISION BECAUSE ABS(X) IS TOO BIG',1,2)
    !
    z = 32.0/y**2 - 1.0
    ampl = (0.75+CSEVL(z,bm0cs,ntm0))/SQRT(y)
    theta = y - pi4 + CSEVL(z,bth0cs,ntth0)/y
    BESJ0 = ampl*COS(theta)
    RETURN
  ENDIF
  !
  BESJ0 = 1.0
  IF ( y>xsml ) BESJ0 = CSEVL(.125*y*y-1.,bj0cs,ntj0)
  RETURN
END FUNCTION BESJ0
