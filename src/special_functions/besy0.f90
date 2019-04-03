!** BESY0
REAL FUNCTION BESY0(X)
  IMPLICIT NONE
  !>
  !***
  !  Compute the Bessel function of the second kind of order
  !            zero.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C10A1
  !***
  ! **Type:**      SINGLE PRECISION (BESY0-S, DBESY0-D)
  !***
  ! **Keywords:**  BESSEL FUNCTION, FNLIB, ORDER ZERO, SECOND KIND,
  !             SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! BESY0(X) calculates the Bessel function of the second kind
  ! of order zero for real argument X.
  !
  ! Series for BY0        on the interval  0.          to  1.60000D+01
  !                                        with weighted error   1.20E-17
  !                                         log weighted error  16.92
  !                               significant figures required  16.15
  !                                    decimal places required  17.48
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
  ! **Routines called:**  BESJ0, CSEVL, INITS, R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)

  REAL ampl, BESJ0, CSEVL, R1MACH, theta, X, y, z
  INTEGER INITS
  INTEGER, SAVE :: nty0, ntm0, ntth0
  REAL, SAVE :: xsml, xmax
  REAL, PARAMETER :: by0cs(13) = [ -.011277839392865573E0, -.12834523756042035E0, &
    -.10437884799794249E0,  .023662749183969695E0, -.002090391647700486E0, &
    .000103975453939057E0, -.000003369747162423E0, .000000077293842676E0, &
    -.000000001324976772E0, .000000000017648232E0, -.000000000000188105E0, &
    .000000000000001641E0, -.000000000000000011E0 ]
  REAL, PARAMETER :: bm0cs(21) = [ .09284961637381644E0, -.00142987707403484E0, &
    .00002830579271257E0, -.00000143300611424E0, .00000012028628046E0, &
    -.00000001397113013E0, .00000000204076188E0, -.00000000035399669E0, &
    .00000000007024759E0, -.00000000001554107E0, .00000000000376226E0, &
    -.00000000000098282E0, .00000000000027408E0, -.00000000000008091E0, &
    .00000000000002511E0, -.00000000000000814E0, .00000000000000275E0, &
    -.00000000000000096E0, .00000000000000034E0, -.00000000000000012E0, &
    .00000000000000004E0 ]
  REAL, PARAMETER :: bth0cs(24) = [ -.24639163774300119E0, .001737098307508963E0, &
    -.000062183633402968E0, .000004368050165742E0, -.000000456093019869E0, &
    .000000062197400101E0, -.000000010300442889E0, .000000001979526776E0, &
    -.000000000428198396E0, .000000000102035840E0, -.000000000026363898E0, &
    .000000000007297935E0, -.000000000002144188E0, .000000000000663693E0, &
    -.000000000000215126E0, .000000000000072659E0, -.000000000000025465E0, &
    .000000000000009229E0, -.000000000000003448E0, .000000000000001325E0, &
    -.000000000000000522E0, .000000000000000210E0, -.000000000000000087E0, &
    .000000000000000036E0 ]
  REAL, PARAMETER :: twodpi = 0.63661977236758134E0
  REAL, PARAMETER ::  pi4 = 0.78539816339744831E0
  LOGICAL :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  BESY0
  IF ( first ) THEN
    nty0 = INITS(by0cs,13,0.1*R1MACH(3))
    ntm0 = INITS(bm0cs,21,0.1*R1MACH(3))
    ntth0 = INITS(bth0cs,24,0.1*R1MACH(3))
    !
    xsml = SQRT(4.0*R1MACH(3))
    xmax = 1.0/R1MACH(4)
    first = .FALSE.
  ENDIF
  !
  IF ( X<=0. ) CALL XERMSG('SLATEC','BESY0','X IS ZERO OR NEGATIVE',1,2)
  IF ( X>4.0 ) THEN
    !
    IF ( X>xmax ) CALL XERMSG('SLATEC','BESY0',&
      'NO PRECISION BECAUSE X IS BIG',2,2)
    !
    z = 32.0/X**2 - 1.0
    ampl = (0.75+CSEVL(z,bm0cs,ntm0))/SQRT(X)
    theta = X - pi4 + CSEVL(z,bth0cs,ntth0)/X
    BESY0 = ampl*SIN(theta)
    RETURN
  ENDIF
  !
  y = 0.
  IF ( X>xsml ) y = X*X
  BESY0 = twodpi*LOG(0.5*X)*BESJ0(X) + .375 + CSEVL(.125*y-1.,by0cs,nty0)
  RETURN
END FUNCTION BESY0
