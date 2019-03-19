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
  
  REAL ampl, bj0cs, bm0cs, bth0cs, CSEVL, pi4, R1MACH, &
    theta, X, xmax, xsml, y, z
  INTEGER INITS, ntj0, ntm0, ntth0
  DIMENSION bj0cs(13), bm0cs(21), bth0cs(24)
  LOGICAL first
  SAVE bj0cs, bm0cs, bth0cs, pi4, ntj0, ntm0, ntth0, xsml, xmax, &
    first
  DATA bj0cs(1)/.100254161968939137E0/
  DATA bj0cs(2)/ - .665223007764405132E0/
  DATA bj0cs(3)/.248983703498281314E0/
  DATA bj0cs(4)/ - .0332527231700357697E0/
  DATA bj0cs(5)/.0023114179304694015E0/
  DATA bj0cs(6)/ - .0000991127741995080E0/
  DATA bj0cs(7)/.0000028916708643998E0/
  DATA bj0cs(8)/ - .0000000612108586630E0/
  DATA bj0cs(9)/.0000000009838650793E0/
  DATA bj0cs(10)/ - .0000000000124235515E0/
  DATA bj0cs(11)/.0000000000001265433E0/
  DATA bj0cs(12)/ - .0000000000000010619E0/
  DATA bj0cs(13)/.0000000000000000074E0/
  DATA bm0cs(1)/.09284961637381644E0/
  DATA bm0cs(2)/ - .00142987707403484E0/
  DATA bm0cs(3)/.00002830579271257E0/
  DATA bm0cs(4)/ - .00000143300611424E0/
  DATA bm0cs(5)/.00000012028628046E0/
  DATA bm0cs(6)/ - .00000001397113013E0/
  DATA bm0cs(7)/.00000000204076188E0/
  DATA bm0cs(8)/ - .00000000035399669E0/
  DATA bm0cs(9)/.00000000007024759E0/
  DATA bm0cs(10)/ - .00000000001554107E0/
  DATA bm0cs(11)/.00000000000376226E0/
  DATA bm0cs(12)/ - .00000000000098282E0/
  DATA bm0cs(13)/.00000000000027408E0/
  DATA bm0cs(14)/ - .00000000000008091E0/
  DATA bm0cs(15)/.00000000000002511E0/
  DATA bm0cs(16)/ - .00000000000000814E0/
  DATA bm0cs(17)/.00000000000000275E0/
  DATA bm0cs(18)/ - .00000000000000096E0/
  DATA bm0cs(19)/.00000000000000034E0/
  DATA bm0cs(20)/ - .00000000000000012E0/
  DATA bm0cs(21)/.00000000000000004E0/
  DATA bth0cs(1)/ - .24639163774300119E0/
  DATA bth0cs(2)/.001737098307508963E0/
  DATA bth0cs(3)/ - .000062183633402968E0/
  DATA bth0cs(4)/.000004368050165742E0/
  DATA bth0cs(5)/ - .000000456093019869E0/
  DATA bth0cs(6)/.000000062197400101E0/
  DATA bth0cs(7)/ - .000000010300442889E0/
  DATA bth0cs(8)/.000000001979526776E0/
  DATA bth0cs(9)/ - .000000000428198396E0/
  DATA bth0cs(10)/.000000000102035840E0/
  DATA bth0cs(11)/ - .000000000026363898E0/
  DATA bth0cs(12)/.000000000007297935E0/
  DATA bth0cs(13)/ - .000000000002144188E0/
  DATA bth0cs(14)/.000000000000663693E0/
  DATA bth0cs(15)/ - .000000000000215126E0/
  DATA bth0cs(16)/.000000000000072659E0/
  DATA bth0cs(17)/ - .000000000000025465E0/
  DATA bth0cs(18)/.000000000000009229E0/
  DATA bth0cs(19)/ - .000000000000003448E0/
  DATA bth0cs(20)/.000000000000001325E0/
  DATA bth0cs(21)/ - .000000000000000522E0/
  DATA bth0cs(22)/.000000000000000210E0/
  DATA bth0cs(23)/ - .000000000000000087E0/
  DATA bth0cs(24)/.000000000000000036E0/
  DATA pi4/0.78539816339744831E0/
  DATA first/.TRUE./
  !* FIRST EXECUTABLE STATEMENT  BESJ0
  IF ( first ) THEN
    ntj0 = INITS(bj0cs,13,0.1*R1MACH(3))
    ntm0 = INITS(bm0cs,21,0.1*R1MACH(3))
    ntth0 = INITS(bth0cs,24,0.1*R1MACH(3))
    !
    xsml = SQRT(8.0*R1MACH(3))
    xmax = 1.0/R1MACH(4)
  ENDIF
  first = .FALSE.
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
