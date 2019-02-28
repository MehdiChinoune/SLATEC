!DECK BESY0
FUNCTION BESY0(X)
  IMPLICIT NONE
  REAL ampl, BESJ0, BESY0, bm0cs, bth0cs, by0cs, CSEVL, pi4, &
    R1MACH, theta, twodpi, X, xmax, xsml, y, z
  INTEGER INITS, ntm0, ntth0, nty0
  !***BEGIN PROLOGUE  BESY0
  !***PURPOSE  Compute the Bessel function of the second kind of order
  !            zero.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C10A1
  !***TYPE      SINGLE PRECISION (BESY0-S, DBESY0-D)
  !***KEYWORDS  BESSEL FUNCTION, FNLIB, ORDER ZERO, SECOND KIND,
  !             SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
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
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  BESJ0, CSEVL, INITS, R1MACH, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !***END PROLOGUE  BESY0
  DIMENSION by0cs(13), bm0cs(21), bth0cs(24)
  LOGICAL first
  SAVE by0cs, bm0cs, bth0cs, twodpi, pi4, nty0, ntm0, ntth0, xsml, &
    xmax, first
  DATA by0cs(1)/ - .011277839392865573E0/
  DATA by0cs(2)/ - .12834523756042035E0/
  DATA by0cs(3)/ - .10437884799794249E0/
  DATA by0cs(4)/.023662749183969695E0/
  DATA by0cs(5)/ - .002090391647700486E0/
  DATA by0cs(6)/.000103975453939057E0/
  DATA by0cs(7)/ - .000003369747162423E0/
  DATA by0cs(8)/.000000077293842676E0/
  DATA by0cs(9)/ - .000000001324976772E0/
  DATA by0cs(10)/.000000000017648232E0/
  DATA by0cs(11)/ - .000000000000188105E0/
  DATA by0cs(12)/.000000000000001641E0/
  DATA by0cs(13)/ - .000000000000000011E0/
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
  DATA twodpi/0.63661977236758134E0/
  DATA pi4/0.78539816339744831E0/
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  BESY0
  IF ( first ) THEN
    nty0 = INITS(by0cs,13,0.1*R1MACH(3))
    ntm0 = INITS(bm0cs,21,0.1*R1MACH(3))
    ntth0 = INITS(bth0cs,24,0.1*R1MACH(3))
    !
    xsml = SQRT(4.0*R1MACH(3))
    xmax = 1.0/R1MACH(4)
  ENDIF
  first = .FALSE.
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
