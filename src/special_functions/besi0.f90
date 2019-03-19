!** BESI0
REAL FUNCTION BESI0(X)
  IMPLICIT NONE
  !>
  !***
  !  Compute the hyperbolic Bessel function of the first kind
  !            of order zero.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C10B1
  !***
  ! **Type:**      SINGLE PRECISION (BESI0-S, DBESI0-D)
  !***
  ! **Keywords:**  FIRST KIND, FNLIB, HYPERBOLIC BESSEL FUNCTION,
  !             MODIFIED BESSEL FUNCTION, ORDER ZERO, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! BESI0(X) computes the modified (hyperbolic) Bessel function
  ! of the first kind of order zero and real argument X.
  !
  ! Series for BI0        on the interval  0.          to  9.00000D+00
  !                                        with weighted error   2.46E-18
  !                                         log weighted error  17.61
  !                               significant figures required  17.90
  !                                    decimal places required  18.15
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  BESI0E, CSEVL, INITS, R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  
  REAL BESI0E, bi0cs, CSEVL, R1MACH, X, xmax, xsml, y
  INTEGER INITS, nti0
  DIMENSION bi0cs(12)
  LOGICAL first
  SAVE bi0cs, nti0, xsml, xmax, first
  DATA bi0cs(1)/ - .07660547252839144951E0/
  DATA bi0cs(2)/1.927337953993808270E0/
  DATA bi0cs(3)/.2282644586920301339E0/
  DATA bi0cs(4)/.01304891466707290428E0/
  DATA bi0cs(5)/.00043442709008164874E0/
  DATA bi0cs(6)/.00000942265768600193E0/
  DATA bi0cs(7)/.00000014340062895106E0/
  DATA bi0cs(8)/.00000000161384906966E0/
  DATA bi0cs(9)/.00000000001396650044E0/
  DATA bi0cs(10)/.00000000000009579451E0/
  DATA bi0cs(11)/.00000000000000053339E0/
  DATA bi0cs(12)/.00000000000000000245E0/
  DATA first/.TRUE./
  !* FIRST EXECUTABLE STATEMENT  BESI0
  IF ( first ) THEN
    nti0 = INITS(bi0cs,12,0.1*R1MACH(3))
    xsml = SQRT(4.5*R1MACH(3))
    xmax = LOG(R1MACH(2))
  ENDIF
  first = .FALSE.
  !
  y = ABS(X)
  IF ( y>3.0 ) THEN
    !
    IF ( y>xmax ) CALL XERMSG('SLATEC','BESI0','ABS(X) SO BIG I0 OVERFLOWS',&
      1,2)
    !
    BESI0 = EXP(y)*BESI0E(X)
    RETURN
  ENDIF
  !
  BESI0 = 1.0
  IF ( y>xsml ) BESI0 = 2.75 + CSEVL(y*y/4.5-1.0,bi0cs,nti0)
  RETURN
END FUNCTION BESI0
