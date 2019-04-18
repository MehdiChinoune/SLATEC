!** BESI1
REAL FUNCTION BESI1(X)
  !>
  !***
  !  Compute the modified (hyperbolic) Bessel function of the
  !            first kind of order one.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C10B1
  !***
  ! **Type:**      SINGLE PRECISION (BESI1-S, DBESI1-D)
  !***
  ! **Keywords:**  FIRST KIND, FNLIB, HYPERBOLIC BESSEL FUNCTION,
  !             MODIFIED BESSEL FUNCTION, ORDER ONE, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! BESI1(X) calculates the modified (hyperbolic) Bessel function
  ! of the first kind of order one for real argument X.
  !
  ! Series for BI1        on the interval  0.          to  9.00000D+00
  !                                        with weighted error   2.40E-17
  !                                         log weighted error  16.62
  !                               significant figures required  16.23
  !                                    decimal places required  17.14
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  BESI1E, CSEVL, INITS, R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  USE service, ONLY : XERMSG, R1MACH
  REAL X, y
  INTEGER, SAVE :: nti1
  REAL, SAVE :: xmin, xsml, xmax
  REAL, PARAMETER :: bi1cs(11) = [ -.001971713261099859E0, .40734887667546481E0, &
    .034838994299959456E0, .001545394556300123E0, .000041888521098377E0, &
    .000000764902676483E0, .000000010042493924E0, .000000000099322077E0, &
    .000000000000766380E0, .000000000000004741E0, .000000000000000024E0 ]
  LOGICAL :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  BESI1
  IF ( first ) THEN
    nti1 = INITS(bi1cs,11,0.1*R1MACH(3))
    xmin = 2.0*R1MACH(1)
    xsml = SQRT(4.5*R1MACH(3))
    xmax = LOG(R1MACH(2))
    first = .FALSE.
  END IF
  !
  y = ABS(X)
  IF ( y>3.0 ) THEN
    !
    IF ( y>xmax ) CALL XERMSG('SLATEC','BESI1','ABS(X) SO BIG I1 OVERFLOWS',&
      2,2)
    !
    BESI1 = EXP(y)*BESI1E(X)
    RETURN
  END IF
  !
  BESI1 = 0.0
  IF ( y==0.0 ) RETURN
  !
  IF ( y<=xmin ) CALL XERMSG('SLATEC','BESI1',&
    'ABS(X) SO SMALL I1 UNDERFLOWS',1,1)
  IF ( y>xmin ) BESI1 = 0.5*X
  IF ( y>xsml ) BESI1 = X*(.875+CSEVL(y*y/4.5-1.,bi1cs,nti1))
  RETURN
END FUNCTION BESI1
