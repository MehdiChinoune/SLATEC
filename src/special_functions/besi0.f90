!** BESI0
REAL(SP) FUNCTION BESI0(X)
  !> Compute the hyperbolic Bessel function of the first kind
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
  USE service, ONLY : XERMSG, R1MACH
  REAL(SP) :: X
  REAL(SP) :: y
  INTEGER, SAVE :: nti0
  REAL(SP), PARAMETER :: xsml = SQRT(4.5_SP*R1MACH(3)), xmax = LOG(R1MACH(2))
  REAL(SP), PARAMETER :: bi0cs(12) = [ -.07660547252839144951_SP, 1.927337953993808270_SP, &
    .2282644586920301339_SP, .01304891466707290428_SP, .00043442709008164874_SP, &
    .00000942265768600193_SP, .00000014340062895106_SP, .00000000161384906966_SP, &
    .00000000001396650044_SP, .00000000000009579451_SP, .00000000000000053339_SP, &
    .00000000000000000245_SP ]
  LOGICAL, SAVE :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  BESI0
  IF( first ) THEN
    nti0 = INITS(bi0cs,12,0.1_SP*R1MACH(3))
    first = .FALSE.
  END IF
  !
  y = ABS(X)
  IF( y>3._SP ) THEN
    !
    IF( y>xmax ) CALL XERMSG('BESI0','ABS(X) SO BIG I0 OVERFLOWS',1,2)
    !
    BESI0 = EXP(y)*BESI0E(X)
    RETURN
  END IF
  !
  BESI0 = 1._SP
  IF( y>xsml ) BESI0 = 2.75_SP + CSEVL(y*y/4.5_SP-1._SP,bi0cs,nti0)
  RETURN
END FUNCTION BESI0
