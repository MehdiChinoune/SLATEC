!** BESK0E
REAL FUNCTION BESK0E(X)
  IMPLICIT NONE
  !>
  !***
  !  Compute the exponentially scaled modified (hyperbolic)
  !            Bessel function of the third kind of order zero.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C10B1
  !***
  ! **Type:**      SINGLE PRECISION (BESK0E-S, DBSK0E-D)
  !***
  ! **Keywords:**  EXPONENTIALLY SCALED, FNLIB, HYPERBOLIC BESSEL FUNCTION,
  !             MODIFIED BESSEL FUNCTION, ORDER ZERO, SPECIAL FUNCTIONS,
  !             THIRD KIND
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! BESK0E(X) computes the exponentially scaled modified (hyperbolic)
  ! Bessel function of third kind of order zero for real argument
  ! X .GT. 0.0, i.e., EXP(X)*K0(X).
  !
  ! Series for BK0        on the interval  0.          to  4.00000D+00
  !                                        with weighted error   3.57E-19
  !                                         log weighted error  18.45
  !                               significant figures required  17.99
  !                                    decimal places required  18.97
  !
  ! Series for AK0        on the interval  1.25000D-01 to  5.00000D-01
  !                                        with weighted error   5.34E-17
  !                                         log weighted error  16.27
  !                               significant figures required  14.92
  !                                    decimal places required  16.89
  !
  ! Series for AK02       on the interval  0.          to  1.25000D-01
  !                                        with weighted error   2.34E-17
  !                                         log weighted error  16.63
  !                               significant figures required  14.67
  !                                    decimal places required  17.20
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  BESI0, CSEVL, INITS, R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)

  REAL BESI0, CSEVL, R1MACH, X, y
  INTEGER INITS
  INTEGER, SAVE :: ntk0, ntak0, ntak02
  REAL, SAVE :: xsml
  REAL, PARAMETER :: bk0cs(11) = [ -.03532739323390276872E0, .3442898999246284869E0, &
    .03597993651536150163E0, .00126461541144692592E0, .00002286212103119451E0, &
    .00000025347910790261E0, .00000000190451637722E0, .00000000001034969525E0, &
    .00000000000004259816E0, .00000000000000013744E0, .00000000000000000035E0 ]
  REAL, PARAMETER :: ak0cs(17) = [ -.07643947903327941E0, -.02235652605699819E0, &
    .00077341811546938E0, -.00004281006688886E0, .00000308170017386E0, &
    -.00000026393672220E0, .00000002563713036E0, -.00000000274270554E0, &
    .00000000031694296E0, -.00000000003902353E0, .00000000000506804E0, &
    -.00000000000068895E0, .00000000000009744E0, -.00000000000001427E0, &
    .00000000000000215E0, -.00000000000000033E0, .00000000000000005E0 ]
  REAL, PARAMETER :: ak02cs(14) = [ -.01201869826307592E0,-.00917485269102569E0, &
    .00014445509317750E0, -.00000401361417543E0, .00000015678318108E0, &
    -.00000000777011043E0, .00000000046111825E0, -.00000000003158592E0, &
    .00000000000243501E0, -.00000000000020743E0, .00000000000001925E0, &
    -.00000000000000192E0, .00000000000000020E0, -.00000000000000002E0 ]
  LOGICAL :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  BESK0E
  IF ( first ) THEN
    ntk0 = INITS(bk0cs,11,0.1*R1MACH(3))
    ntak0 = INITS(ak0cs,17,0.1*R1MACH(3))
    ntak02 = INITS(ak02cs,14,0.1*R1MACH(3))
    xsml = SQRT(4.0*R1MACH(3))
    first = .FALSE.
  END IF
  !
  IF ( X<=0. ) CALL XERMSG('SLATEC','BESK0E','X IS ZERO OR NEGATIVE',2,2)
  IF ( X>2. ) THEN
    !
    IF ( X<=8. ) THEN
      BESK0E = (1.25+CSEVL((16./X-5.)/3.,ak0cs,ntak0))/SQRT(X)
    ELSE
      BESK0E = (1.25+CSEVL(16./X-1.,ak02cs,ntak02))/SQRT(X)
    END IF
    RETURN
  END IF
  !
  y = 0.
  IF ( X>xsml ) y = X*X
  BESK0E = EXP(X)*(-LOG(0.5*X)*BESI0(X)-.25+CSEVL(.5*y-1.,bk0cs,ntk0))
  RETURN
END FUNCTION BESK0E
