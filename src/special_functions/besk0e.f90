!** BESK0E
REAL(SP) ELEMENTAL FUNCTION BESK0E(X)
  !> Compute the exponentially scaled modified (hyperbolic)  Bessel function
  !  of the third kind of order zero.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C10B1
  !***
  ! **Type:**      SINGLE PRECISION (BESK0E-S, DBSK0E-D)
  !***
  ! **Keywords:**  EXPONENTIALLY SCALED, FNLIB, HYPERBOLIC BESSEL FUNCTION,
  !             MODIFIED BESSEL FUNCTION, ORDER ZERO, SPECIAL FUNCTIONS, THIRD KIND
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! BESK0E(X) computes the exponentially scaled modified (hyperbolic)
  ! Bessel function of third kind of order zero for real argument
  ! X > 0.0, i.e., EXP(X)*K0(X).
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
  !   900326  Removed duplicate information from DESCRIPTION section.  (WRB)
  USE service, ONLY : R1MACH
  REAL(SP), INTENT(IN) :: X
  REAL(SP) :: y
  INTEGER, PARAMETER :: ntk0 = 6, ntak0 = 7, ntak02 = 6
  REAL(SP), PARAMETER :: xsml = SQRT(4._SP*R1MACH(3))
  REAL(SP), PARAMETER :: bk0cs(11) = [ -.03532739323390276872_SP, .3442898999246284869_SP, &
    .03597993651536150163_SP, .00126461541144692592_SP, .00002286212103119451_SP, &
    .00000025347910790261_SP, .00000000190451637722_SP, .00000000001034969525_SP, &
    .00000000000004259816_SP, .00000000000000013744_SP, .00000000000000000035_SP ]
  REAL(SP), PARAMETER :: ak0cs(17) = [ -.07643947903327941_SP, -.02235652605699819_SP, &
    .00077341811546938_SP, -.00004281006688886_SP, .00000308170017386_SP, &
    -.00000026393672220_SP, .00000002563713036_SP, -.00000000274270554_SP, &
    .00000000031694296_SP, -.00000000003902353_SP, .00000000000506804_SP, &
    -.00000000000068895_SP, .00000000000009744_SP, -.00000000000001427_SP, &
    .00000000000000215_SP, -.00000000000000033_SP, .00000000000000005_SP ]
  REAL(SP), PARAMETER :: ak02cs(14) = [ -.01201869826307592_SP,-.00917485269102569_SP, &
    .00014445509317750_SP, -.00000401361417543_SP, .00000015678318108_SP, &
    -.00000000777011043_SP, .00000000046111825_SP, -.00000000003158592_SP, &
    .00000000000243501_SP, -.00000000000020743_SP, .00000000000001925_SP, &
    -.00000000000000192_SP, .00000000000000020_SP, -.00000000000000002_SP ]
  !* FIRST EXECUTABLE STATEMENT  BESK0E
  ! ntk0 = INITS(bk0cs,0.1_SP*R1MACH(3))
  ! ntak0 = INITS(ak0cs,0.1_SP*R1MACH(3))
  ! ntak02 = INITS(ak02cs,0.1_SP*R1MACH(3))
  !
  IF( X<=0. ) THEN
    ERROR STOP 'BESK0E : X <= 0'
  ELSEIF( X<=2._SP ) THEN
    y = 0.
    IF( X>xsml ) y = X*X
    BESK0E = EXP(X)*(-LOG(0.5_SP*X)*BESI0(X)-.25_SP+CSEVL(.5_SP*y-1._SP,bk0cs(1:ntk0)))
  ELSEIF( X<=8._SP ) THEN
    BESK0E = (1.25_SP+CSEVL((16._SP/X-5._SP)/3._SP,ak0cs(1:ntak0)))/SQRT(X)
  ELSE
    BESK0E = (1.25_SP+CSEVL(16._SP/X-1._SP,ak02cs(1:ntak02)))/SQRT(X)
  END IF

  RETURN
END FUNCTION BESK0E