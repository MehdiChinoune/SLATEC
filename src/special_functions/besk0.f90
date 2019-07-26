!** BESK0
REAL(SP) ELEMENTAL FUNCTION BESK0(X)
  !> Compute the modified (hyperbolic) Bessel function of the third kind of order zero.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C10B1
  !***
  ! **Type:**      SINGLE PRECISION (BESK0-S, DBESK0-D)
  !***
  ! **Keywords:**  FNLIB, HYPERBOLIC BESSEL FUNCTION,
  !             MODIFIED BESSEL FUNCTION, ORDER ZERO, SPECIAL FUNCTIONS,
  !             THIRD KIND
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! BESK0(X) calculates the modified (hyperbolic) Bessel function
  ! of the third kind of order zero for real argument X > 0.0.
  !
  ! Series for BK0        on the interval  0.          to  4.00000D+00
  !                                        with weighted error   3.57E-19
  !                                         log weighted error  18.45
  !                               significant figures required  17.99
  !                                    decimal places required  18.97
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  BESI0, BESK0E, CSEVL, INITS, R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.  (WRB)
  USE service, ONLY : eps_2_sp, tiny_sp
  !
  REAL(SP), INTENT(IN) :: X
  !
  REAL(SP) :: y
  INTEGER, PARAMETER :: ntk0 = 6
  REAL(SP), PARAMETER :: xsml = SQRT(4._SP*eps_2_sp), xmaxt = -LOG(tiny_sp), &
    xmax = xmaxt - 0.5_SP*xmaxt*LOG(xmaxt)/(xmaxt+0.5_SP) - 0.01_SP
  REAL(SP), PARAMETER :: bk0cs(11) = [ -.03532739323390276872_SP, .3442898999246284869_SP, &
    .03597993651536150163_SP, .00126461541144692592_SP, .00002286212103119451_SP, &
    .00000025347910790261_SP, .00000000190451637722_SP, .00000000001034969525_SP, &
    .00000000000004259816_SP, .00000000000000013744_SP, .00000000000000000035_SP ]
  !* FIRST EXECUTABLE STATEMENT  BESK0
  ! ntk0 = INITS(bk0cs,0.1_SP*eps_2_sp)
  !
  IF( X<=0. ) THEN
    ERROR STOP 'BESK0 : X <= 0'
  ELSEIF( X<=2._SP ) THEN
    y = 0.
    IF( X>xsml ) y = X*X
    BESK0 = -LOG(0.5_SP*X)*BESI0(X) - .25_SP + CSEVL(.5_SP*y-1._SP,bk0cs(1:ntk0))
  ELSEIF( X<=xmax ) THEN
    BESK0 = EXP(-X)*BESK0E(X)
  ELSE
    BESK0 = 0.
    ! CALL XERMSG('BESK0 : X SO BIG K0 UNDERFLOWS',1,1)
  END IF

  RETURN
END FUNCTION BESK0