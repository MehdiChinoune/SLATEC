!** BESK1E
REAL FUNCTION BESK1E(X)
  IMPLICIT NONE
  !>
  !***
  !  Compute the exponentially scaled modified (hyperbolic)
  !            Bessel function of the third kind of order one.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C10B1
  !***
  ! **Type:**      SINGLE PRECISION (BESK1E-S, DBSK1E-D)
  !***
  ! **Keywords:**  EXPONENTIALLY SCALED, FNLIB, HYPERBOLIC BESSEL FUNCTION,
  !             MODIFIED BESSEL FUNCTION, ORDER ONE, SPECIAL FUNCTIONS,
  !             THIRD KIND
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! BESK1E(X) computes the exponentially scaled modified (hyperbolic)
  ! Bessel function of third kind of order one for real argument
  ! X .GT. 0.0, i.e., EXP(X)*K1(X).
  !
  ! Series for BK1        on the interval  0.          to  4.00000D+00
  !                                        with weighted error   7.02E-18
  !                                         log weighted error  17.15
  !                               significant figures required  16.73
  !                                    decimal places required  17.67
  !
  ! Series for AK1        on the interval  1.25000D-01 to  5.00000D-01
  !                                        with weighted error   6.06E-17
  !                                         log weighted error  16.22
  !                               significant figures required  15.41
  !                                    decimal places required  16.83
  !
  ! Series for AK12       on the interval  0.          to  1.25000D-01
  !                                        with weighted error   2.58E-17
  !                                         log weighted error  16.59
  !                               significant figures required  15.22
  !                                    decimal places required  17.16
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  BESI1, CSEVL, INITS, R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)

  REAL BESI1, CSEVL, R1MACH, X, y
  INTEGER INITS
  INTEGER, SAVE :: ntk1, ntak1, ntak12
  REAL, SAVE :: xmin, xsml
  REAL, PARAMETER :: bk1cs(11) = [ .0253002273389477705E0,-.353155960776544876E0, &
    -.122611180822657148E0, -.0069757238596398643E0,-.0001730288957513052E0, &
    -.0000024334061415659E0, -.0000000221338763073E0,-.0000000001411488392E0, &
    -.0000000000006666901E0, -.0000000000000024274E0,-.0000000000000000070E0 ]
  REAL, PARAMETER :: ak1cs(17) = [ .2744313406973883E0, .07571989953199368E0, &
    -.00144105155647540E0, .00006650116955125E0, -.00000436998470952E0, &
    .00000035402774997E0, -.00000003311163779E0, .00000000344597758E0, &
    -.00000000038989323E0, .00000000004720819E0, -.00000000000604783E0, &
    .00000000000081284E0, -.00000000000011386E0, .00000000000001654E0, &
    -.00000000000000248E0, .00000000000000038E0, -.00000000000000006E0 ]
  REAL, PARAMETER :: ak12cs(14) = [ .06379308343739001E0, .02832887813049721E0, &
    -.00024753706739052E0, .00000577197245160E0, -.00000020689392195E0, &
    .00000000973998344E0, -.00000000055853361E0, .00000000003732996E0, &
    -.00000000000282505E0, .00000000000023720E0, -.00000000000002176E0, &
    .00000000000000215E0, -.00000000000000022E0, .00000000000000002E0 ]
  LOGICAL :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  BESK1E
  IF ( first ) THEN
    ntk1 = INITS(bk1cs,11,0.1*R1MACH(3))
    ntak1 = INITS(ak1cs,17,0.1*R1MACH(3))
    ntak12 = INITS(ak12cs,14,0.1*R1MACH(3))
    !
    xmin = EXP(MAX(LOG(R1MACH(1)),-LOG(R1MACH(2)))+.01)
    xsml = SQRT(4.0*R1MACH(3))
    first = .FALSE.
  ENDIF
  !
  IF ( X<=0. ) CALL XERMSG('SLATEC','BESK1E','X IS ZERO OR NEGATIVE',2,2)
  IF ( X>2.0 ) THEN
    !
    IF ( X<=8. ) THEN
      BESK1E = (1.25+CSEVL((16./X-5.)/3.,ak1cs,ntak1))/SQRT(X)
    ELSE
      BESK1E = (1.25+CSEVL(16./X-1.,ak12cs,ntak12))/SQRT(X)
    ENDIF
    RETURN
  ENDIF
  !
  IF ( X<xmin ) CALL XERMSG('SLATEC','BESK1E','X SO SMALL K1 OVERFLOWS',3,2)
  y = 0.
  IF ( X>xsml ) y = X*X
  BESK1E = EXP(X)*(LOG(0.5*X)*BESI1(X)+(0.75+CSEVL(.5*y-1.,bk1cs,ntk1))/X)
  RETURN
END FUNCTION BESK1E
