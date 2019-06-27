!** BESK1E
REAL(SP) ELEMENTAL FUNCTION BESK1E(X)
  !> Compute the exponentially scaled modified (hyperbolic) Bessel function
  !  of the third kind of order one.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C10B1
  !***
  ! **Type:**      SINGLE PRECISION (BESK1E-S, DBSK1E-D)
  !***
  ! **Keywords:**  EXPONENTIALLY SCALED, FNLIB, HYPERBOLIC BESSEL FUNCTION,
  !             MODIFIED BESSEL FUNCTION, ORDER ONE, SPECIAL FUNCTIONS, THIRD KIND
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! BESK1E(X) computes the exponentially scaled modified (hyperbolic)
  ! Bessel function of third kind of order one for real argument
  ! X > 0.0, i.e., EXP(X)*K1(X).
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
  !   900326  Removed duplicate information from DESCRIPTION section.  (WRB)
  USE service, ONLY : R1MACH
  REAL(SP), INTENT(IN) :: X
  REAL(SP) :: y
  INTEGER, PARAMETER :: ntk1 = 7, ntak1 = 7, ntak12 = 6
  REAL(SP), PARAMETER :: xmin = EXP(MAX(LOG(R1MACH(1)),-LOG(R1MACH(2)))+.01_SP), &
    xsml = SQRT(4._SP*R1MACH(3))
  REAL(SP), PARAMETER :: bk1cs(11) = [ .0253002273389477705_SP,-.353155960776544876_SP, &
    -.122611180822657148_SP, -.0069757238596398643_SP,-.0001730288957513052_SP, &
    -.0000024334061415659_SP, -.0000000221338763073_SP,-.0000000001411488392_SP, &
    -.0000000000006666901_SP, -.0000000000000024274_SP,-.0000000000000000070_SP ]
  REAL(SP), PARAMETER :: ak1cs(17) = [ .2744313406973883_SP, .07571989953199368_SP, &
    -.00144105155647540_SP, .00006650116955125_SP, -.00000436998470952_SP, &
    .00000035402774997_SP, -.00000003311163779_SP, .00000000344597758_SP, &
    -.00000000038989323_SP, .00000000004720819_SP, -.00000000000604783_SP, &
    .00000000000081284_SP, -.00000000000011386_SP, .00000000000001654_SP, &
    -.00000000000000248_SP, .00000000000000038_SP, -.00000000000000006_SP ]
  REAL(SP), PARAMETER :: ak12cs(14) = [ .06379308343739001_SP, .02832887813049721_SP, &
    -.00024753706739052_SP, .00000577197245160_SP, -.00000020689392195_SP, &
    .00000000973998344_SP, -.00000000055853361_SP, .00000000003732996_SP, &
    -.00000000000282505_SP, .00000000000023720_SP, -.00000000000002176_SP, &
    .00000000000000215_SP, -.00000000000000022_SP, .00000000000000002_SP ]
  !* FIRST EXECUTABLE STATEMENT  BESK1E
  ! ntk1 = INITS(bk1cs,0.1_SP*R1MACH(3))
  ! ntak1 = INITS(ak1cs,0.1_SP*R1MACH(3))
  ! ntak12 = INITS(ak12cs,0.1_SP*R1MACH(3))
  !
  IF( X<=0. ) THEN
    ERROR STOP 'BESK1E : X <= 0'
  ELSEIF( X<xmin ) THEN
    ERROR STOP 'BESK1E : X SO SMALL K1 OVERFLOWS'
  ELSEIF( X<=2._SP ) THEN
    y = 0.
    IF( X>xsml ) y = X*X
    BESK1E = EXP(X)*(LOG(0.5_SP*X)*BESI1(X)+(0.75_SP+CSEVL(.5_SP*y-1._SP,bk1cs(1:ntk1)))/X)
  ELSEIF( X<=8. ) THEN
    BESK1E = (1.25_SP+CSEVL((16._SP/X-5._SP)/3._SP,ak1cs(1:ntak1)))/SQRT(X)
  ELSE
    BESK1E = (1.25_SP+CSEVL(16._SP/X-1._SP,ak12cs(1:ntak12)))/SQRT(X)
  END IF

  RETURN
END FUNCTION BESK1E