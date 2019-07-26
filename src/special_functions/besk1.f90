!** BESK1
REAL(SP) ELEMENTAL FUNCTION BESK1(X)
  !> Compute the modified (hyperbolic) Bessel function of the third kind of order one.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C10B1
  !***
  ! **Type:**      SINGLE PRECISION (BESK1-S, DBESK1-D)
  !***
  ! **Keywords:**  FNLIB, HYPERBOLIC BESSEL FUNCTION,
  !             MODIFIED BESSEL FUNCTION, ORDER ONE, SPECIAL FUNCTIONS, THIRD KIND
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! BESK1(X) computes the modified (hyperbolic) Bessel function of third
  ! kind of order one for real argument X, where X > 0.
  !
  ! Series for BK1        on the interval  0.          to  4.00000D+00
  !                                        with weighted error   7.02E-18
  !                                         log weighted error  17.15
  !                               significant figures required  16.73
  !                                    decimal places required  17.67
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  BESI1, BESK1E, CSEVL, INITS, R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTIONsection.  (WRB)
  USE service, ONLY : tiny_sp, huge_sp, eps_2_sp
  !
  REAL(SP), INTENT(IN) :: X
  !
  REAL(SP) :: y
  INTEGER, PARAMETER :: ntk1 = 7
  REAL(SP), PARAMETER :: xmin = EXP(MAX(LOG(tiny_sp),-LOG(huge_sp))+.01_SP), &
    xsml = SQRT(4._SP*eps_2_sp), xmaxt = -LOG(tiny_sp), &
    xmax = xmaxt - 0.5_SP*xmaxt*LOG(xmaxt)/(xmaxt+0.5_SP)
  REAL(SP), PARAMETER :: bk1cs(11) = [ .0253002273389477705_SP,-.353155960776544876_SP, &
    -.122611180822657148_SP, -.0069757238596398643_SP,-.0001730288957513052_SP, &
    -.0000024334061415659_SP, -.0000000221338763073_SP,-.0000000001411488392_SP, &
    -.0000000000006666901_SP, -.0000000000000024274_SP,-.0000000000000000070_SP ]
  !* FIRST EXECUTABLE STATEMENT  BESK1
  ! ntk1 = INITS(bk1cs,0.1_SP*eps_2_sp)
  !
  IF( X<=0. ) THEN
    ERROR STOP 'BESK1 : X <= 0'
  ELSEIF( X<xmin ) THEN
    ERROR STOP 'BESK1 : X SO SMALL K1 OVERFLOWS'
  ELSEIF( X<=2._SP ) THEN
    y = 0.
    IF( X>xsml ) y = X*X
    BESK1 = LOG(0.5_SP*X)*BESI1(X) + (0.75_SP+CSEVL(.5_SP*y-1._SP,bk1cs(1:ntk1)))/X
  ELSEIF( X<=xmax ) THEN
    BESK1 = EXP(-X)*BESK1E(X)
  ELSE
    BESK1 = 0.
    ! CALL XERMSG('BESK1 : X SO BIG K1 UNDERFLOWS',1,1)
  END IF

  RETURN
END FUNCTION BESK1