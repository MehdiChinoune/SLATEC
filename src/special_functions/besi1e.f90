!** BESI1E
REAL(SP) ELEMENTAL FUNCTION BESI1E(X)
  !> Compute the exponentially scaled modified (hyperbolic) Bessel function of
  !  the first kind of order one.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C10B1
  !***
  ! **Type:**      SINGLE PRECISION (BESI1E-S, DBSI1E-D)
  !***
  ! **Keywords:**  EXPONENTIALLY SCALED, FIRST KIND, FNLIB,
  !             HYPERBOLIC BESSEL FUNCTION, MODIFIED BESSEL FUNCTION,
  !             ORDER ONE, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! BESI1E(X) calculates the exponentially scaled modified (hyperbolic)
  ! Bessel function of the first kind of order one for real argument X;
  ! i.e., EXP(-ABS(X))*I1(X).
  !
  ! Series for BI1        on the interval  0.          to  9.00000D+00
  !                                        with weighted error   2.40E-17
  !                                         log weighted error  16.62
  !                               significant figures required  16.23
  !                                    decimal places required  17.14
  !
  ! Series for AI1        on the interval  1.25000D-01 to  3.33333D-01
  !                                        with weighted error   6.98E-17
  !                                         log weighted error  16.16
  !                               significant figures required  14.53
  !                                    decimal places required  16.82
  !
  ! Series for AI12       on the interval  0.          to  1.25000D-01
  !                                        with weighted error   3.55E-17
  !                                         log weighted error  16.45
  !                               significant figures required  14.69
  !                                    decimal places required  17.12
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  CSEVL, INITS, R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   890210  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.  (WRB)
  !   920618  Removed space from variable names.  (RWC, WRB)
  USE service, ONLY : R1MACH
  REAL(SP), INTENT(IN) :: X
  REAL(SP) :: y
  INTEGER, PARAMETER :: nti1 = 7, ntai1 = 9, ntai12 = 6
  REAL(SP), PARAMETER :: xmin = 2._SP*R1MACH(1), xsml = SQRT(4.5_SP*R1MACH(3))
  REAL(SP), PARAMETER :: bi1cs(11) = [ -.001971713261099859_SP, .40734887667546481_SP, &
    .034838994299959456_SP, .001545394556300123_SP, .000041888521098377_SP, &
    .000000764902676483_SP, .000000010042493924_SP, .000000000099322077_SP, &
    .000000000000766380_SP, .000000000000004741_SP, .000000000000000024_SP ]
  REAL(SP), PARAMETER :: ai1cs(21) = [ -.02846744181881479_SP,-.01922953231443221_SP, &
    -.00061151858579437_SP,-.00002069971253350_SP, .00000858561914581_SP, &
    .00000104949824671_SP, -.00000029183389184_SP,-.00000001559378146_SP, &
    .00000001318012367_SP, -.00000000144842341_SP,-.00000000029085122_SP, &
    .00000000012663889_SP, -.00000000001664947_SP,-.00000000000166665_SP, &
    .00000000000124260_SP, -.00000000000027315_SP, .00000000000002023_SP, &
    .00000000000000730_SP, -.00000000000000333_SP, .00000000000000071_SP, &
    -.00000000000000006_SP ]
  REAL(SP), PARAMETER :: ai12cs(22) = [ .02857623501828014_SP,-.00976109749136147_SP, &
    -.00011058893876263_SP,-.00000388256480887_SP,-.00000025122362377_SP, &
    -.00000002631468847_SP,-.00000000383538039_SP,-.00000000055897433_SP, &
    -.00000000001897495_SP, .00000000003252602_SP, .00000000001412580_SP, &
    .00000000000203564_SP, -.00000000000071985_SP,-.00000000000040836_SP, &
    -.00000000000002101_SP, .00000000000004273_SP, .00000000000001041_SP, &
    -.00000000000000382_SP,-.00000000000000186_SP, .00000000000000033_SP, &
    .00000000000000028_SP, -.00000000000000003_SP ]
  !* FIRST EXECUTABLE STATEMENT  BESI1E
  ! nti1 = INITS(bi1cs,0.1_SP*R1MACH(3))
  ! ntai1 = INITS(ai1cs,0.1_SP*R1MACH(3))
  ! ntai12 = INITS(ai12cs,0.1_SP*R1MACH(3))
  !
  y = ABS(X)
  IF( y>3._SP ) THEN
    IF( y<=8. ) THEN
      BESI1E = (0.375_SP+CSEVL((48._SP/y-11._SP)/5._SP,ai1cs(1:ntai1)))/SQRT(y)
    ELSE
      BESI1E = (0.375_SP+CSEVL(16._SP/y-1._SP,ai12cs(1:ntai12)))/SQRT(y)
    END IF
    BESI1E = SIGN(BESI1E,X)
  ELSEIF( y>xmin ) THEN
    IF( y<=xsml ) THEN
      BESI1E = 0.5_SP*X
    ELSE
      BESI1E = X*(0.875_SP+CSEVL(y*y/4.5_SP-1._SP,bi1cs(1:nti1)))
    END IF
    BESI1E = EXP(-y)*BESI1E
  ELSE
    !CALL XERMSG('BESI1E : ABS(X) SO SMALL I1 UNDERFLOWS',1,1)
    BESI1E = 0._SP
  END IF

  RETURN
END FUNCTION BESI1E