!** DBESI1
REAL(DP) ELEMENTAL FUNCTION DBESI1(X)
  !> Compute the modified (hyperbolic) Bessel function of the first kind of order one.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C10B1
  !***
  ! **Type:**      DOUBLE PRECISION (BESI1-S, DBESI1-D)
  !***
  ! **Keywords:**  FIRST KIND, FNLIB, HYPERBOLIC BESSEL FUNCTION,
  !             MODIFIED BESSEL FUNCTION, ORDER ONE, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! DBESI1(X) calculates the double precision modified (hyperbolic)
  ! Bessel function of the first kind of order one and double precision
  ! argument X.
  !
  ! Series for BI1        on the interval  0.          to  9.00000E+00
  !                                        with weighted error   1.44E-32
  !                                         log weighted error  31.84
  !                               significant figures required  31.45
  !                                    decimal places required  32.46
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  D1MACH, DBSI1E, DCSEVL, INITDS, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  USE service, ONLY : D1MACH
  REAL(DP), INTENT(IN) :: X
  REAL(DP) :: y
  INTEGER, PARAMETER :: nti1 = 11
  REAL(DP), PARAMETER :: xmin = 2._DP*D1MACH(1), xsml = SQRT(4.5_DP*D1MACH(3)), &
    xmax = LOG(D1MACH(2))
  REAL(DP), PARAMETER :: bi1cs(17) = [ -.19717132610998597316138503218149E-2_DP, &
    +.40734887667546480608155393652014E+0_DP, +.34838994299959455866245037783787E-1_DP, &
    +.15453945563001236038598401058489E-2_DP, +.41888521098377784129458832004120E-4_DP, &
    +.76490267648362114741959703966069E-6_DP, +.10042493924741178689179808037238E-7_DP, &
    +.99322077919238106481371298054863E-10_DP, +.76638017918447637275200171681349E-12_DP, &
    +.47414189238167394980388091948160E-14_DP, +.24041144040745181799863172032000E-16_DP, &
    +.10171505007093713649121100799999E-18_DP, +.36450935657866949458491733333333E-21_DP, &
    +.11205749502562039344810666666666E-23_DP, +.29875441934468088832000000000000E-26_DP, &
    +.69732310939194709333333333333333E-29_DP, +.14367948220620800000000000000000E-31_DP ]
  !* FIRST EXECUTABLE STATEMENT  DBESI1
  ! nti1 = INITDS(bi1cs,0.1_SP*D1MACH(3))
  !
  y = ABS(X)
  IF( y>xmax ) THEN
    ERROR STOP 'DBESI1 : ABS(X) SO BIG I1 OVERFLOWS'
  ELSEIF( y>3._DP ) THEN
    DBESI1 = EXP(y)*DBSI1E(X)
  ELSEIF( y>xsml ) THEN
    DBESI1 = X*(0.875_DP+DCSEVL(y*y/4.5_DP-1._DP,bi1cs(1:nti1)))
  ELSEIF( y>xmin ) THEN
    DBESI1 = 0.5_DP*X
  ELSE
    ! IF( y<=xmin ) CALL XERMSG('DBESI1','ABS(X) SO SMALL I1 UNDERFLOWS',1,1)
    DBESI1 = 0._DP
  END IF

  RETURN
END FUNCTION DBESI1