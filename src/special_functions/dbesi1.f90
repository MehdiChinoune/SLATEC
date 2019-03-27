!** DBESI1
REAL(8) FUNCTION DBESI1(X)
  IMPLICIT NONE
  !>
  !***
  !  Compute the modified (hyperbolic) Bessel function of the
  !            first kind of order one.
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

  INTEGER INITDS, nti1
  REAL(8) :: X, xmax, xmin, xsml, y, D1MACH, DCSEVL, DBSI1E
  SAVE nti1, xmin, xsml, xmax
  REAL(8), PARAMETER :: bi1cs(17) = [ -.19717132610998597316138503218149D-2, &
    +.40734887667546480608155393652014D+0, +.34838994299959455866245037783787D-1, &
    +.15453945563001236038598401058489D-2, +.41888521098377784129458832004120D-4, &
    +.76490267648362114741959703966069D-6, +.10042493924741178689179808037238D-7, &
    +.99322077919238106481371298054863D-10, +.76638017918447637275200171681349D-12, &
    +.47414189238167394980388091948160D-14, +.24041144040745181799863172032000D-16, &
    +.10171505007093713649121100799999D-18, +.36450935657866949458491733333333D-21, &
    +.11205749502562039344810666666666D-23, +.29875441934468088832000000000000D-26, &
    +.69732310939194709333333333333333D-29, +.14367948220620800000000000000000D-31 ]
  LOGICAL :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  DBESI1
  IF ( first ) THEN
    nti1 = INITDS(bi1cs,17,0.1*REAL(D1MACH(3)))
    xmin = 2.0D0*D1MACH(1)
    xsml = SQRT(4.5D0*D1MACH(3))
    xmax = LOG(D1MACH(2))
    first = .FALSE.
  ENDIF
  !
  y = ABS(X)
  IF ( y>3.0D0 ) THEN
    !
    IF ( y>xmax ) CALL XERMSG('SLATEC','DBESI1','ABS(X) SO BIG I1 OVERFLOWS',2,2)
    !
    DBESI1 = EXP(y)*DBSI1E(X)
    RETURN
  ENDIF
  !
  DBESI1 = 0.D0
  IF ( y==0.D0 ) RETURN
  !
  IF ( y<=xmin ) CALL XERMSG('SLATEC','DBESI1',&
    'ABS(X) SO SMALL I1 UNDERFLOWS',1,1)
  IF ( y>xmin ) DBESI1 = 0.5D0*X
  IF ( y>xsml ) DBESI1 = X*(0.875D0+DCSEVL(y*y/4.5D0-1.D0,bi1cs,nti1))
  RETURN
END FUNCTION DBESI1
