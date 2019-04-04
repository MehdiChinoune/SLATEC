!** BESI0E
REAL FUNCTION BESI0E(X)
  IMPLICIT NONE
  !>
  !***
  !  Compute the exponentially scaled modified (hyperbolic)
  !            Bessel function of the first kind of order zero.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C10B1
  !***
  ! **Type:**      SINGLE PRECISION (BESI0E-S, DBSI0E-D)
  !***
  ! **Keywords:**  EXPONENTIALLY SCALED, FIRST KIND, FNLIB,
  !             HYPERBOLIC BESSEL FUNCTION, MODIFIED BESSEL FUNCTION,
  !             ORDER ZERO, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! BESI0E(X) calculates the exponentially scaled modified (hyperbolic)
  ! Bessel function of the first kind of order zero for real argument X;
  ! i.e., EXP(-ABS(X))*I0(X).
  !
  !
  ! Series for BI0        on the interval  0.          to  9.00000D+00
  !                                        with weighted error   2.46E-18
  !                                         log weighted error  17.61
  !                               significant figures required  17.90
  !                                    decimal places required  18.15
  !
  !
  ! Series for AI0        on the interval  1.25000D-01 to  3.33333D-01
  !                                        with weighted error   7.87E-17
  !                                         log weighted error  16.10
  !                               significant figures required  14.69
  !                                    decimal places required  16.76
  !
  !
  ! Series for AI02       on the interval  0.          to  1.25000D-01
  !                                        with weighted error   3.79E-17
  !                                         log weighted error  16.42
  !                               significant figures required  14.86
  !                                    decimal places required  17.09
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  CSEVL, INITS, R1MACH

  !* REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890313  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)

  REAL CSEVL, R1MACH, X, y
  INTEGER INITS
  INTEGER, SAVE :: nti0, ntai0, ntai02
  REAL, SAVE :: xsml
  REAL, PARAMETER :: bi0cs(12) = [ -.07660547252839144951E0, 1.927337953993808270E0, &
    .2282644586920301339E0, .01304891466707290428E0, .00043442709008164874E0, &
    .00000942265768600193E0, .00000014340062895106E0, .00000000161384906966E0, &
    .00000000001396650044E0, .00000000000009579451E0, .00000000000000053339E0, &
    .00000000000000000245E0 ]
  REAL, PARAMETER :: ai0cs(21) = [ .07575994494023796E0, .00759138081082334E0, &
    .00041531313389237E0,  .00001070076463439E0,-.00000790117997921E0, &
    -.00000078261435014E0, .00000027838499429E0, .00000000825247260E0, &
    -.00000001204463945E0, .00000000155964859E0, .00000000022925563E0, &
    -.00000000011916228E0, .00000000001757854E0, .00000000000112822E0, &
    -.00000000000114684E0, .00000000000027155E0,-.00000000000002415E0, &
    -.00000000000000608E0, .00000000000000314E0,-.00000000000000071E0, &
    .00000000000000007E0 ]
  REAL, PARAMETER :: ai02cs(22) = [ .05449041101410882E0, .00336911647825569E0, &
    .00006889758346918E0,  .00000289137052082E0, .00000020489185893E0, &
    .00000002266668991E0,  .00000000339623203E0, .00000000049406022E0, &
    .00000000001188914E0, -.00000000003149915E0,-.00000000001321580E0, &
    -.00000000000179419E0, .00000000000071801E0, .00000000000038529E0, &
    .00000000000001539E0, -.00000000000004151E0,-.00000000000000954E0, &
    .00000000000000382E0, .00000000000000176E0, -.00000000000000034E0, &
    -.00000000000000027E0, .00000000000000003E0 ]
  LOGICAL :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  BESI0E
  IF ( first ) THEN
    nti0 = INITS(bi0cs,12,0.1*R1MACH(3))
    ntai0 = INITS(ai0cs,21,0.1*R1MACH(3))
    ntai02 = INITS(ai02cs,22,0.1*R1MACH(3))
    xsml = SQRT(4.5*R1MACH(3))
    first = .FALSE.
  END IF
  !
  y = ABS(X)
  IF ( y>3.0 ) THEN
    !
    IF ( y<=8. ) BESI0E = (.375+CSEVL((48./y-11.)/5.,ai0cs,ntai0))/SQRT(y)
    IF ( y>8. ) BESI0E = (.375+CSEVL(16./y-1.,ai02cs,ntai02))/SQRT(y)
    RETURN
  END IF
  !
  BESI0E = 1.0 - X
  IF ( y>xsml ) BESI0E = EXP(-y)*(2.75+CSEVL(y*y/4.5-1.0,bi0cs,nti0))
  RETURN
END FUNCTION BESI0E
