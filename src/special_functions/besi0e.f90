!** BESI0E
REAL(SP) ELEMENTAL FUNCTION BESI0E(X)
  !> Compute the exponentially scaled modified (hyperbolic) Bessel function of
  !  the first kind of order zero.
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
  USE service, ONLY : R1MACH
  REAL(SP), INTENT(IN) :: X
  REAL(SP) :: y
  INTEGER, PARAMETER :: nti0 = 7, ntai0 = 9, ntai02 = 6
  REAL(SP), PARAMETER :: xsml = SQRT(4.5_SP*R1MACH(3))
  REAL(SP), PARAMETER :: bi0cs(12) = [ -.07660547252839144951_SP, 1.927337953993808270_SP, &
    .2282644586920301339_SP, .01304891466707290428_SP, .00043442709008164874_SP, &
    .00000942265768600193_SP, .00000014340062895106_SP, .00000000161384906966_SP, &
    .00000000001396650044_SP, .00000000000009579451_SP, .00000000000000053339_SP, &
    .00000000000000000245_SP ]
  REAL(SP), PARAMETER :: ai0cs(21) = [ .07575994494023796_SP, .00759138081082334_SP, &
    .00041531313389237_SP,  .00001070076463439_SP,-.00000790117997921_SP, &
    -.00000078261435014_SP, .00000027838499429_SP, .00000000825247260_SP, &
    -.00000001204463945_SP, .00000000155964859_SP, .00000000022925563_SP, &
    -.00000000011916228_SP, .00000000001757854_SP, .00000000000112822_SP, &
    -.00000000000114684_SP, .00000000000027155_SP,-.00000000000002415_SP, &
    -.00000000000000608_SP, .00000000000000314_SP,-.00000000000000071_SP, &
    .00000000000000007_SP ]
  REAL(SP), PARAMETER :: ai02cs(22) = [ .05449041101410882_SP, .00336911647825569_SP, &
    .00006889758346918_SP,  .00000289137052082_SP, .00000020489185893_SP, &
    .00000002266668991_SP,  .00000000339623203_SP, .00000000049406022_SP, &
    .00000000001188914_SP, -.00000000003149915_SP,-.00000000001321580_SP, &
    -.00000000000179419_SP, .00000000000071801_SP, .00000000000038529_SP, &
    .00000000000001539_SP, -.00000000000004151_SP,-.00000000000000954_SP, &
    .00000000000000382_SP, .00000000000000176_SP, -.00000000000000034_SP, &
    -.00000000000000027_SP, .00000000000000003_SP ]
  !* FIRST EXECUTABLE STATEMENT  BESI0E
  ! nti0 = INITS(bi0cs,0.1_SP*R1MACH(3))
  ! ntai0 = INITS(ai0cs,0.1_SP*R1MACH(3))
  ! ntai02 = INITS(ai02cs,0.1_SP*R1MACH(3))
  !
  y = ABS(X)
  IF( y<=xsml ) THEN
    BESI0E = 1._SP - X
  ELSEIF( y<=3._SP ) THEN
    BESI0E = EXP(-y)*(2.75_SP+CSEVL(y*y/4.5_SP-1._SP,bi0cs(1:nti0)))
  ELSEIF( y<=8. ) THEN
    BESI0E = (0.375_SP+CSEVL((48._SP/y-11._SP)/5._SP,ai0cs(1:ntai0)))/SQRT(y)
  ELSE
    BESI0E = (0.375_SP+CSEVL(16._SP/y-1._SP,ai02cs(1:ntai02)))/SQRT(y)
  END IF

  RETURN
END FUNCTION BESI0E