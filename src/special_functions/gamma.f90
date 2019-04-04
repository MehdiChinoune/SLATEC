!** GAMMA
REAL FUNCTION GAMMA(X)
  IMPLICIT NONE
  !>
  !***
  !  Compute the complete Gamma function.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C7A
  !***
  ! **Type:**      SINGLE PRECISION (GAMMA-S, DGAMMA-D, CGAMMA-C)
  !***
  ! **Keywords:**  COMPLETE GAMMA FUNCTION, FNLIB, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! GAMMA computes the gamma function at X, where X is not 0, -1, -2, ....
  ! GAMMA and X are single precision.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  CSEVL, GAMLIM, INITS, R1MACH, R9LGMC, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)

  REAL CSEVL, R1MACH, R9LGMC, sinpiy, X, y
  INTEGER i, INITS, n
  INTEGER, SAVE :: ngcs
  REAL, SAVE :: xmin, xmax, dxrel
  REAL, PARAMETER :: gcs(23) = [ .008571195590989331E0, .004415381324841007E0, &
    .05685043681599363E0,   -.004219835396418561E0,  .001326808181212460E0, &
    -.0001893024529798880E0, .0000360692532744124E0,-.0000060567619044608E0, &
    .0000010558295463022E0, -.0000001811967365542E0, .0000000311772496471E0, &
    -.0000000053542196390E0, .0000000009193275519E0,-.0000000001577941280E0, &
    .0000000000270798062E0, -.0000000000046468186E0, .0000000000007973350E0, &
    -.0000000000001368078E0, .0000000000000234731E0,-.0000000000000040274E0, &
    .0000000000000006910E0, -.0000000000000001185E0, .0000000000000000203E0 ]
  REAL, PARAMETER :: pi = 3.14159265358979324E0
  ! SQ2PIL IS LOG (SQRT (2.*PI) )
  REAL, PARAMETER :: sq2pil = 0.91893853320467274E0
  LOGICAL :: first = .TRUE.
  !
  ! LANL DEPENDENT CODE REMOVED 81.02.04
  !
  !* FIRST EXECUTABLE STATEMENT  GAMMA
  IF ( first ) THEN
    !
    ! ---------------------------------------------------------------------
    ! INITIALIZE.  FIND LEGAL BOUNDS FOR X, AND DETERMINE THE NUMBER OF
    ! TERMS IN THE SERIES REQUIRED TO ATTAIN AN ACCURACY TEN TIMES BETTER
    ! THAN MACHINE PRECISION.
    !
    ngcs = INITS(gcs,23,0.1*R1MACH(3))
    !
    CALL GAMLIM(xmin,xmax)
    dxrel = SQRT(R1MACH(4))
    !
    ! ---------------------------------------------------------------------
    ! FINISH INITIALIZATION.  START EVALUATING GAMMA(X).
    !
    first = .FALSE.
  END IF
  !
  y = ABS(X)
  IF ( y>10.0 ) THEN
    !
    ! COMPUTE GAMMA(X) FOR ABS(X) .GT. 10.0.  RECALL Y = ABS(X).
    !
    IF ( X>xmax ) CALL XERMSG('SLATEC','GAMMA','X SO BIG GAMMA OVERFLOWS',3,2)
    !
    GAMMA = 0.
    IF ( X<xmin ) CALL XERMSG('SLATEC','GAMMA','X SO SMALL GAMMA UNDERFLOWS'&
      ,2,1)
    IF ( X<xmin ) RETURN
    !
    GAMMA = EXP((y-0.5)*LOG(y)-y+sq2pil+R9LGMC(y))
    IF ( X>0. ) RETURN
    !
    IF ( ABS((X-AINT(X-0.5))/X)<dxrel ) CALL XERMSG('SLATEC','GAMMA',&
      'ANSWER LT HALF PRECISION, X TOO NEAR NEGATIVE INTEGER',1,1)
    !
    sinpiy = SIN(pi*y)
    IF ( sinpiy==0. ) CALL XERMSG('SLATEC','GAMMA','X IS A NEGATIVE INTEGER'&
      ,4,2)
    !
    GAMMA = -pi/(y*sinpiy*GAMMA)
    RETURN
  ELSE
    !
    ! COMPUTE GAMMA(X) FOR ABS(X) .LE. 10.0.  REDUCE INTERVAL AND
    ! FIND GAMMA(1+Y) FOR 0. .LE. Y .LT. 1. FIRST OF ALL.
    !
    n = INT( X )
    IF ( X<0. ) n = n - 1
    y = X - n
    n = n - 1
    GAMMA = 0.9375 + CSEVL(2.*y-1.,gcs,ngcs)
    IF ( n==0 ) RETURN
    !
    IF ( n<=0 ) THEN
      !
      ! COMPUTE GAMMA(X) FOR X .LT. 1.
      !
      n = -n
      IF ( X==0. ) CALL XERMSG('SLATEC','GAMMA','X IS 0',4,2)
      IF ( X<0..AND.X+n-2==0. )&
        CALL XERMSG('SLATEC','GAMMA','X IS A NEGATIVE INTEGER',4,2)
      IF ( X<(-0.5).AND.ABS((X-AINT(X-0.5))/X)<dxrel )&
        CALL XERMSG('SLATEC','GAMMA',&
        'ANSWER LT HALF PRECISION BECAUSE X TOO NEAR NEGATIVE INTEGER',1,1)
      !
      DO i = 1, n
        GAMMA = GAMMA/(X+i-1)
      END DO
      RETURN
    END IF
  END IF
  !
  ! GAMMA(X) FOR X .GE. 2.
  !
  DO i = 1, n
    GAMMA = (y+i)*GAMMA
  END DO
  RETURN
END FUNCTION GAMMA
