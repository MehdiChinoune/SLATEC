!** D9LGMC
REAL(8) FUNCTION D9LGMC(X)
  !>
  !  Compute the log Gamma correction factor so that
  !            LOG(GAMMA(X)) = LOG(SQRT(2*PI)) + (X-5.)*LOG(X) - X
  !            + D9LGMC(X).
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C7E
  !***
  ! **Type:**      DOUBLE PRECISION (R9LGMC-S, D9LGMC-D, C9LGMC-C)
  !***
  ! **Keywords:**  COMPLETE GAMMA FUNCTION, CORRECTION TERM, FNLIB,
  !             LOG GAMMA, LOGARITHM, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! Compute the log gamma correction factor for X .GE. 10. so that
  ! LOG (GAMMA(X)) = LOG(SQRT(2*PI)) + (X-.5)*LOG(X) - X + D9lGMC(X)
  !
  ! Series for ALGM       on the interval  0.          to  1.00000E-02
  !                                        with weighted error   1.28E-31
  !                                         log weighted error  30.89
  !                               significant figures required  29.81
  !                                    decimal places required  31.48
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  D1MACH, DCSEVL, INITDS, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900720  Routine changed from user-callable to subsidiary.  (WRB)
  USE service, ONLY : XERMSG, D1MACH
  REAL(8) :: X
  INTEGER, SAVE :: nalgm
  REAL(8), PARAMETER :: xbig = 1.0D0/SQRT(D1MACH(3)), &
    xmax = EXP(MIN(LOG(D1MACH(2)/12.D0),-LOG(12.D0*D1MACH(1))))
  REAL(8), PARAMETER :: algmcs(15) = [ +.1666389480451863247205729650822D+0, &
    -.1384948176067563840732986059135D-4, +.9810825646924729426157171547487D-8, &
    -.1809129475572494194263306266719D-10, +.6221098041892605227126015543416D-13, &
    -.3399615005417721944303330599666D-15, +.2683181998482698748957538846666D-17, &
    -.2868042435334643284144622399999D-19, +.3962837061046434803679306666666D-21, &
    -.6831888753985766870111999999999D-23, +.1429227355942498147573333333333D-24, &
    -.3547598158101070547199999999999D-26, +.1025680058010470912000000000000D-27, &
    -.3401102254316748799999999999999D-29, +.1276642195630062933333333333333D-30 ]
  LOGICAL, SAVE :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  D9LGMC
  IF ( first ) THEN
    nalgm = INITDS(algmcs,15,D1MACH(3))
    first = .FALSE.
  END IF
  !
  IF ( X<10.D0 ) CALL XERMSG('D9LGMC','X MUST BE GE 10',1,2)
  IF ( X>=xmax ) THEN
    !
    D9LGMC = 0.D0
    CALL XERMSG('D9LGMC','X SO BIG D9LGMC UNDERFLOWS',2,1)
    RETURN
  END IF
  !
  D9LGMC = 1.D0/(12.D0*X)
  IF ( X<xbig ) D9LGMC = DCSEVL(2.0D0*(10.D0/X)**2-1.D0,algmcs,nalgm)/X
  RETURN
END FUNCTION D9LGMC
