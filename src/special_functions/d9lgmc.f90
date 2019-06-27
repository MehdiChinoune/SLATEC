!** D9LGMC
REAL(DP) ELEMENTAL FUNCTION D9LGMC(X)
  !> Compute the log Gamma correction factor so that
  !   LOG(GAMMA(X)) = LOG(SQRT(2*PI)) + (X-5.)*LOG(X) - X + D9LGMC(X).
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
  ! Compute the log gamma correction factor for X >= 10. so that
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
  USE service, ONLY : D1MACH
  REAL(DP), INTENT(IN) :: X
  INTEGER, PARAMETER :: nalgm = 6
  REAL(DP), PARAMETER :: xbig = 1._DP/SQRT(D1MACH(3)), &
    xmax = EXP(MIN(LOG(D1MACH(2)/12._DP),-LOG(12._DP*D1MACH(1))))
  REAL(DP), PARAMETER :: algmcs(15) = [ +.1666389480451863247205729650822E+0_DP, &
    -.1384948176067563840732986059135E-4_DP, +.9810825646924729426157171547487E-8_DP, &
    -.1809129475572494194263306266719E-10_DP, +.6221098041892605227126015543416E-13_DP, &
    -.3399615005417721944303330599666E-15_DP, +.2683181998482698748957538846666E-17_DP, &
    -.2868042435334643284144622399999E-19_DP, +.3962837061046434803679306666666E-21_DP, &
    -.6831888753985766870111999999999E-23_DP, +.1429227355942498147573333333333E-24_DP, &
    -.3547598158101070547199999999999E-26_DP, +.1025680058010470912000000000000E-27_DP, &
    -.3401102254316748799999999999999E-29_DP, +.1276642195630062933333333333333E-30_DP ]
  !* FIRST EXECUTABLE STATEMENT  D9LGMC
  ! nalgm = INITDS(algmcs,D1MACH(3))
  !
  IF( X<10._DP ) THEN
    ERROR STOP 'D9LGMC : X MUST BE >= 10'
  ELSEIF( X<xbig ) THEN
    D9LGMC = DCSEVL(2._DP*(10._DP/X)**2-1._DP,algmcs(1:nalgm))/X
  ELSEIF( X<xmax ) THEN
    D9LGMC = 1._DP/(12._DP*X)
  ELSEIF( X>=xmax ) THEN
    D9LGMC = 0._DP
    ! CALL XERMSG('D9LGMC : X SO BIG D9LGMC UNDERFLOWS',2,1)
  END IF

  RETURN
END FUNCTION D9LGMC