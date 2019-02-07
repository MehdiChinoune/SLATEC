!*==D9LGMC.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK D9LGMC
REAL(8) FUNCTION D9LGMC(X)
  IMPLICIT NONE
  !*--D9LGMC5
  !*** Start of declarations inserted by SPAG
  INTEGER INITDS, nalgm
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  D9LGMC
  !***SUBSIDIARY
  !***PURPOSE  Compute the log Gamma correction factor so that
  !            LOG(DGAMMA(X)) = LOG(SQRT(2*PI)) + (X-5.)*LOG(X) - X
  !            + D9LGMC(X).
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C7E
  !***TYPE      DOUBLE PRECISION (R9LGMC-S, D9LGMC-D, C9LGMC-C)
  !***KEYWORDS  COMPLETE GAMMA FUNCTION, CORRECTION TERM, FNLIB,
  !             LOG GAMMA, LOGARITHM, SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! Compute the log gamma correction factor for X .GE. 10. so that
  ! LOG (DGAMMA(X)) = LOG(SQRT(2*PI)) + (X-.5)*LOG(X) - X + D9lGMC(X)
  !
  ! Series for ALGM       on the interval  0.          to  1.00000E-02
  !                                        with weighted error   1.28E-31
  !                                         log weighted error  30.89
  !                               significant figures required  29.81
  !                                    decimal places required  31.48
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  D1MACH, DCSEVL, INITDS, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900720  Routine changed from user-callable to subsidiary.  (WRB)
  !***END PROLOGUE  D9LGMC
  REAL(8) :: X, algmcs(15), xbig, xmax, DCSEVL, D1MACH
  LOGICAL first
  SAVE algmcs, nalgm, xbig, xmax, first
  DATA algmcs(1)/ + .1666389480451863247205729650822D+0/
  DATA algmcs(2)/ - .1384948176067563840732986059135D-4/
  DATA algmcs(3)/ + .9810825646924729426157171547487D-8/
  DATA algmcs(4)/ - .1809129475572494194263306266719D-10/
  DATA algmcs(5)/ + .6221098041892605227126015543416D-13/
  DATA algmcs(6)/ - .3399615005417721944303330599666D-15/
  DATA algmcs(7)/ + .2683181998482698748957538846666D-17/
  DATA algmcs(8)/ - .2868042435334643284144622399999D-19/
  DATA algmcs(9)/ + .3962837061046434803679306666666D-21/
  DATA algmcs(10)/ - .6831888753985766870111999999999D-23/
  DATA algmcs(11)/ + .1429227355942498147573333333333D-24/
  DATA algmcs(12)/ - .3547598158101070547199999999999D-26/
  DATA algmcs(13)/ + .1025680058010470912000000000000D-27/
  DATA algmcs(14)/ - .3401102254316748799999999999999D-29/
  DATA algmcs(15)/ + .1276642195630062933333333333333D-30/
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  D9LGMC
  IF ( first ) THEN
    nalgm = INITDS(algmcs,15,REAL(D1MACH(3)))
    xbig = 1.0D0/SQRT(D1MACH(3))
    xmax = EXP(MIN(LOG(D1MACH(2)/12.D0),-LOG(12.D0*D1MACH(1))))
  ENDIF
  first = .FALSE.
  !
  IF ( X<10.D0 ) CALL XERMSG('SLATEC','D9LGMC','X MUST BE GE 10',1,2)
  IF ( X>=xmax ) THEN
    !
    D9LGMC = 0.D0
    CALL XERMSG('SLATEC','D9LGMC','X SO BIG D9LGMC UNDERFLOWS',2,1)
    GOTO 99999
  ENDIF
  !
  D9LGMC = 1.D0/(12.D0*X)
  IF ( X<xbig ) D9LGMC = DCSEVL(2.0D0*(10.D0/X)**2-1.D0,algmcs,nalgm)/X
  RETURN
  !
  99999 CONTINUE
  END FUNCTION D9LGMC
