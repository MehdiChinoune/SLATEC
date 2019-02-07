!*==R9LGMC.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK R9LGMC
FUNCTION R9LGMC(X)
  IMPLICIT NONE
  !*--R9LGMC5
  !*** Start of declarations inserted by SPAG
  REAL algmcs, CSEVL, R1MACH, R9LGMC, X, xbig, xmax
  INTEGER INITS, nalgm
  !*** End of declarations inserted by SPAG
  !***BEGIN PROLOGUE  R9LGMC
  !***SUBSIDIARY
  !***PURPOSE  Compute the log Gamma correction factor so that
  !            LOG(GAMMA(X)) = LOG(SQRT(2*PI)) + (X-.5)*LOG(X) - X
  !            + R9LGMC(X).
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C7E
  !***TYPE      SINGLE PRECISION (R9LGMC-S, D9LGMC-D, C9LGMC-C)
  !***KEYWORDS  COMPLETE GAMMA FUNCTION, CORRECTION TERM, FNLIB,
  !             LOG GAMMA, LOGARITHM, SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! Compute the log gamma correction factor for X .GE. 10.0 so that
  !  LOG (GAMMA(X)) = LOG(SQRT(2*PI)) + (X-.5)*LOG(X) - X + R9LGMC(X)
  !
  ! Series for ALGM       on the interval  0.          to  1.00000D-02
  !                                        with weighted error   3.40E-16
  !                                         log weighted error  15.47
  !                               significant figures required  14.39
  !                                    decimal places required  15.86
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  CSEVL, INITS, R1MACH, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770801  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900720  Routine changed from user-callable to subsidiary.  (WRB)
  !***END PROLOGUE  R9LGMC
  DIMENSION algmcs(6)
  LOGICAL first
  SAVE algmcs, nalgm, xbig, xmax, first
  DATA algmcs(1)/.166638948045186E0/
  DATA algmcs(2)/ - .0000138494817606E0/
  DATA algmcs(3)/.0000000098108256E0/
  DATA algmcs(4)/ - .0000000000180912E0/
  DATA algmcs(5)/.0000000000000622E0/
  DATA algmcs(6)/ - .0000000000000003E0/
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  R9LGMC
  IF ( first ) THEN
    nalgm = INITS(algmcs,6,R1MACH(3))
    xbig = 1.0/SQRT(R1MACH(3))
    xmax = EXP(MIN(LOG(R1MACH(2)/12.0),-LOG(12.0*R1MACH(1))))
  ENDIF
  first = .FALSE.
  !
  IF ( X<10.0 ) CALL XERMSG('SLATEC','R9LGMC','X MUST BE GE 10',1,2)
  IF ( X>=xmax ) THEN
    !
    R9LGMC = 0.0
    CALL XERMSG('SLATEC','R9LGMC','X SO BIG R9LGMC UNDERFLOWS',2,1)
    GOTO 99999
  ENDIF
  !
  R9LGMC = 1.0/(12.0*X)
  IF ( X<xbig ) R9LGMC = CSEVL(2.0*(10./X)**2-1.,algmcs,nalgm)/X
  RETURN
  !
  99999 CONTINUE
  END FUNCTION R9LGMC
