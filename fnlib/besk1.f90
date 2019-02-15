!DECK BESK1
FUNCTION BESK1(X)
  IMPLICIT NONE
  REAL BESI1, BESK1, BESK1E, bk1cs, CSEVL, R1MACH, X, xmax, xmaxt, &
    xmin, xsml, y
  INTEGER INITS, ntk1
  !***BEGIN PROLOGUE  BESK1
  !***PURPOSE  Compute the modified (hyperbolic) Bessel function of the
  !            third kind of order one.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C10B1
  !***TYPE      SINGLE PRECISION (BESK1-S, DBESK1-D)
  !***KEYWORDS  FNLIB, HYPERBOLIC BESSEL FUNCTION,
  !             MODIFIED BESSEL FUNCTION, ORDER ONE, SPECIAL FUNCTIONS,
  !             THIRD KIND
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! BESK1(X) computes the modified (hyperbolic) Bessel function of third
  ! kind of order one for real argument X, where X .GT. 0.
  !
  ! Series for BK1        on the interval  0.          to  4.00000D+00
  !                                        with weighted error   7.02E-18
  !                                         log weighted error  17.15
  !                               significant figures required  16.73
  !                                    decimal places required  17.67
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  BESI1, BESK1E, CSEVL, INITS, R1MACH, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !***END PROLOGUE  BESK1
  DIMENSION bk1cs(11)
  LOGICAL first
  SAVE bk1cs, ntk1, xmin, xsml, xmax, first
  DATA bk1cs(1)/.0253002273389477705E0/
  DATA bk1cs(2)/ - .353155960776544876E0/
  DATA bk1cs(3)/ - .122611180822657148E0/
  DATA bk1cs(4)/ - .0069757238596398643E0/
  DATA bk1cs(5)/ - .0001730288957513052E0/
  DATA bk1cs(6)/ - .0000024334061415659E0/
  DATA bk1cs(7)/ - .0000000221338763073E0/
  DATA bk1cs(8)/ - .0000000001411488392E0/
  DATA bk1cs(9)/ - .0000000000006666901E0/
  DATA bk1cs(10)/ - .0000000000000024274E0/
  DATA bk1cs(11)/ - .0000000000000000070E0/
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  BESK1
  IF ( first ) THEN
    ntk1 = INITS(bk1cs,11,0.1*R1MACH(3))
    xmin = EXP(MAX(LOG(R1MACH(1)),-LOG(R1MACH(2)))+.01)
    xsml = SQRT(4.0*R1MACH(3))
    xmaxt = -LOG(R1MACH(1))
    xmax = xmaxt - 0.5*xmaxt*LOG(xmaxt)/(xmaxt+0.5)
  ENDIF
  first = .FALSE.
  !
  IF ( X<=0. ) CALL XERMSG('SLATEC','BESK1','X IS ZERO OR NEGATIVE',2,2)
  IF ( X>2.0 ) THEN
    !
    BESK1 = 0.
    IF ( X>xmax ) CALL XERMSG('SLATEC','BESK1','X SO BIG K1 UNDERFLOWS',1,1)
    IF ( X>xmax ) RETURN
    !
    BESK1 = EXP(-X)*BESK1E(X)
    GOTO 99999
  ENDIF
  !
  IF ( X<xmin ) CALL XERMSG('SLATEC','BESK1','X SO SMALL K1 OVERFLOWS',3,2)
  y = 0.
  IF ( X>xsml ) y = X*X
  BESK1 = LOG(0.5*X)*BESI1(X) + (0.75+CSEVL(.5*y-1.,bk1cs,ntk1))/X
  RETURN
  !
  99999 CONTINUE
  END FUNCTION BESK1
