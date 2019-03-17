!DECK ATANH
REAL FUNCTION ATANH(X)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  ATANH
  !***PURPOSE  Compute the arc hyperbolic tangent.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C4C
  !***TYPE      SINGLE PRECISION (ATANH-S, DATANH-D, CATANH-C)
  !***KEYWORDS  ARC HYPERBOLIC TANGENT, ATANH, ELEMENTARY FUNCTIONS,
  !             FNLIB, INVERSE HYPERBOLIC TANGENT
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! ATANH(X) computes the arc hyperbolic tangent of X.
  !
  ! Series for ATNH       on the interval  0.          to  2.50000D-01
  !                                        with weighted error   6.70E-18
  !                                         log weighted error  17.17
  !                               significant figures required  16.01
  !                                    decimal places required  17.76
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  CSEVL, INITS, R1MACH, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !***END PROLOGUE  ATANH
  REAL atnhcs, CSEVL, dxrel, R1MACH, sqeps, X, y
  INTEGER INITS, nterms
  DIMENSION atnhcs(15)
  LOGICAL first
  SAVE atnhcs, nterms, dxrel, sqeps, first
  DATA atnhcs(1)/.094395102393195492E0/
  DATA atnhcs(2)/.049198437055786159E0/
  DATA atnhcs(3)/.002102593522455432E0/
  DATA atnhcs(4)/.000107355444977611E0/
  DATA atnhcs(5)/.000005978267249293E0/
  DATA atnhcs(6)/.000000350506203088E0/
  DATA atnhcs(7)/.000000021263743437E0/
  DATA atnhcs(8)/.000000001321694535E0/
  DATA atnhcs(9)/.000000000083658755E0/
  DATA atnhcs(10)/.000000000005370503E0/
  DATA atnhcs(11)/.000000000000348665E0/
  DATA atnhcs(12)/.000000000000022845E0/
  DATA atnhcs(13)/.000000000000001508E0/
  DATA atnhcs(14)/.000000000000000100E0/
  DATA atnhcs(15)/.000000000000000006E0/
  DATA first/.TRUE./
  !***FIRST EXECUTABLE STATEMENT  ATANH
  IF ( first ) THEN
    nterms = INITS(atnhcs,15,0.1*R1MACH(3))
    dxrel = SQRT(R1MACH(4))
    sqeps = SQRT(3.0*R1MACH(3))
  ENDIF
  first = .FALSE.
  !
  y = ABS(X)
  IF ( y>=1.0 ) CALL XERMSG('SLATEC','ATANH','ABS(X) GE 1',2,2)
  !
  IF ( 1.0-y<dxrel ) CALL XERMSG('SLATEC','ATANH',&
    'ANSWER LT HALF PRECISION BECAUSE ABS(X) TOO NEAR 1',1,1)
  !
  ATANH = X
  IF ( y>sqeps.AND.y<=0.5 ) ATANH = X*(1.0+CSEVL(8.*X*X-1.,atnhcs,nterms))
  IF ( y>0.5 ) ATANH = 0.5*LOG((1.0+X)/(1.0-X))
  !
END FUNCTION ATANH
