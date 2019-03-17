!DECK BETA
REAL FUNCTION BETA(A,B)
  IMPLICIT NONE
  !***BEGIN PROLOGUE  BETA
  !***PURPOSE  Compute the complete Beta function.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C7B
  !***TYPE      SINGLE PRECISION (BETA-S, DBETA-D, CBETA-C)
  !***KEYWORDS  COMPLETE BETA FUNCTION, FNLIB, SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! BETA computes the complete beta function.
  !
  ! Input Parameters:
  !       A   real and positive
  !       B   real and positive
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  ALBETA, GAMLIM, GAMMA, R1MACH, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !   900727  Added EXTERNAL statement.  (WRB)
  !***END PROLOGUE  BETA
  REAL A, ALBETA, alnsml, B, GAMMA, R1MACH, xmax, xmin
  EXTERNAL GAMMA
  SAVE xmax, alnsml
  DATA xmax, alnsml/0., 0./
  !***FIRST EXECUTABLE STATEMENT  BETA
  IF ( alnsml==0.0 ) THEN
    CALL GAMLIM(xmin,xmax)
    alnsml = LOG(R1MACH(1))
  ENDIF
  !
  IF ( A<=0..OR.B<=0. ) CALL XERMSG('SLATEC','BETA',&
    'BOTH ARGUMENTS MUST BE GT 0',2,2)
  !
  IF ( A+B<xmax ) THEN
    BETA = GAMMA(A)*GAMMA(B)/GAMMA(A+B)
    RETURN
  ENDIF
  !
  BETA = ALBETA(A,B)
  IF ( BETA<alnsml ) CALL XERMSG('SLATEC','BETA',&
    'A AND/OR B SO BIG BETA UNDERFLOWS',1,2)
  !
  BETA = EXP(BETA)
  !
END FUNCTION BETA