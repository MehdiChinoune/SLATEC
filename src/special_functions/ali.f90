!DECK ALI
FUNCTION ALI(X)
  IMPLICIT NONE
  REAL ALI, EI, X
  !***BEGIN PROLOGUE  ALI
  !***PURPOSE  Compute the logarithmic integral.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C5
  !***TYPE      SINGLE PRECISION (ALI-S, DLI-D)
  !***KEYWORDS  FNLIB, LOGARITHMIC INTEGRAL, SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! ALI(X) computes the logarithmic integral; i.e., the
  ! integral from 0.0 to X of (1.0/ln(t))dt.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  EI, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !***END PROLOGUE  ALI
  !***FIRST EXECUTABLE STATEMENT  ALI
  IF ( X<=0.0 ) CALL XERMSG('SLATEC','ALI',&
    'LOG INTEGRAL UNDEFINED FOR X LE 0',1,2)
  IF ( X==1.0 ) CALL XERMSG('SLATEC','ALI',&
    'LOG INTEGRAL UNDEFINED FOR X = 1',2,2)
  !
  ALI = EI(LOG(X))
  !
END FUNCTION ALI
