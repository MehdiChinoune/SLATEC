!DECK ACOSH
FUNCTION ACOSH(X)
  IMPLICIT NONE
  REAL ACOSH, aln2, R1MACH, X, xmax
  !***BEGIN PROLOGUE  ACOSH
  !***PURPOSE  Compute the arc hyperbolic cosine.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C4C
  !***TYPE      SINGLE PRECISION (ACOSH-S, DACOSH-D, CACOSH-C)
  !***KEYWORDS  ACOSH, ARC HYPERBOLIC COSINE, ELEMENTARY FUNCTIONS, FNLIB,
  !             INVERSE HYPERBOLIC COSINE
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! ACOSH(X) computes the arc hyperbolic cosine of X.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  R1MACH, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  !***END PROLOGUE  ACOSH
  SAVE aln2, xmax
  DATA aln2/0.69314718055994530942E0/
  DATA xmax/0./
  !***FIRST EXECUTABLE STATEMENT  ACOSH
  IF ( xmax==0. ) xmax = 1.0/SQRT(R1MACH(3))
  !
  IF ( X<1.0 ) CALL XERMSG('SLATEC','ACOSH','X LESS THAN 1',1,2)
  !
  IF ( X<xmax ) ACOSH = LOG(X+SQRT(X*X-1.0))
  IF ( X>=xmax ) ACOSH = aln2 + LOG(X)
  !
END FUNCTION ACOSH
