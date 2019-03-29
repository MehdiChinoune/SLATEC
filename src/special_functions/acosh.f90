!** ACOSH
REAL FUNCTION ACOSH(X)
  IMPLICIT NONE
  !>
  !***
  !  Compute the arc hyperbolic cosine.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C4C
  !***
  ! **Type:**      SINGLE PRECISION (ACOSH-S, DACOSH-D, CACOSH-C)
  !***
  ! **Keywords:**  ACOSH, ARC HYPERBOLIC COSINE, ELEMENTARY FUNCTIONS, FNLIB,
  !             INVERSE HYPERBOLIC COSINE
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! ACOSH(X) computes the arc hyperbolic cosine of X.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  R1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770401  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)

  REAL R1MACH, X
  REAL, PARAMETER :: aln2 = 0.69314718055994530942E0
  REAL :: xmax = 0.
  !* FIRST EXECUTABLE STATEMENT  ACOSH
  IF ( xmax==0. ) xmax = 1.0/SQRT(R1MACH(3))
  !
  IF ( X<1.0 ) CALL XERMSG('SLATEC','ACOSH','X LESS THAN 1',1,2)
  !
  IF ( X<xmax ) ACOSH = LOG(X+SQRT(X*X-1.0))
  IF ( X>=xmax ) ACOSH = aln2 + LOG(X)
  !
END FUNCTION ACOSH
