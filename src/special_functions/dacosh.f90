!** DACOSH
REAL(8) FUNCTION DACOSH(X)
  IMPLICIT NONE
  !>
  !***
  !  Compute the arc hyperbolic cosine.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C4C
  !***
  ! **Type:**      DOUBLE PRECISION (ACOSH-S, DACOSH-D, CACOSH-C)
  !***
  ! **Keywords:**  ACOSH, ARC HYPERBOLIC COSINE, ELEMENTARY FUNCTIONS, FNLIB,
  !             INVERSE HYPERBOLIC COSINE
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! DACOSH(X) calculates the double precision arc hyperbolic cosine for
  ! double precision argument X.  The result is returned on the
  ! positive branch.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  D1MACH, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  
  REAL(8) :: X, dln2, xmax, D1MACH
  SAVE dln2, xmax
  DATA dln2/0.69314718055994530941723212145818D0/
  DATA xmax/0.D0/
  !* FIRST EXECUTABLE STATEMENT  DACOSH
  IF ( xmax==0.D0 ) xmax = 1.0D0/SQRT(D1MACH(3))
  !
  IF ( X<1.D0 ) CALL XERMSG('SLATEC','DACOSH','X LESS THAN 1',1,2)
  !
  IF ( X<xmax ) DACOSH = LOG(X+SQRT(X*X-1.0D0))
  IF ( X>=xmax ) DACOSH = dln2 + LOG(X)
  !
END FUNCTION DACOSH
