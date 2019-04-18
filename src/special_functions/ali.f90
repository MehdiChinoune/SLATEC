!** ALI
REAL FUNCTION ALI(X)
  !>
  !***
  !  Compute the logarithmic integral.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C5
  !***
  ! **Type:**      SINGLE PRECISION (ALI-S, DLI-D)
  !***
  ! **Keywords:**  FNLIB, LOGARITHMIC INTEGRAL, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! ALI(X) computes the logarithmic integral; i.e., the
  ! integral from 0.0 to X of (1.0/ln(t))dt.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  EI, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770601  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900326  Removed duplicate information from DESCRIPTION section.
  !           (WRB)
  USE service, ONLY : XERMSG
  REAL X
  !* FIRST EXECUTABLE STATEMENT  ALI
  IF ( X<=0.0 ) CALL XERMSG('SLATEC','ALI', 'LOG INTEGRAL UNDEFINED FOR X LE 0',1,2)
  IF ( X==1.0 ) CALL XERMSG('SLATEC','ALI', 'LOG INTEGRAL UNDEFINED FOR X = 1',2,2)
  !
  ALI = EI(LOG(X))
  !
END FUNCTION ALI
