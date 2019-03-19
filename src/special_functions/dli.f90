!** DLI
REAL(8) FUNCTION DLI(X)
  IMPLICIT NONE
  !>
  !***
  !  Compute the logarithmic integral.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C5
  !***
  ! **Type:**      DOUBLE PRECISION (ALI-S, DLI-D)
  !***
  ! **Keywords:**  FNLIB, LOGARITHMIC INTEGRAL, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! DLI(X) calculates the double precision logarithmic integral
  ! for double precision argument X.
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  DEI, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  
  REAL(8) :: X, DEI
  !* FIRST EXECUTABLE STATEMENT  DLI
  IF ( X<=0.D0 ) CALL XERMSG('SLATEC','DLI',&
    'LOG INTEGRAL UNDEFINED FOR X LE 0',1,2)
  IF ( X==1.D0 ) CALL XERMSG('SLATEC','DLI',&
    'LOG INTEGRAL UNDEFINED FOR X = 0',2,2)
  !
  DLI = DEI(LOG(X))
  !
END FUNCTION DLI
