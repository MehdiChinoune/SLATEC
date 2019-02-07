!*==DLI.f90  processed by SPAG 6.72Dc at 10:56 on  6 Feb 2019
!DECK DLI
REAL(8) FUNCTION DLI(X)
  IMPLICIT NONE
  !*--DLI5
  !***BEGIN PROLOGUE  DLI
  !***PURPOSE  Compute the logarithmic integral.
  !***LIBRARY   SLATEC (FNLIB)
  !***CATEGORY  C5
  !***TYPE      DOUBLE PRECISION (ALI-S, DLI-D)
  !***KEYWORDS  FNLIB, LOGARITHMIC INTEGRAL, SPECIAL FUNCTIONS
  !***AUTHOR  Fullerton, W., (LANL)
  !***DESCRIPTION
  !
  ! DLI(X) calculates the double precision logarithmic integral
  ! for double precision argument X.
  !
  !***REFERENCES  (NONE)
  !***ROUTINES CALLED  DEI, XERMSG
  !***REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !***END PROLOGUE  DLI
  REAL(8) :: X , DEI
  !***FIRST EXECUTABLE STATEMENT  DLI
  IF ( X<=0.D0 ) CALL XERMSG('SLATEC','DLI',&
    'LOG INTEGRAL UNDEFINED FOR X LE 0',1,2)
  IF ( X==1.D0 ) CALL XERMSG('SLATEC','DLI',&
    'LOG INTEGRAL UNDEFINED FOR X = 0',2,2)
  !
  DLI = DEI(LOG(X))
  !
END FUNCTION DLI
